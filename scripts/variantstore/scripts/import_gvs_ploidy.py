import os
import json
import gzip

from collections import namedtuple, defaultdict, abc

import hail as hl

from avro.datafile import DataFileReader
from avro.io import DatumReader
from hail.utils.java import info

def import_ploidy(*avros) -> dict[str, hl.Struct]:
    """
    Parameters
    ----------
    avros :
        Path(s) of ploidy data
    """
    PloidyRecord = namedtuple("PloidyRecord", "location sample_name ploidy")

    # the implementation of GCS for Hadoop doesn't allow seeking to the end of a file
    # so I'm monkey patching DataFileReader
    def patched_determine_file_length(self) -> int:
        remember_pos = self.reader.tell()
        self.reader.seek(-1, 2)
        file_length = self.reader.tell() + 1
        self.reader.seek(remember_pos)
        return file_length

    original_determine_file_length = DataFileReader.determine_file_length
    DataFileReader.determine_file_length = patched_determine_file_length

    fs = hl.current_backend().fs
    ploidy_table = defaultdict(dict)
    for file in avros:
        with fs.open(file, "rb") as data:
            for record in DataFileReader(data, DatumReader()):
                location, sample_name, ploidy = PloidyRecord(**record)
                if sample_name in ploidy_table[location]:
                    raise ValueError(
                        f"duplicate key `{sample_name}` for location {location}"
                    )
                ploidy_table[location][sample_name] = ploidy

    # undo our monkey patch
    DataFileReader.determine_file_length = original_determine_file_length

    hg38 = hl.get_reference("GRCh38")
    xy_contigs = set(hg38.x_contigs + hg38.y_contigs)
    ploidy_table = {
        contig: ploidy_table[key]
        for contig, key in zip(hg38.contigs, sorted(ploidy_table))
        if contig in xy_contigs
    }
    print(f"ploidy table keys are {', '.join(ploidy_table.keys()}")
    x_table = ploidy_table["chrX"]
    y_table = ploidy_table["chrY"]
    assert set(x_table) == set(y_table)

    return {
        sample_name: hl.Struct(
            x_ploidy=x_table[sample_name], y_ploidy=y_table[sample_name]
        )
        for sample_name in x_table
    }


def patch_reference_data(rd, ploidy) -> hl.MatrixTable:
    """
    Parameters
    ----------
    rd : MatrixTable
        vds reference data
    ploidy : dict[str, dict[str, int]]
        table of ploidy information. Keys of outer dict are contigs. Keys of inner dict are sample names.
        Values of inner dict are the ploidy to use for the reference genotype in nonpar regions.
    """
    rd = rd.annotate_cols(ploidy_data=hl.literal(ploidy)[rd.s])
    rd = rd.annotate_rows(autosome_or_par=rd.locus.in_autosome_or_par(), is_y=rd.locus.contig == 'chrY')
    rd = rd.annotate_entries(
        GT=hl.if_else(
            rd.autosome_or_par,
            hl.call(0, 0),
            hl.rbind(
                hl.if_else(rd.is_y, rd.ploidy_data.y_ploidy, rd.ploidy_data.x_ploidy),
                lambda ploidy: hl.switch(ploidy)
                .when(1, hl.call(0))
                .when(2, hl.call(0, 0))
                .or_error(
                    "expected 1 or 2 for ploidy information, found: " + hl.str(ploidy)
                ),
            ),
        )
    )

    return rd.drop("ploidy_data", "autosome_or_par", "is_y")
