"""
utilities for patching AoU VDS, add new filters, update reference data with ploidy information

The 'entry point' is ``patch_vds``.

The other functions can be used with a subset of the data for testing. For example:

    vds = hl.vds.read_vds(ECHO_PATH)

    vd = vds.variant_data
    rd = vds.reference_data

    # to patch the reference data
    rd = hl.filter_intervals(REF_TESTING_INTERVALS)
    ploidy = import_ploidy(PLOIDY_AVROS)

    rd = patch_reference_data(rd, ploidy)

    # run checks to for reference ploidy/genotypes here...

    # to patch the variant data
    vd = hl.filter_intervals(VAR_TESTING_INTERVALS)

    site_path = os.path.join(TMP_DIR, 'site_filters.ht')
    site = import_site_filters(SITE_AVROS, site_path, VAR_TESTING_INTERVALS)

    vets_path = os.path.join(TMP_DIR, 'vets.ht')
    vets = import_vets(vets_AVROS, vets_path, VAR_TESTING_INTERVALS)

    vd = patch_variant_data(vd, site=site, vets=vets)

    # run checks for variant/filtering data here...
"""
import os
import json
import gzip

from collections import namedtuple, defaultdict, abc

import hail as hl

from avro.datafile import DataFileReader
from avro.io import DatumReader
from hail.utils.java import info
from hailtop.fs.router_fs import RouterFS


def patch_vds(
    *,
    vds_path: str,
    site_filtering_data: abc.Sequence[str],
    vets_filtering_data: abc.Sequence[str],
    ploidy_data: str | abc.Sequence[str],
    output_path: str,
    tmp_dir: str,
    overwrite: bool = False,
) -> hl.vds.VariantDataset:
    """
    Parameters
    ----------
    vds_path : str
        Path to the current vds
    site_filtering_data : list[str]
        Paths to site filtering files.
    vets_filtering_data : list[str]
        Paths to VETS/VQSR filtering files.
    ploidy_data : str | list[str]
        Path(s) to ploidy data file(s).
    output_path : str
        Path to the new vds
    tmp_dir : str
        Temporary directory
    overwrite : bool
        Overwrite ``output_path``?
    """
    vds_intervals = extract_vds_intervals(vds_path)
    vds = hl.vds.read_vds(vds_path)

    site_path = os.path.join(tmp_dir, "site_filters.ht")
    vets_path = os.path.join(tmp_dir, "vets.ht")

    site = import_site_filters(site_filtering_data, site_path, vds_intervals)
    vets = import_vets(vets_filtering_data, vets_path, vds_intervals)

    if isinstance(ploidy_data, str):
        ploidy_data = (ploidy_data,)
    ploidy = import_ploidy(*ploidy_data)

    variant_data = patch_variant_data(vds.variant_data, site, vets)
    reference_data = patch_reference_data(vds.reference_data, ploidy)


    return hl.vds.VariantDataset(
        reference_data=reference_data, variant_data=variant_data
    ).checkpoint(output_path, overwrite=overwrite)


def extract_vds_intervals(path: str) -> list[hl.Interval]:
    """Extracts the partition bounds from a VDS path"""
    fs = hl.current_backend().fs
    md_file = os.path.join(path, "reference_data", "rows", "rows", "metadata.json.gz")
    with fs.open(md_file, "rb") as md_gz_stream:
        with gzip.open(md_gz_stream) as md_stream:
            metadata = json.load(md_stream)
    interval_list_type = hl.tarray(
        hl.tinterval(hl.tstruct(locus=hl.tlocus(reference_genome="GRCh38")))
    )
    return interval_list_type._convert_from_json(metadata["_jRangeBounds"])


_GRCH38 = None


def translate_locus(location):
    """Translates an int64-encoded locus into a locus object."""
    global _GRCH38
    if _GRCH38 is None:
        _GRCH38 = hl.get_reference("GRCh38")
    factor = 1_000_000_000_000
    chrom = hl.literal(_GRCH38.contigs[:26])[hl.int32(location / factor) - 1]
    pos = hl.int32(location % factor)
    return hl.locus(chrom, pos, reference_genome=_GRCH38)


def import_site_filters(
    avros: abc.Sequence[str], site_path: str, intervals: list[hl.Interval], force=False
) -> hl.Table:
    """
    Parameters
    ----------
    avros : Sequence[str]
        List of paths for raw site filtering data
    site_path : str
        Path to site filters table where, if a hail table exists, it will be read, unless ``force`` is true
    intervals : list[Interval]
        a list of intervals to read the table
    force : bool
        always import the filtering data?
    """
    fs = hl.current_backend().fs

    if force or not fs.exists(os.path.join(site_path, "_SUCCESS")):
        info("Importing and writing site filters to temporary storage")
        site = hl.import_avro(avros)
        site = site.transmute(
            locus=translate_locus(site.location),
            filters=hl.set(site.filters.split(",")),
        )
        site = site.key_by("locus")
        site.write(site_path, overwrite=True)

    return hl.read_table(site_path, _intervals=intervals)


def import_vets(
    avros: abc.Sequence[str], vets_path: str, intervals: list[hl.Interval], force=False
) -> hl.Table:
    """
    Parameters
    ----------
    avros : Sequence[str]
        List of paths for vets/vets data
    vets_path : str
        Path to variant filters table where, if a hail table exists, it will be read, unless ``force`` is true
    intervals : list[Interval]
        a list of intervals to read the table
    force : bool
        always import the filtering data?
    """
    fs = hl.current_backend().fs

    if force or not fs.exists(os.path.join(vets_path, "_SUCCESS")):
        info("Importing and writing vets filter data to temporary storage")
        vets = hl.import_avro(avros)
        vets = vets.transmute(locus=translate_locus(vets.location))
        vets = vets.key_by("locus")
        vets.write(vets_path, overwrite=True)

    return hl.read_table(vets_path, _intervals=intervals)


def import_ploidy(*avros) -> dict[str, dict[str, int]]:
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
    return {
        contig: ploidy_table[key]
        for contig, key in zip(hg38.contigs, sorted(ploidy_table))
    }


def patch_variant_data(vd, site, vets) -> hl.MatrixTable:
    """
    Parameters
    ----------
    vd : MatrixTable
        vds variant data
    site : Table
        site filtering table
    vets : Table
        vets filtering table
    """
    vd = vd.annotate_rows(
        filters=hl.coalesce(site[vd.locus].filters, hl.empty_set(hl.tstr))
    )

    # vets ref/alt come in normalized individually, so need to renormalize to the dataset ref allele
    vd = vd.annotate_rows(
        as_vets=hl.dict(
            vets.index(vd.locus, all_matches=True).map(
                lambda record: (
                    record.alt + vd.alleles[0][hl.len(record.ref) :],
                    record.drop("ref", "alt"),
                )
            )
        )
    )

    is_snp = vd.alleles[1:].map(lambda alt: hl.is_snp(vd.alleles[0], alt))
    vd = vd.annotate_rows(
        allele_NO=vd.alleles[1:].map(
            lambda allele: hl.coalesce(vd.as_vets.get(allele).yng_status == "N", False)
        ),
        allele_YES=vd.alleles[1:].map(
            lambda allele: hl.coalesce(vd.as_vets.get(allele).yng_status == "Y", True)
        ),
        allele_is_snp=is_snp,
        allele_OK=hl._zip_func(
            is_snp,
            vd.alleles[1:],
            f=lambda is_snp, alt: hl.coalesce(
                vd.as_vets.get(alt).calibration_sensitivity
                <= hl.if_else(
                    is_snp,
                    vd.truth_sensitivity_snp_threshold,
                    vd.truth_sensitivity_indel_threshold,
                ),
                True,
            ),
        ),
        as_vets=vd.as_vets.map_values(lambda value: value.drop("yng_status")),
    )
    lgt = vd.LGT
    la = vd.LA
    allele_NO = vd.allele_NO
    allele_YES = vd.allele_YES
    allele_OK = vd.allele_OK
    allele_is_snp = vd.allele_is_snp
    ft = (
        hl.range(lgt.ploidy)
        .map(lambda idx: la[lgt[idx]])
        .filter(lambda x: x != 0)
        .fold(
            lambda acc, called_idx: hl.struct(
                any_no=acc.any_no | allele_NO[called_idx - 1],
                any_yes=acc.any_yes | allele_YES[called_idx - 1],
                any_snp=acc.any_snp | allele_is_snp[called_idx - 1],
                any_indel=acc.any_indel | ~allele_is_snp[called_idx - 1],
                any_snp_ok=acc.any_snp_ok
                | (allele_is_snp[called_idx - 1] & allele_OK[called_idx - 1]),
                any_indel_ok=acc.any_indel_ok
                | (~allele_is_snp[called_idx - 1] & allele_OK[called_idx - 1]),
            ),
            hl.struct(
                any_no=False,
                any_yes=False,
                any_snp=False,
                any_indel=False,
                any_snp_ok=False,
                any_indel_ok=False,
            ),
        )
    )

    vd = vd.annotate_entries(
        FT=~ft.any_no
        & (
            ft.any_yes
            | ((~ft.any_snp | ft.any_snp_ok) & (~ft.any_indel | ft.any_indel_ok))
        )
    )

    vd = vd.drop("allele_NO", "allele_YES", "allele_is_snp", "allele_OK")
    return vd


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
    rd = rd.annotate_globals(ploidy_data=hl.literal(ploidy))
    rd = rd.annotate_rows(autosome_or_par=rd.locus.in_autosome_or_par())
    rd = rd.annotate_entries(
        GT=hl.if_else(
            rd.autosome_or_par,
            hl.call(0, 0),
            hl.rbind(
                rd.ploidy_data[rd.locus.contig].get(rd.s, 2),
                lambda ploidy: hl.switch(ploidy)
                .when(1, hl.call(0))
                .when(2, hl.call(0, 0))
                .or_error(
                    "expected 1 or 2 for ploidy information, found: " + hl.str(ploidy)
                ),
            ),
        )
    )

    return rd.drop("ploidy_data", "autosome_or_par")
