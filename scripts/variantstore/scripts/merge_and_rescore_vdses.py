import argparse
import os
from collections import abc
from typing import Sequence

import hail as hl
from hail.utils.java import info


def merge_vdses(vds1, vds2, output_path, hail_temp_path):
    combiner = hl.vds.new_combiner(
        output_path=output_path,
        temp_path=hail_temp_path,
        vds_paths=[vds1, vds2],
        use_genome_default_intervals=True)

    combiner.run()


def translate_locus(location):
    """Translates an int64-encoded locus into a locus object."""
    global _GRCH38
    if _GRCH38 is None:
        _GRCH38 = hl.get_reference("GRCh38")
    factor = 1_000_000_000_000
    chrom = hl.literal(_GRCH38.contigs[:26])[hl.int32(location / factor) - 1]
    pos = hl.int32(location % factor)
    return hl.locus(chrom, pos, reference_genome=_GRCH38)


def find_site_and_vets_filter_data(avro_path: str) -> (Sequence[str], Sequence[str]):
    from hail_gvs_util import gcs_generate_avro_args, gcs_pattern
    from google.cloud import storage

    avro_bucket_name, avro_object_prefix = gcs_pattern.match(avro_path).groups()
    avro_bucket = storage.Client().get_bucket(avro_bucket_name)
    ret = [gcs_generate_avro_args(avro_bucket, avro_object_prefix, key)
           for key in ["site_filtering_data", "vets_filtering_data"]]
    return ret


def import_site_filters(
        avros: abc.Sequence[str],
        site_path: str
) -> hl.Table:
    """
Parameters
----------
avros : Sequence[str]
    List of paths for raw site filtering data
site_path : str
    Path to where site filters table will be written
"""
    info("Importing and writing site filters to temporary storage")
    site = hl.import_avro(avros)
    site = site.transmute(
        locus=translate_locus(site.location),
        filters=hl.set(site.filters.split(",")),
    )
    site = site.key_by("locus")
    site.write(site_path, overwrite=True)
    return hl.read_table(site_path)


def import_vets_filters(
        avros: abc.Sequence[str], vets_path: str
) -> hl.Table:
    """
    Parameters
    ----------
    avros : Sequence[str]
        List of paths for vets allele/variant filtering data
    vets_path : str
        Path to where variant filters table will be written

    """

    info("Importing and writing vets filter data to temporary storage")
    vets = hl.import_avro(avros)
    vets = vets.transmute(locus=translate_locus(vets.location))
    vets = vets.key_by("locus")
    vets.write(vets_path, overwrite=True)
    return hl.read_table(vets_path)


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Given two Hail VDSes and a path to exported Avro files containing the latest filter data, merge the two VDSes to produce a single output VDSes with updated filter information.')
    parser.add_argument('--input-vds-first-path', type=str,
                        help='Path to first input VDS', required=True)
    parser.add_argument('--input-vds-second-path', type=str,
                        help='Path to second input VDS', required=True)
    parser.add_argument('--avro-path', type=str,
                        help='Path at which exported GVS Avro files are found, including the filter data to apply in the output merged VDS.',
                        required=True)
    parser.add_argument('--output-vds-path', type=str,
                        help='Path to write output VDS', required=True)

    parser.add_argument('--temp-path', type=str, help='Path to temporary directory', required=True)

    args = parser.parse_args()
    tmp_dir = f'{args.temp_path}/hail_tmp_general'
    hl.init(
        idempotent=True,
        tmp_dir=tmp_dir
    )
    hl.default_reference('GRCh38')

    site_path = os.path.join(tmp_dir, "site_filters.ht")
    vets_path = os.path.join(tmp_dir, "vets_filters.ht")

    site_filtering_data, vets_filtering_data = find_site_and_vets_filter_data(args.avro_path)

    site = import_site_filters(site_filtering_data, site_path)
    vets = import_vets_filters(vets_filtering_data, vets_path)

    merge_vdses(args.input_vds_first_path,
                args.input_vds_second_path,
                args.output_vds_path,
                tmp_dir)

    merged_vds = hl.vds.read_vds(args.output_vds_path)

    patch_variant_data(merged_vds, site, vets)
