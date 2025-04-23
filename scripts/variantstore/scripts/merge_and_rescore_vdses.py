import argparse
import os
from collections import abc
from typing import Sequence

import hail as hl
from hail.utils.java import info


_GRCH38 = None


def merge_vdses(general_temp_path: str, output_merged_vds_path: str, *input_vds):
    combiner = hl.vds.new_combiner(
        output_path=output_merged_vds_path,
        temp_path=general_temp_path,
        vds_paths=list(input_vds),
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
    info("Importing and writing site filters to temporary storage at {}".format(site_path))
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

    info("Importing and writing vets filter data to temporary storage at {}".format(vets_path))
    vets = hl.import_avro(avros)
    vets = vets.transmute(locus=translate_locus(vets.location))
    vets = vets.key_by("locus")
    vets.write(vets_path, overwrite=True)
    return hl.read_table(vets_path)


def patch_variant_data(vd: hl.MatrixTable, site_filters: hl.Table, vets_filters: hl.Table) -> hl.MatrixTable:
    """
    Parameters
    ----------
    vd : MatrixTable
        vds variant data
    site_filters : Table
        site filtering table
    vets_filters : Table
        vets filtering table
    """
    vd = vd.annotate_rows(
        filters=hl.coalesce(site_filters[vd.locus].filters, hl.empty_set(hl.tstr))
    )

    # vets ref/alt come in normalized individually, so need to renormalize to the dataset ref allele
    vd = vd.annotate_rows(
        as_vets=hl.dict(
            vets_filters.index(vd.locus, all_matches=True).map(
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
                all_snps_ok=acc.all_snps_ok
                            & (~allele_is_snp[called_idx - 1] | allele_OK[called_idx - 1]),
                all_indels_ok=acc.all_indels_ok
                              & (allele_is_snp[called_idx - 1] | allele_OK[called_idx - 1]),
            ),
            hl.struct(
                any_no=False,
                any_yes=False,
                all_snps_ok=True,
                all_indels_ok=True,
            ),
        )
    )

    vd = vd.annotate_entries(
        FT=~ft.any_no
           & (
                   ft.any_yes | (ft.all_snps_ok & ft.all_indels_ok)
           )
    )

    vd = vd.drop("allele_NO", "allele_YES", "allele_is_snp", "allele_OK")
    return vd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Given two Hail VDSes and a path to exported Avro files containing the latest filter data, merge the two VDSes to produce a single output VDS with updated filter.')

    parser.add_argument('--input-echo-vds', help="Echo VDS with Echo filter, will not be overwritten", required=True)

    parser.add_argument('--input-unmerged-foxtrot-vds', help="Foxtrot VDS with Foxtrot filter, contains only samples new for Foxtrot, will not be overwritten", required=True)

    parser.add_argument('--input-foxtrot-avro-path', type=str,
                        help='Path at which exported Foxtrot GVS Avro files are found, including the filter data to apply in the output merged VDS.',
                        required=True)

    parser.add_argument('--output-vds-path', type=str,
                        help='Path to write output VDS', required=True)

    parser.add_argument('--truth-sensitivity-snp-threshold', type=float,
                        help='Sensitivity threshold for SNP filtering.',
                        default=0.997)

    parser.add_argument('--truth-sensitivity-indel-threshold', type=float,
                        help='Sensitivity threshold for indel filtering.',
                        default=0.99)

    parser.add_argument('--samples-to-remove-path', type=str,
                        help="File with ids of samples to remove, one sample id per line, header should be 'research_id'.",
                        required=False)

    parser.add_argument('--temp-path', type=str, help='Path to temporary directory', required=True)

    args = parser.parse_args()
    tmp_dir = f'{args.temp_path}/hail_tmp_general'
    tmp_merged_vds_path = f'{args.temp_path}/hail_tmp_merge/tmp_merged.vds'

    hl.init(
        idempotent=True,
        tmp_dir=tmp_dir
    )
    hl.default_reference('GRCh38')
    samples_to_remove_table = None

    # Do this first to fail fast if there's a problem.
    if args.samples_to_remove_path:
        samples_to_remove_table = hl.import_table(args.samples_to_remove_path, delimiter=',')
        samples_to_remove_table = samples_to_remove_table.key_by(s=samples_to_remove_table.research_id)

    site_path = os.path.join(tmp_dir, "site_filters.ht")
    vets_path = os.path.join(tmp_dir, "vets_filters.ht")

    site_filtering_data, vets_filtering_data = find_site_and_vets_filter_data(args.input_foxtrot_avro_path)

    site = import_site_filters(site_filtering_data, site_path)
    vets = import_vets_filters(vets_filtering_data, vets_path)

    # First merge the Echo and new-to-Foxtrot VDSes, then rescore. Although the new-to-Foxtrot VDS already has correct
    # filter data, we don't get correct AC/AN/AF results during our VDS tieout tests if we score just the Foxtrot-only
    # VDS and then merge.
    merge_vdses(tmp_dir,
                tmp_merged_vds_path,
                args.input_echo_vds, args.input_unmerged_foxtrot_vds)

    tmp_merged_vds = hl.vds.read_vds(tmp_merged_vds_path)
    # Drop any samples that need dropping before patching.
    if samples_to_remove_table:
        tmp_merged_vds = hl.vds.filter_samples(tmp_merged_vds, samples_to_remove_table, keep=False,
                                               remove_dead_alleles=True)

    # These globals seem to get dropped after the merge. Add them back as they are required for rescoring.
    vd = tmp_merged_vds.variant_data
    vd = vd.annotate_globals(truth_sensitivity_snp_threshold=args.truth_sensitivity_snp_threshold,
                             truth_sensitivity_indel_threshold=args.truth_sensitivity_indel_threshold)

    rescored_vd = patch_variant_data(vd, site, vets)
    merged_and_rescored_vds = hl.vds.VariantDataset(tmp_merged_vds.reference_data, rescored_vd)
    merged_and_rescored_vds.write(args.output_vds_path)
