import argparse
import hail as hl
from typing import Union
from datetime import datetime

import create_vat_inputs

###
# VAT preparation:
# filtering:
# hard filter bad sites: hard_filter_non_passing_sites()
# hard filter based on FT flag:
# get GT, replace_lgt_with_gt(), swap to no-calls: failing_gts_to_no_call()
# track how many sites have more than 50 alt alleles TODO: we currently aren't tracking this, are we?
# drop 50+ alternate alleles
# calculate the AC, AN, AF, SC, for the full population and for the subpopulations
# split multi allelic sites to be one variant per row (note that we want to do this after the above calculations)
# drop all spanning deletions
### The rest of these things we're gonna do with bcftools?
# make a tsv for custom annotations (do we want to do this with bcftools or just like grep in a WDL?)
# drop sites where AC=0??!?!?!??!! (We probably want to do this _after_ we join this back with the original VDS?)




def hard_filter_non_passing_sites(vds):
    """
    Hard filter out non-passing sites
    note: no AC/AN and AF for filtered out positions
    """
    vd = vds.variant_data
    filtered_vd = vd.filter_rows(hl.len(vd.filters) == 0)
    # now we apply it back to the vds
    filtered_vds = hl.vds.VariantDataset(vds.reference_data, filtered_vd)
    return filtered_vds


def replace_lgt_with_gt(vds):
    """
    Replace LGT with GT for easier calculations later.
    """
    filtered_vd = vds.variant_data
    filtered_vd = filtered_vd.annotate_entries(GT=hl.vds.lgt_to_gt(filtered_vd.LGT, filtered_vd.LA))
    return hl.vds.VariantDataset(vds.reference_data, filtered_vd)


def failing_gts_to_no_call(vds):
    """
     Respect the FT flag by setting all failing GTs to a no call

     Logic for assigning non passing GTs as no-calls to ensure that AC,AN and AF respect the FT flag.

     FT is True => GT keeps its current value
     FT is False => GT assigned no-call
     FT is missing => GT keeps its current value
    """
    gt_vds = replace_lgt_with_gt(vds)
    vd = gt_vds.variant_data
    filtered_vd = vd.annotate_entries(GT=hl.or_missing(hl.coalesce(vd.FT, True), vd.GT))
    # GT will not line up with the LGT but this version of the VDS data will be sliced and not fully moved anywhere
    return hl.vds.VariantDataset(vds.reference_data, filtered_vd)


def remove_too_many_alt_allele_sites(vds):
    """
     Remove sites with more than 50 alternate alleles
    """
    vd = vds.variant_data
    vd_50_aa_cutoff = vd.filter_rows(hl.len(vd.alleles) <= 50)

    return hl.vds.VariantDataset(vds.reference_data, vd_50_aa_cutoff)



def annotate_adj(
        mt: hl.MatrixTable,
        adj_gq: int = 30,
        adj_ab: float = 0.2,
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria (assumes diploid).
    Defaults similar to gnomAD values, but GQ >= 20 changed to GQ >= 30 to make up for lack of DP filter.
    """
    if "LGT" in mt.entry and "LAD" in mt.entry:
        gt_expr = mt.LGT
        ad_expr = mt.LAD
    else:
        assert "GT" in mt.entry and "AD" in mt.entry
        gt_expr = mt.GT
        ad_expr = mt.AD
    return mt.annotate_entries(
        adj=get_adj_expr(
            gt_expr, mt.GQ, ad_expr, adj_gq, adj_ab
        )
    )
def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 30,
    adj_ab: float = 0.2,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.
    Defaults similar to gnomAD values, but GQ >= 20 changed to GQ >= 30 to make up for lack of DP filter.
    """
    total_ad = hl.sum(ad_expr)
    return (
        (gq_expr >= adj_gq)
        & (
            hl.case()
            .when(~gt_expr.is_het(), True)
            .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / total_ad >= adj_ab)
            .default(
                (ad_expr[gt_expr[0]] / total_ad >= adj_ab)
                & (ad_expr[gt_expr[1]] / total_ad >= adj_ab)
            )
        )
    )
def matrix_table_ac_an_af(mt, ancestry_file):
    """
    Create a DENSE MATRIX TABLE to calculate AC, AN, AF and Sample Count
    TODO: test sample_count
    """
    sample_id_to_sub_population = create_vat_inputs.parse_ancestry_file(ancestry_file)
    mt = mt.annotate_cols(pop=hl.literal(sample_id_to_sub_population)[mt.s])
    mt = annotate_adj(mt)
    ac_an_af_mt = mt.select_rows(
        ac_an_af=hl.agg.call_stats(mt.GT, mt.alleles),
        call_stats_by_pop=hl.agg.group_by(mt.pop, hl.agg.call_stats(mt.GT, mt.alleles)),
        ac_an_af_adj=hl.agg.filter(mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)),
        call_stats_by_pop_adj=hl.agg.filter(mt.adj, hl.agg.group_by(mt.pop, hl.agg.call_stats(mt.GT, mt.alleles))),
    )

    return hl.methods.split_multi(ac_an_af_mt, left_aligned=True) # split each alternate allele onto it's own row. This will also remove all spanning delstions for us


def vds_ac_an_af(mt, vds):
    """
    Join the dense matrix table input (just AC, AN and AF) back to the original VDS. Note any hard filtering
    """
    qc_data = mt.rows()
    vd = vds.variant_data
    vd = vd.annotate_rows(ac_an_af=qc_data[vd.row_key])
    return hl.vds.VariantDataset(vds.reference_data, vd)


def write_sites_only_vcf(ac_an_af_split, sites_only_vcf_path):
    # CHROM	POS	REF	ALT	AC	AN	AF	Hom	AC_afr	AN_afr	AF_afr	Hom_afr	AC_amr	AN_amr	AF_amr	Hom_amr	AC_eas	AN_eas	AF_eas	AC_Hom_eas	AC_Het_eas	AC_Hemi_eas	AC_eur	AN_eur	AF_eur	AC_Hom_eur	AC_Het_eur	AC_Hemi_eur	AC_mid	AN_mid	AF_mid	AC_Hom_mid	AC_Het_mid	AC_Hemi_mid	AC_oth	AN_oth	AF_oth	AC_Hom_oth	AC_Het_oth	AC_Hemi_oth	AC_sas	AN_sas	AF_sas	AC_Hom_sas	AC_Het_sas	AC_Hemi_sas

    pop_info_fields = {}
    for pop in ['afr', 'amr', 'eas', 'eur', 'mid', 'oth', 'sas']: ## TODO double check that this is all of them
        pop_info_fields[f'AC_{pop}'] = ac_an_af_split.call_stats_by_pop.get(pop).AC[ac_an_af_split.row.a_index]
        pop_info_fields[f'AN_{pop}'] = ac_an_af_split.call_stats_by_pop.get(pop).AN
        pop_info_fields[f'AF_{pop}'] = ac_an_af_split.call_stats_by_pop.get(pop).AF[ac_an_af_split.row.a_index]
        pop_info_fields[f'Hom_{pop}'] = ac_an_af_split.call_stats_by_pop.get(pop).homozygote_count[ac_an_af_split.row.a_index]


    ac_an_af_rows = ac_an_af_split.annotate(
        info = hl.struct(
            AC=ac_an_af_split.row.ac_an_af.AC[ac_an_af_split.row.a_index],
            AN=ac_an_af_split.row.ac_an_af.AN,
            AF=ac_an_af_split.row.ac_an_af.AF[ac_an_af_split.row.a_index],
            Hom=ac_an_af_split.row.ac_an_af.homozygote_count[ac_an_af_split.row.a_index],
            **pop_info_fields
        )
    )

    # note that SC = AC - homozygote_count

    ht = ac_an_af_rows.filter(ac_an_af_rows.alleles[1] != "*") # remove spanning deletions
    # create a filtered sites only VCF
    hl.export_vcf(ht, sites_only_vcf_path)


def add_variant_tracking_info(mt, sites_only_vcf_path):
    # only need the table of row fields and leaves this as the only field
    var_ids_path = sites_only_vcf_path.replace(r".sites-only.vcf", ".var_ids.tsv.bgz")
    t = mt.rows()
    t.select(var_origin_id=hl.format('%s-%s-%s-%s', t.locus.contig, t.locus.position, t.alleles[0], t.alleles[1])).export(var_ids_path, parallel='separate_header')


def main(vds, ancestry_file_location, sites_only_vcf_path, dry_run_n_parts=1):
    n_parts = vds.variant_data.n_partitions()
    if dry_run_n_parts:
        n_rounds = 1
        parts_per_round = dry_run_n_parts
        ht_paths = [sites_only_vcf_path.replace(r".sites-only.vcf.bgz", f'_dryrun.ht')]
        sites_only_vcf_path = sites_only_vcf_path.replace(r".vcf.bgz", f'_dryrun.vcf.bgz')
    else:
        n_rounds = 5
        # Add in 'n_rounds - 1' to include all of the partitions in the set of groups, otherwise we would omit the final
        # n_parts % n_rounds partitions.
        parts_per_round = (n_parts + n_rounds - 1) // n_rounds
        ht_paths = [sites_only_vcf_path.replace(r".sites-only.vcf.bgz", f'_{i}.ht') for i in range(n_rounds)]
    for i in range(n_rounds):
        part_range = range(i*parts_per_round, min((i+1)*parts_per_round, n_parts))
        vds_part = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(part_range),
            vds.variant_data._filter_partitions(part_range),
        )

        transforms = [
            remove_too_many_alt_allele_sites,
            hard_filter_non_passing_sites,
            failing_gts_to_no_call,
        ]
        transformed_vds=vds_part
        for transform in transforms:
            transformed_vds = transform(transformed_vds)

        print(f"densifying dense matrix table index {i}")
        mt = hl.vds.to_dense_mt(transformed_vds)

        with open(ancestry_file_location, 'r') as ancestry_file:
            mt = matrix_table_ac_an_af(mt, ancestry_file) # this adds subpopulation information and splits our multi-allelic rows

        ht = mt.rows()
        ht = ht.select('call_stats_by_pop', 'a_index', 'ac_an_af', 'ac_an_af_adj', 'call_stats_by_pop_adj')
        ht.write(ht_paths[i])

        # potentially in the future: merge AC, AN, AF back to the original VDS with: vds = vds_ac_an_af(mt, vds)

        # for debugging information -- remove for now to get us through Echo
        # add_variant_tracking_info(mt, sites_only_vcf_path)

    # create a sites only VCF (that is hard filtered!) and that can be made into a custom annotations TSV for Nirvana to use with AC, AN, AF, SC for all subpopulations and populations
    ht_list = [hl.read_table(ht_path).naive_coalesce(5000) for ht_path in ht_paths] # repartition each table to 5k partitions before we union them
    ht_all = ht_list[0].union(*ht_list[1:])
    write_sites_only_vcf(ht_all, sites_only_vcf_path)


def annotate_entry_filter_flag(mt):
    """
    Annotate FT flag for entries--this is needed so that tools like bcftools can handle the output vcf
    """
    fail_case = 'FAIL'
    return mt.annotate_entries(FT=hl.if_else(mt.FT, 'PASS', fail_case))


def write_tie_out_vcf(vds, vcf_output_path):
    """
    This is for writing tieout VCFs for toy-sized data only. Do not use for AoU-sized data as it would write a giant
    MatrixTable!
    """
    mt = hl.vds.to_dense_mt(vds)
    mt = annotate_entry_filter_flag(mt)

    # Write VCF
    hl.export_vcf(mt, vcf_output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Create VAT inputs TSV')
    parser.add_argument('--ancestry_input_path', type=str, help='Input ancestry file path', required=True)
    parser.add_argument('--vds_input_path', type=str, help='Input VDS path')
    parser.add_argument('--sites_only_output_path', type=str, help='Output sites-only VCF path'),
    parser.add_argument('--temp_path', type=str, help='Path to temporary directory', required=True)

    args = parser.parse_args()

    try:
        temp_path = args.temp_path if not args.temp_path.endswith('/') else args.temp_path[:-1]
        time_stamp = datetime.today().strftime('%Y-%m-%d_%H-%M-%S')
        hl.init(tmp_dir=f'{temp_path}/hail_tmp_create_vat_inputs_{time_stamp}')

        vds = hl.vds.read_vds(args.vds_input_path)
        local_ancestry_file = create_vat_inputs.download_ancestry_file(args.ancestry_input_path)

        main(vds, local_ancestry_file, args.sites_only_output_path)
    except Exception:
      hl.copy_log(args.sites_only_output_path.replace(r".sites-only.vcf", ".log"))
      raise