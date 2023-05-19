import argparse
import hail as hl

import create_vat_inputs


###
# VAT preparation:
# filtering:
# hard filter bad sites: hard_filter_non_passing_sites()
# hard filter based on FT flag:
# get GT, ~replace_lgt_with_gt()~, swap to no-calls: failing_gts_to_no_call(), (no longer need to replace_lgt_with_gt now that lgt is included in the VDS)
# turn the GQ0s into no calls so that ANs are correct -- hmmm maybe we should check with Lee about this because we do have GQ0s in the VDS
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
    vd = vds.variant_data
    filtered_vd = vd.annotate_entries(GT=hl.or_missing(hl.coalesce(vd.FT, True), vd.GT))
    # GT will not line up with the LGT but this version of the VDS data will be sliced and not fully moved anywhere
    return hl.vds.VariantDataset(vds.reference_data, filtered_vd)


def remove_too_many_alt_allele_sites(vds):
    """
     Remove sites with more than 50 alternate alleles (and print how many)
    """
    vd = vds.variant_data
    print(vd.aggregate_rows(hl.agg.count_where(hl.len(vd.alleles) > 50)))
    vd_50_aa_cutoff = vd.filter_rows(hl.len(vd.alleles) <= 50)

    return hl.vds.VariantDataset(vds.reference_data, vd_50_aa_cutoff)


def gq0_to_no_call(vds):
    """
    Turn the GQ0s into no calls so that ANs are correct
    """
    rd = vds.reference_data
    # TODO would be better to drop these once its a dense mt?
    rd = rd.filter_entries(rd.GQ > 0)
    return hl.vds.VariantDataset(rd, vds.variant_data)


def matrix_table_ac_an_af(mt, ancestry_file):
    """
    Create a DENSE MATRIX TABLE to calculate AC, AN, AF and Sample Count
    TODO: test sample_count
    """

    sample_id_to_sub_population = create_vat_inputs.parse_ancestry_file(ancestry_file)

    mt = mt.annotate_cols(pop=hl.literal(sample_id_to_sub_population)[mt.s])
    ac_an_af_mt = mt.select_rows(
        ac_an_af=hl.agg.call_stats(mt.GT, mt.alleles),
        call_stats_by_pop=hl.agg.group_by(mt.pop, hl.agg.call_stats(mt.GT, mt.alleles)),
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


    ac_an_af_rows = ac_an_af_split.annotate_rows(
        info = hl.struct(
            AC=ac_an_af_split.row.ac_an_af.AC[ac_an_af_split.row.a_index],
            AN=ac_an_af_split.row.ac_an_af.AN,
            AF=ac_an_af_split.row.ac_an_af.AF[ac_an_af_split.row.a_index],
            Hom=ac_an_af_split.row.ac_an_af.homozygote_count[ac_an_af_split.row.a_index],
            **pop_info_fields
        )
    )

    # note that SC = AC - homozygote_count

    ht = ac_an_af_rows.rows()
    ht = ht.filter(ht.alleles[1] != "*") # remove spanning deletions
    # create a filtered sites only VCF
    hl.export_vcf(ht, sites_only_vcf_path)



def main(vds, ancestry_file_location, sites_only_vcf_path):
    transforms = [
        hard_filter_non_passing_sites,
        failing_gts_to_no_call,
        remove_too_many_alt_allele_sites,
        gq0_to_no_call
    ]
    transformed_vds=vds
    for transform in transforms:
        transformed_vds = transform(transformed_vds)

    mt = hl.vds.to_dense_mt(transformed_vds)

    with open(ancestry_file_location, 'r') as ancestry_file:
        mt = matrix_table_ac_an_af(mt, ancestry_file) # this adds subpopulation information and splits our multi-allelic rows

    # potentially in the future: merge AC, AN, AF back to the original VDS with: vds = vds_ac_an_af(mt, vds)

    # create a sites only VCF (that is hard filtered!) and that can be made into a custom annotations TSV for Nirvana to use with AC, AN, AF, SC for all subpopulations and populations
    write_sites_only_vcf(mt, sites_only_vcf_path)


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
    parser.add_argument('--vds_input_path', type=str, help='Input VDS path', default="@VDS_INPUT_PATH@")
    parser.add_argument('--sites_only_output_path', type=str, help='Output sites-only VCF path',
                        default="@SITES_ONLY_VCF_OUTPUT_PATH@")

    args = parser.parse_args()

    vds = hl.vds.read_vds(args.vds_path)
    local_ancestry_file = create_vat_inputs.download_ancestry_file(args.ancestry_file)

    main(vds, local_ancestry_file, args.sites_only_vcf)
