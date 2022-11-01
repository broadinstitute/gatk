import argparse
import hail as hl

import create_vat_inputs


###
# VAT preparation:
# filtering:
# hard filter bad sites: hard_filter_non_passing_sites()
# hard filter based on FT flag:
# get GT, ~replace_lgt_with_gt()~, swap to no-calls: failing_gts_to_no_call(), TODO: do we still need this with Tim's change?
# turn the GQ0s into no calls so that ANs are correct -- hmmm maybe we should check with Lee about this because we do have GQ0s in the VDS
# track how many sites have more than 100 alt alleles
# drop 100+ alternate alleles (TODO: or DONT!!!!)
# calculate the AC, AN, AF, SC, for the full population and for the subpopulations
### I think the rest of these things we're gonna do with bcftools?
# make a tsv that's one variant per row (do we want to do this with bcftools or just like grep in a WDL?)
# drop any * rows in that tsv (I dont think these exist!)
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
    TODO note that this may no longer be needed depending on Tim's changes
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
    return mt.select_rows(
        ac_an_af=hl.agg.call_stats(mt.GT, mt.alleles),
        call_stats_by_pop=hl.agg.group_by(mt.pop, hl.agg.call_stats(mt.GT, mt.alleles)),
    )


def vds_ac_an_af(mt, vds):
    """
    Join the dense matrix table input (just AC, AN and AF) back to the original VDS.
    """
    qc_data = mt.rows() ### TODO: make sure this join does not hard filter!!!!)
    vd = vds.variant_data
    vd = vd.annotate_rows(ac_an_af=qc_data[vd.row_key])
    return hl.vds.VariantDataset(vds.reference_data, vd) ## TODO: IT LOOKS LIKE IT DOES HARD FILTER!


def write_sites_only_vcf(vds, sites_only_vcf_path):
    # TODO we will want to drop some cols because there's no reason to carry around some of this stuff
    split_filtered_vds = hl.vds.split_multi(vds)
    hl.export_vcf(split_filtered_vds.variant_data.rows(), sites_only_vcf_path)


def write_vat_custom_annotations(mt, vat_custom_annotations_tsv_path):
    """
    Create a VAT TSV file with subpopulation data from the input dense matrix table.
    """
    # TODO:
    # Do we need to track the dropped sites?
    # filter out sites with too many alt alleles
    # filter out variants with an AC of 0


    # split multi allelic sites to new lines
    ac_an_af_split = hl.methods.split_multi(mt)

    #CHROM	POS	REF	ALT	AC	AN	AF	AC_Hom	AC_Het	AC_Hemi	AC_afr	AN_afr	AF_afr	AC_Hom_afr	AC_Het_afr	AC_Hemi_afr	AC_amr	AN_amr	AF_amr	AC_Hom_amr	AC_Het_amr	AC_Hemi_amr	AC_eas	AN_eas	AF_eas	AC_Hom_eas	AC_Het_eas	AC_Hemi_eas	AC_eur	AN_eur	AF_eur	AC_Hom_eur	AC_Het_eur	AC_Hemi_eur	AC_mid	AN_mid	AF_mid	AC_Hom_mid	AC_Het_mid	AC_Hemi_mid	AC_oth	AN_oth	AF_oth	AC_Hom_oth	AC_Het_oth	AC_Hemi_oth	AC_sas	AN_sas	AF_sas	AC_Hom_sas	AC_Het_sas	AC_Hemi_sas
    ac_an_af_rows = ac_an_af_split.select_rows(
        CHROM=ac_an_af_split.row.locus.contig,
        POS=ac_an_af_split.row.locus.position,
        REF=ac_an_af_split.row.alleles[0],
        ALT=ac_an_af_split.row.alleles[1],
        AC=ac_an_af_split.row.ac_an_af.AC[1],
        AN=ac_an_af_split.row.ac_an_af.AN,
        AF=ac_an_af_split.row.ac_an_af.AF[1],
        homozygote_count=ac_an_af_split.row.ac_an_af.homozygote_count[1],
        eas_AC = ac_an_af_split.call_stats_by_pop.get('eas').AC[1],
        eas_AN = ac_an_af_split.call_stats_by_pop.get('eas').AN,
        eas_AF = ac_an_af_split.call_stats_by_pop.get('eas').AF[1],
        eas_homozygote_count = ac_an_af_split.call_stats_by_pop.get('eas').homozygote_count[1],
    )

    # note that SC = AC - homozygote_count

    ac_an_af_rows=ac_an_af_rows.drop('tranche_data')
    ac_an_af_rows=ac_an_af_rows.drop('truth_sensitivity_snp_threshold')
    ac_an_af_rows=ac_an_af_rows.drop('truth_sensitivity_indel_threshold')
    ac_an_af_rows=ac_an_af_rows.drop('snp_vqslod_threshold')
    ac_an_af_rows=ac_an_af_rows.drop('indel_vqslod_threshold')

    print("Now lets export the custom annotations tsv")

    hl.Table.export(ac_an_af_rows.rows(), vat_custom_annotations_tsv_path)







def main(vds, ancestry_file_location, sites_only_vcf_path, vat_custom_annotations_tsv_path):
    transforms = [
        hard_filter_non_passing_sites,
        replace_lgt_with_gt,
        failing_gts_to_no_call,
        gq0_to_no_call
    ]

    for transform in transforms:
        transformed_vds = transform(vds)

    mt = hl.vds.to_dense_mt(transformed_vds)

    with open(ancestry_file_location, 'r') as ancestry_file:
        mt = matrix_table_ac_an_af(mt, ancestry_file)

    # merge AC, AN, AF back to the original VDS (TODO: make sure this join does not hard filter!!!!)
    vds = vds_ac_an_af(mt, vds)
    # TODO: write this vds somewhere now that is has AC, AN and AF?

    # create a sites only VCF (that is hard filtered!)
    write_sites_only_vcf(transformed_vds, sites_only_vcf_path)
    # create a custom annotations TSV for Nirvana to use with AC, AN, AF, SC for all subpopulations and populations
    write_vat_custom_annotations(mt, vat_custom_annotations_tsv_path)


def annotate_entry_filter_flag(mt):
    """
    Annotate FT flag for entries
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
    parser.add_argument('--ancestry_file', type=str, help='Input ancestry file', required=True)
    parser.add_argument('--vds_path', type=str, help='Input VDS Path', default="@VDS_INPUT_PATH@") #TODO: names feel confusing to me that this is an input and the next two are outputs?
    parser.add_argument('--sites_only_vcf', type=str, help='Output sites-only VCF file',
                        default="@SITES_ONLY_VCF_OUTPUT_PATH@")
    parser.add_argument('--vat_custom_annotations', type=str, help='Output VAT custom annotations file',
                        default="@VAT_CUSTOM_ANNOTATIONS_OUTPUT_PATH@")

    args = parser.parse_args()

    vds = hl.vds.read_vds(args.vds_path)
    # write_sites_only_vcf(vds, args.sites_only_vcf)
    local_ancestry_file = create_vat_inputs.download_ancestry_file(args.ancestry_file)

    main(vds, local_ancestry_file, args.sites_only_vcf, args.vat_custom_annotations)
