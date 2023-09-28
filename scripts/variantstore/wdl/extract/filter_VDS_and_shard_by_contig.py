import hail as hl
import argparse


# requires multiallelic mt as input
def write_chr_vcf(mt, chrom, full_path, metadata, tabix=True):
    print(f'writing VCF for chr={chrom}')
    mt = mt.filter_rows(mt.locus.contig == chrom)
    hl.export_vcf(mt, full_path, tabix=tabix, metadata=metadata)


def make_dense_mt(vds, max_alt_alleles=None, is_keep_as_vqsr_fields=False):
    """
    Given a VDS, follow a standard method to densify into a MatrixTable.
    This process includes:
        - Rendering genotypes (GT), genotype allelic depth (AD)
        - Transmute FT to from a boolean to a string
        - Drop an artifact of Hail: the gvcf_info field
        - Put AC, AF, and AN into the standard MatrixTable convention (ie, in info) for rendering as VCF.
        - Drop local (eg, LAD, LGT) and other extraneous fields.  Please note that this code is smart enough to silently
            skip any fields that are not present without throwing an error.

    :param vds: input VDS as a gs URL
    :param max_alt_alleles (None):  Drop any variant sites that have the number of alt alleles exceeding this value.
        Specify None to not prune any variant sites by this criteria.
    :param is_keep_as_vqsr_fields (False):  If True, keep the AS_VQSR fields, but move the fields to the info field.
    :return: a dense MatrixTable
    """
    vd_gt = vds.variant_data

    # if going to filter alleles, do it here
    if max_alt_alleles is not None:
        vd_gt = vd_gt.filter_rows(hl.len(vd_gt.alleles) < max_alt_alleles)

    vd_gt = vd_gt.annotate_entries(AD = hl.vds.local_to_global(vd_gt.LAD,
                                                               vd_gt.LA, n_alleles=hl.len(vd_gt.alleles),
                                                               fill_value=0, number='R'))

    # Update the variant data to include GT and convert FT to a string.
    vd_gt = vd_gt.transmute_entries(
        GT = hl.vds.lgt_to_gt(vd_gt.LGT, vd_gt.LA),
    )

    # Some VDS do not have FT and Hail will throw an exception if you try to modify it w/ transmute.
    if 'FT' in vd_gt.entry:
        vd_gt = vd_gt.transmute_entries(FT = hl.if_else(vd_gt.FT, "PASS", "FAIL"))

    # Drop gvcf_info (if it exists) since it causes issues in test data we have seen.
    if 'gvcf_info' in vd_gt.entry:
        vd_gt = vd_gt.drop('gvcf_info')

    # Densify and calculate AC, AN, etc
    d_callset = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, vd_gt))

    # ## Add variant_qc
    # Make sure to prune duplicate fields
    d_callset = hl.variant_qc(d_callset)

    d_callset = d_callset.annotate_rows(info = hl.struct(AC=d_callset.variant_qc.AC[1:],
                                                             AF=d_callset.variant_qc.AF[1:],
                                                             AN=d_callset.variant_qc.AN,
                                                             homozygote_count=d_callset.variant_qc.homozygote_count))

    if is_keep_as_vqsr_fields and ('as_vqsr' in d_callset.row):
        d_callset = d_callset.annotate_rows(info = d_callset.info.annotate(AS_VQSLOD = d_callset.as_vqsr.values().vqslod,
                                         AS_YNG = d_callset.as_vqsr.values().yng_status))

    # Drop duplicate nested fields that are already in INFO field rendered by call_stats()
    d_callset = d_callset.annotate_rows(variant_qc = d_callset.variant_qc.drop("AC", "AF", "AN", "homozygote_count"))

    # ## Drop a bunch of unused fields
    #  'LGT', 'LA' have already been removed
    # Please note that this cell is specific to AoU VDS and fields specified here may not exist in test data.
    # This will not fail if a field is missing.
    fields_to_drop = ['as_vqsr', 'LAD',
                'tranche_data', 'truth_sensitivity_snp_threshold',
             'truth_sensitivity_indel_threshold','snp_vqslod_threshold','indel_vqslod_threshold']
    d_callset = d_callset.drop(*(f for f in fields_to_drop if f in d_callset.entry or f in d_callset.row or f in d_callset.col))
    return d_callset


# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('--vds_url')
parser.add_argument('--bed_url')
parser.add_argument('--vcf_header_url')
parser.add_argument('--contig')
parser.add_argument('--output_gs_url')
args = parser.parse_args()

# Load the arguments
vds_url = args.vds_url
bed_url = args.bed_url
vcf_header_url = args.vcf_header_url
contig = args.contig
output_gs_url = args.output_gs_url

# Load only the contig of interest.
vds1 = hl.vds.read_vds(vds_url)
vds1 = hl.vds.filter_chromosomes(vds1, keep=contig)

# Load bed file
bed_intervals = hl.import_bed(bed_url, reference_genome='GRCh38')

# Filter by bed file
callset = hl.vds.filter_intervals(vds1, bed_intervals, split_reference_blocks=True)

# Densify into a MatrixTable
mt = make_dense_mt(vds1, max_alt_alleles=50)

# Write this contig as a VCF
write_chr_vcf(mt, contig, output_gs_url, hl.get_vcf_metadata(vcf_header_url), tabix=True)
