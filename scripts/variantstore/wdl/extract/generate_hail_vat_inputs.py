import argparse


def generate_hail_vat_inputs_script(vds_input_path,
                                    vcf_output_path,
                                    sites_only_vcf_output_path,
                                    vat_custom_annotations_tsv_path,
                                    ancestry_file_path):
    """
    What exactly does this do?
    """

    hail_script = f"""

## Get the original VDS    
vds = hl.vds.read_vds('{vds_input_path}')

from datetime import datetime
start = datetime.now()
current_time = start.strftime("%H:%M:%S")
print("Start Time =", current_time)

# Now parse the ancestry file to get it ready for the subpopulation work

import csv

# gross unicode escape of curly braces in an f string
# https://stackoverflow.com/q/42521230
sample_id_to_sub_population_map = \u007b \u007d
# TODO this is currently being treated as if it were a local file while all other parameters are GCS files
with open('{ancestry_file_path}', 'r') as file:
  reader = csv.reader(file, delimiter='\t')
  next(reader) # skip header
  for row in reader:
    key = row[0]
    value = row[4]
    sample_id_to_sub_population_map[key] = value

# print(sample_id_to_sub_population_map)
# need to make it looks like: sample_id_to_sub_population_map = {"ERS4367795":"eur","ERS4367796":"eas","ERS4367797":"eur","ERS4367798":"afr","ERS4367799":"oth","ERS4367800":"oth","ERS4367801":"oth","ERS4367802":"oth","ERS4367803":"oth","ERS4367804":"oth","ERS4367805":"oth"}

## * Hard filter out non-passing sites !!!!!TODO ok wait, DO WE want to do this? I guess it depends on what we are handing off and when.
# note: no AC/AN and AF for filtered out positions
vd = vds.variant_data
filtered_vd = vd.filter_rows(hl.len(vd.filters)==0)
filtered_vds = hl.vds.VariantDataset(vds.reference_data, filtered_vd) # now we apply it back to the vds 

## * Replace LGT with GT ( for easier calculations later )
filtered_vd = filtered_vds.variant_data
filtered_vd = filtered_vd.annotate_entries(GT=hl.vds.lgt_to_gt(filtered_vd.LGT, filtered_vd.LA) )
filtered_vds = hl.vds.VariantDataset(filtered_vds.reference_data, filtered_vd)

## * Respect the FT flag by setting all failing GTs to a no call
# Logic for assigning non passing GTs as no-calls to ensure that AC,AN and AF respect the FT flag
# filtered_vd.FT is True => GT keeps its current value
# filtered_vd.FT is False => GT assigned no-call
# filtered_vd.FT is missing => GT keeps its current value

filtered_vd = filtered_vd.annotate_entries(GT=hl.or_missing(hl.coalesce(filtered_vd.FT, True), filtered_vd.GT))
# TODO drop GT after it is used since it will not line up with the LGT


## * Turn the GQ0s into no calls so that ANs are correct
rd = filtered_vds.reference_data
rd = rd.filter_entries(rd.GQ > 0) ## would be better to drop these once its a dense mt? 
filtered_vds = hl.vds.VariantDataset(rd, filtered_vd)

## * Create a DENSE MATRIX TABLE to calculate AC, AN, AF and TODO: Sample Count
mt = hl.vds.to_dense_mt(filtered_vds)
mt = mt.annotate_cols(pop = hl.literal(sample_id_to_sub_population_map)[mt.s])
mt = mt.select_rows(
    ac_an_af =  hl.agg.call_stats(mt.GT, mt.alleles), 
    call_stats_by_pop = hl.agg.group_by(mt.pop, hl.agg.call_stats(mt.GT, mt.alleles))
)

qc_data = mt.rows() ## this is what we will use to create the TSV for Nirvana custom annotations

## Now we join this back to the VDS
filtered_vd = filtered_vds.variant_data
filtered_vd = filtered_vd.annotate_rows(ac_an_af = qc_data[filtered_vd.row_key])
final_vds = hl.vds.VariantDataset(rd, filtered_vd)

## save the VDS to GCP
(unclear that we even need this if we are only going to create a VAT with what is left)
# ???

## Create the VAT inputs:

# 1. Create a Sites only VCF
# TODO we will want to drop some cols because there's no reason to carry around some of this stuff
hl.export_vcf(final_vds.variant_data.rows(), '{sites_only_vcf_output_path}')

# 2. Create a VAT TSV file with subpopulation data
# TODO create the desired TSV-- We're gonna pull in Tim for this
# Do we need to track the dropped sites?
# filter out sites with too many alt alleles and trim extraneous INFO and FORMAT fields
# normalize, left align and split multi allelic sites to new lines, remove duplicate lines
# filter out spanning deletions and variants with an AC of 0
# drop GTs if we have not already since LGT doesn't match
# figure out how to calculate the Sample Count (is it just AC_ref + AC_var?)
# 
hl.export_table(qc_data, '{vat_custom_annotations_tsv_path}')

## The following can be used for extracting a VCF
mt = hl.vds.to_dense_mt(vds)
fail_case = 'FAIL'
mt = mt.annotate_entries(FT=hl.if_else(mt.FT, 'PASS', fail_case))
hl.export_vcf(mt, '{vcf_output_path}')
"""
    return hail_script


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Generate VAT inputs script from Hail VDS and an ancestry file.')
    parser.add_argument('--vds_path', type=str, help='GCS location for input VDS', required=True)
    parser.add_argument('--ancestry_file_path', type=str, help='Local filesystem location for input ancestry file',
                        required=True)
    parser.add_argument('--vcf_output_path', type=str, help='GCS location for VCF output generated from VDS',
                        required=True)
    parser.add_argument('--sites_only_vcf_output_path', type=str,
                        help='GCS location for Sites Only VCF output generated from VDS',
                        required=True)
    parser.add_argument('--vat_custom_annotations_tsv_path', type=str,
                        help='GCS location for output annotations TSV',
                        required=True)

    args = parser.parse_args()

    vat_script = generate_hail_vat_inputs_script(args.vds_path,
                                                 args.ancestry_file_path,
                                                 args.vcf_output_path,
                                                 args.sites_only_vcf_output_path,
                                                 args.vat_custom_annotations_tsv_path
                                                 )
    print(vat_script)
