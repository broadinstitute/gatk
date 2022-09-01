import argparse


def generate_hail_export_tieout_vcf_script(vds_input_path, vcf_output_path):
    hail_script = f"""

import hail as hl
vds = hl.vds.read_vds('{vds_input_path}')
    
from datetime import datetime
start = datetime.now()
current_time = start.strftime("%H:%M:%S")
print("Start Time =", current_time)

## * Hard filter out non-passing sites !!!!!TODO ok wait, DO WE want to do this? I guess it depends on what we are handing off and when.
# note: no AC/AN and AF for filtered out positions
vd = vds.variant_data
filtered_vd = vd.filter_rows(hl.len(vd.filters)==0)
filtered_vds = hl.vds.VariantDataset(vds.reference_data, filtered_vd) # now we apply it back to the vds 

## * Replace LGT with GT ( for easier calculations later )
filtered_vd = filtered_vds.variant_data
filtered_vd = filtered_vd.annotate_entries(GT=hl.vds.lgt_to_gt(filtered_vd.LGT, filtered_vd.LA))
filtered_vds = hl.vds.VariantDataset(filtered_vds.reference_data, filtered_vd)

## * Respect the FT flag by setting all failing GTs to a no call
# TODO We dont seem to be using the dense matrix table here (TODO do we need to?)

# Logic for assigning non passing GTs as no-calls to ensure that AC,AN and AF respect the FT flag
# filtered_vd.FT is True => GT keeps its current value
# filtered_vd.FT is False => GT assigned no-call
# filtered_vd.FT is missing => GT keeps its current value
filtered_vd = filtered_vd.annotate_entries(GT=hl.or_missing(hl.coalesce(filtered_vd.FT, True), filtered_vd.GT))

# TODO drop LGT now that it will be different than the GT

## * Turn the GQ0s into no calls so that ANs are correct
rd = filtered_vds.reference_data
rd = rd.filter_entries(rd.GQ > 0) ## would be better to drop these once its a dense mt? 
filtered_vds = hl.vds.VariantDataset(rd, filtered_vd)

## * Create a DENSE MATRIX TABLE to calculate AC, AN, AF and TODO: Sample Count
mt = hl.vds.to_dense_mt(filtered_vds)
mt = hl.variant_qc(mt)
mt = mt.annotate_rows(AC=mt.variant_qc.AC, AN=mt.variant_qc.AN, AF=mt.variant_qc.AF)
mt = mt.drop('variant_qc')

mt = hl.vds.to_dense_mt(vds)
fail_case = 'FAIL'
mt = mt.annotate_entries(FT=hl.if_else(mt.FT, 'PASS', fail_case))
hl.export_vcf(mt, '{vcf_output_path}')    
        
"""
    return hail_script


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Generate script to write tieout VCF from an input Hail VDS.')
    parser.add_argument('--vds_path', type=str, help='GCS location for input VDS', required=True)
    parser.add_argument('--vcf_output_path', type=str, help='GCS location for VCF output generated from VDS',
                        required=True)

    args = parser.parse_args()

    script = generate_hail_export_tieout_vcf_script(args.vds_path, args.vcf_output_path)
    print(script)