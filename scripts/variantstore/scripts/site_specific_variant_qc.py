import argparse
import hail as hl

def mt_for_site(vds, chrom, pos):
    vd=vds.variant_data
    fvd=vd.filter_rows((vd.locus.contig == chrom) & (vd.locus.position == pos))
    ffvd=fvd.transmute_entries(
        GT = hl.vds.lgt_to_gt(fvd.LGT, fvd.LA),
    )
    f_vds=hl.vds.VariantDataset(vds.reference_data, ffvd)
    mt = hl.vds.to_dense_mt(f_vds)
    mt = hl.variant_qc(mt)
    mt = hl.methods.split_multi(mt, left_aligned=True)
    return mt

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Site specific variant QC')
    parser.add_argument('--vds-path', type=str, help='Input VDS Path', required=True)
    parser.add_argument('--chrom', type=str, help='HG38 chromosome', required=True)
    parser.add_argument('--position', type=int, help='position', required=True)
    parser.add_argument('--temp-path', type=str, help='Path to temporary directory',
                        required=True)

    args = parser.parse_args()

    hl.init(tmp_dir=f'{args.temp_path}/hail_tmp_general')

    vds = hl.vds.read_vds(args.vds_path)
    chrom=args.chrom
    pos=args.position
    mt = mt_for_site(vds, chrom, pos)

    # integration_vds="gs://fc-358fc8a2-f84e-4d45-846a-a533d08f6103/submissions/b5622f6d-4c03-4e54-9bf4-a4b6386754f9/GvsQuickstartIntegration/ff599838-c254-4fdd-8059-1387786c5ffc/call-GvsQuickstartHailVETSIntegration/GvsQuickstartHailIntegration/822b8d15-5f20-449d-a2b5-049a33ff6576/call-GvsExtractAvroFilesForHail/GvsExtractAvroFilesForHail/6d264919-3cfa-4768-82dd-5f10f93067cd/call-OutputPath/avro/gvs_export.vds"
    # anvil_vds="gs://fc-34a42a2b-29aa-43f8-a6dd-2a50f7f3408b/submissions/15db34b0-320a-4275-8cac-9e2d19ef62af/NHGRI_AnVIL_ZERO/"

    mt.variant_qc.AC.show()
    mt.variant_qc.AN.show()
    # mt.variant_qc.AF.show()
