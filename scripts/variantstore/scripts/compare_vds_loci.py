import argparse
import hail as hl


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Count loci present in a "new" VDS but absent from an "old" VDS.')
    parser.add_argument('--new-vds-path', type=str, required=True,
                        help='Path to the newer VDS whose novel loci will be counted.')
    parser.add_argument('--old-vds-path', type=str, required=True,
                        help='Path to the older VDS used as the reference set of loci.')
    parser.add_argument('--temp-path', type=str, required=True,
                        help='Path to a GCS directory for Hail temporary files.')

    args = parser.parse_args()

    hl.init(tmp_dir=f'{args.temp_path}/hail_tmp_general')

    new_vds = hl.vds.read_vds(args.new_vds_path)
    old_vds = hl.vds.read_vds(args.old_vds_path)

    new_loci = new_vds.variant_data.rows().select().key_by('locus')
    old_loci = old_vds.variant_data.rows().select().key_by('locus')

    new_only = new_loci.anti_join(old_loci)
    count = new_only.count()

    print(f'Loci in new VDS ({args.new_vds_path}) but not in old VDS ({args.old_vds_path}): {count}')

