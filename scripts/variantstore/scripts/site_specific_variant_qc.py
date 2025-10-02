import argparse
import hail as hl


def gts_for_loci(vds, loci):
    sites = [hl.locus_interval(l[0], l[1], l[1], reference_genome='GRCh38', includes_end=True) for l in loci]
    vds = hl.vds.filter_intervals(vds, sites)

    vd=vds.variant_data
    fvd=vd.transmute_entries(
        GT = hl.vds.lgt_to_gt(vd.LGT, vd.LA),
    )
    ffvd=fvd.select_entries('GT')
    f_vds=hl.vds.VariantDataset(vds.reference_data, ffvd)

    mt = hl.vds.to_dense_mt(f_vds)
    ht=mt.entries()
    gts = ht.aggregate(hl.agg.counter(ht.GT))
    return gts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Site specific variant QC')
    parser.add_argument('--input-vds-path', type=str, help='Input VDS Path', required=True)
    parser.add_argument('--loci', type=str, help='Loci to include, format is chrA:posA,chrB:posB,...', required=True)
    parser.add_argument('--temp-path', type=str, help='Path to temporary directory', required=True)

    args = parser.parse_args()

    hl.init(tmp_dir=f'{args.temp_path}/hail_tmp_general')

    input_vds = hl.vds.read_vds(args.input_vds_path)
    loci = [l[0].split(':') for l in args.loci.split(',')]
    gts = gts_for_loci(input_vds, loci)
    gts.show()
