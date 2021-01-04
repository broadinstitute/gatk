import argparse

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pysam


def plot_depth(record, small_posteriors_vcf, large_posteriors_vcf, rd_tabix_small, rd_tabix_large, sample_depth, args):
    if record.stop - record.pos < args.min_length:
        return
    #carriers = [s for s in record.samples if sum(record.samples[s]['GT']) > 0]
    carriers = [s for s in record.samples if record.samples[s]['RC'] == 1]
    allele_freq = len(carriers) / float(len(record.samples))
    if allele_freq > args.max_af:
        return
    if allele_freq <= args.min_af:
        return

    chrom = record.chrom
    start = record.pos
    end = record.stop

    length = end - start
    if length < args.model_size_cutoff:
        posteriors_vcf = small_posteriors_vcf
        rd_tabix = rd_tabix_small
        bin_size = 200
    else:
        posteriors_vcf = large_posteriors_vcf
        rd_tabix = rd_tabix_large
        bin_size = 2000

    header = rd_tabix.header[0].split('\t')
    padding = max([1000, 1. * length])
    query_start = max([start - padding, 1])
    query_end = end + padding

    posteriors_records = list(posteriors_vcf.fetch(contig=chrom, start=query_start, end=query_end))
    posteriors_start = [record.pos for record in posteriors_records]
    posteriors_phi = [record.info['PHI_RD'] for record in posteriors_records]
    posteriors_erd = [record.info['ERD'] for record in posteriors_records]
    posteriors_phw_g = [record.info['PHW_G'] for record in posteriors_records]
    posteriors_phw_l = [record.info['PHW_L'] for record in posteriors_records]

    posteriors_cnlp = np.asarray([[record.samples[s]['CNLP'] for s in record.samples] for record in posteriors_records])

    posteriors_carrier_cn = [[record.samples[s]['CN'] for s in carriers] for record in posteriors_records]
    posteriors_carrier_cnlp = [[sorted(record.samples[s]['CNLP']) for s in carriers] for record in posteriors_records]
    posteriors_carrier_gq = [[min(x[1] - x[0], 31) for x in y] for y in posteriors_carrier_cnlp]

    posteriors_start.append(posteriors_start[-1] + bin_size)
    posteriors_phi.append(posteriors_phi[-1])
    posteriors_erd.append(posteriors_erd[-1])
    posteriors_phw_g.append(posteriors_phw_g[-1])
    posteriors_phw_l.append(posteriors_phw_l[-1])
    posteriors_carrier_cn.append(posteriors_carrier_cn[-1])
    posteriors_carrier_cnlp.append(posteriors_carrier_cnlp[-1])
    posteriors_carrier_gq.append(posteriors_carrier_gq[-1])

    non_carriers = list(set(record.samples) - set(carriers))

    mt = pd.DataFrame(data=[x.split('\t') for x in rd_tabix.fetch(reference=chrom, start=query_start, end=query_end)], columns=header) \
        .set_index(keys=['#Chr', 'Start', 'End']) \
        .reset_index(level=['#Chr', 'End'])\
        .drop(labels=['#Chr', 'End'], axis=1)\
        .astype(dtype=np.int32)
    mt.index = mt.index.map(int)

    mt = mt.div(sample_depth.iloc[0]) * 100. / bin_size
    mt_median = mt.median(axis=1)

    rd_carriers = mt[carriers]
    rd_non_carriers = mt[non_carriers]
    rd_median = mt_median.to_numpy()
    rd_positions = mt_median.index.values

    rd_carriers_padded = np.append(rd_carriers.values, np.expand_dims(rd_carriers.values[-1, :], axis=0), axis=0)
    rd_non_carriers_padded = np.append(rd_non_carriers.values, np.expand_dims(rd_non_carriers.values[-1, :], axis=0), axis=0)
    rd_median_padded = np.append(rd_median, rd_median[-1])
    rd_positions_padded = np.append(rd_positions, rd_positions[-1] + bin_size)

    out_path = "{}.{}.depth.png".format(args.out_name, record.id)
    figure_width = 12
    figure_height = 8

    linewidth = 0.5
    alpha_carrier = 1 #max(0.5 / max(len(carriers), 1), 0.1)
    alpha_non_carrier = max(min([1, 5. / max(1, len(non_carriers))]), 0.1)
    title = "{} {} {}:{}-{}".format(record.id, record.info['SVTYPE'], chrom, start, end)

    rows = 5
    cols = 2

    grid_alpha = 0.2
    span_color = (0.9, 0.9, 0.9)
    log_thresh = 5

    fig, axes = plt.subplots(rows, cols, sharex=True, figsize=(figure_width, figure_height))

    ax = axes[0, 0]
    if len(non_carriers) > 0:
        ax.step(rd_positions_padded, rd_non_carriers_padded, color='r', linewidth=linewidth, where='post', alpha=alpha_non_carrier)
    ax.step(rd_positions_padded, rd_carriers_padded, color='b', linewidth=linewidth, where='post', alpha=alpha_carrier)
    ax.step(rd_positions_padded, rd_median_padded, color='k', linewidth=linewidth, where='post')
    ax.set_title(title)
    ax.axvspan(xmin=start, xmax=end, color=span_color)
    ax.set_xticks(posteriors_start, minor=False)
    ax.xaxis.grid(True, which='major', color='k', alpha=grid_alpha, linestyle='dotted')
    ax.tick_params(top=False, bottom=False, left=True, right=False, labelleft=True, labelbottom=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_xticks([])
    ax.set_xlabel('')
    ax.set_xlim([query_start, query_end])
    ax.set_ylabel('Norm depth')
    #if mt_median.max() > log_thresh:
    #    plt.yscale('log')

    ax = axes[1, 0]
    ax.axhline(y=2, linestyle='dashed', alpha=0.5, color='k')
    ax.axvspan(xmin=start, xmax=end, color=span_color)
    ax.set_xticks(posteriors_start, minor=False)
    ax.xaxis.grid(True, which='major', color='k', alpha=grid_alpha, linestyle='dotted')
    ax.step(posteriors_start, posteriors_carrier_cn, color='b', alpha=alpha_carrier, where='post', marker='o')
    ax.set_ylim([-0.1, 4.1])
    ax.set_yticks(np.arange(5))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_ylabel('Carrier CN')
    ax.tick_params(top=False, bottom=False, left=True, right=False, labelleft=True, labelbottom=False)

    ax = axes[2, 0]
    ax.axvspan(xmin=start, xmax=end, color=span_color)
    ax.set_xticks(posteriors_start, minor=False)
    ax.xaxis.grid(True, which='major', color='k', alpha=grid_alpha, linestyle='dotted')
    ax.step(posteriors_start, posteriors_carrier_gq, color='r', alpha=alpha_carrier, where='post', marker='o')
    ax.set_ylim([0, 32])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_ylabel('Carrier GQ')
    ax.tick_params(top=False, bottom=False, left=True, right=False, labelleft=True, labelbottom=False)

    ax = axes[3, 0]
    ax.axvspan(xmin=start, xmax=end, color=span_color)
    ax.step(posteriors_start, posteriors_phi, where='post', color='magenta', marker='o')
    ax.step(posteriors_start, posteriors_erd, where='post', color='teal', marker='o')
    ax.step(rd_positions_padded, rd_median_padded, color='k', where='post')
    ax.set_xticks(posteriors_start, minor=False)
    ax.xaxis.grid(True, which='major', color='k', alpha=grid_alpha, linestyle='dotted')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.legend(['phi', 'erd', 'rd_median'], loc='center left', bbox_to_anchor=(0.9, 0.5))
    ax.tick_params(top=False, bottom=False, left=True, right=False, labelleft=True, labelbottom=False)
    #if mt_median.max() > log_thresh:
    #    plt.yscale('log')

    ax = axes[4, 0]
    ax.axvspan(xmin=start, xmax=end, color=span_color)
    ax.step(posteriors_start, posteriors_phw_g, where='post', color='seagreen', marker='o')
    ax.step(posteriors_start, posteriors_phw_l, where='post', color='red', marker='o')
    ax.set_xticks(posteriors_start, minor=False)
    ax.set_ylim([0, 1])
    ax.xaxis.grid(True, which='major', color='k', alpha=grid_alpha, linestyle='dotted')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.legend(['phw_g', 'phw_l'], loc='center left', bbox_to_anchor=(0.9, 0.5))
    ax.tick_params(top=False, bottom=False, left=True, right=False, labelleft=True, labelbottom=False)

    j = 0
    for i in range(5):
        if i != 2:
            ax = axes[j, 1]
            ax.imshow((posteriors_cnlp[..., i] - posteriors_cnlp[..., 2]).transpose(),
                      interpolation='none', cmap='seismic', aspect='auto', vmin=-99, vmax=99,
                      extent=[posteriors_start[0], posteriors_start[-1], 0, 1])
            j += 1

    ax = axes[4, 1]
    ax.imshow((mt.to_numpy()).transpose(),
              interpolation='none', cmap='coolwarm', aspect='auto', vmin=0, vmax=2.,
              extent=[rd_positions[0], rd_positions[-1] + bin_size, 0, 1])

    plt.savefig(out_path)
    fig = ax.get_figure()
    plt.close(fig)


def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='VCF path', required=True)
    parser.add_argument('--small-depth-file', help='Small event depth file path', required=True)
    parser.add_argument('--large-depth-file', help='Large event depth file path', required=True)
    parser.add_argument('--mean-depth-file', help='Mean sample depth file path', required=True)
    parser.add_argument('--small-posteriors-vcf', help='Small CNV posteriors vcf path', required=True)
    parser.add_argument('--large-posteriors-vcf', help='Large CNV posteriors vcf path', required=True)
    parser.add_argument('--out-name', help='Output base name', required=True)

    parser.add_argument('--sites', nargs='+', help='Plot only these site(s)')

    parser.add_argument('--max-af', type=float, help='Max allele frequency (inclusive)', default=1)
    parser.add_argument('--min-af', type=float, help='Min allele frequency (exclusive)', default=0)
    parser.add_argument('--min-length', type=float, help='Min SV length', default=0)
    parser.add_argument('--model-size-cutoff', type=int, help='Size cutoff for switching depth models', default=5000)

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    vcf = pysam.VariantFile(args.vcf)
    small_posteriors_vcf = pysam.VariantFile(args.small_posteriors_vcf)
    large_posteriors_vcf = pysam.VariantFile(args.large_posteriors_vcf)
    rd_tabix_small = pysam.TabixFile(args.small_depth_file)
    rd_tabix_large = pysam.TabixFile(args.large_depth_file)
    sample_depth = pd.read_csv(args.mean_depth_file, sep='\t', index_col=0, header=None).transpose()
    for record in vcf:
        if args.sites is not None and record.id not in args.sites:
            continue
        if record.info['SVTYPE'] in ['DEL', 'DUP', 'CNV'] \
                or (record.info['SVTYPE'] == 'BND' and record.chrom == record.info['CHR2']):
            plot_depth(record, small_posteriors_vcf, large_posteriors_vcf, rd_tabix_small, rd_tabix_large, sample_depth, args)


if __name__ == "__main__":
    main()
