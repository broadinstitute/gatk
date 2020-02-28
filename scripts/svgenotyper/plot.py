import argparse
from matplotlib import pyplot as plt
import numpy as np
import os
import pandas as pd
from pysam import VariantFile
import ternary


SVTYPES_LIST = ['DEL', 'DUP', 'INS', 'INV', 'BND']
INFOS_LIST = ['PPE', 'PSR1', 'PSR2', 'PRD', 'EPE', 'ESR1', 'ESR2', 'PHI_PE', 'PHI_SR1', 'PHI_SR2']
SAMPLE_SITE_PLOT_DIR = 'sites'


def plot_sites_per_genome_summary(vcf_path: str, out_path: str):
    vcf = VariantFile(vcf_path)
    samples_list = list(vcf.header.samples)
    site_counts = {sample: {svtype: 0 for svtype in SVTYPES_LIST} for sample in samples_list}
    for record in vcf.fetch():
        svtype = record.info['SVTYPE']
        for sample in samples_list:
            if sum(record.samples[sample]['GT']) > 0:
                site_counts[sample][svtype] += 1
    vcf.close()

    rows = 1
    cols = 2
    figure_width = 8
    figure_height = 4

    pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

    total_counts = [sum([x[svtype] for svtype in SVTYPES_LIST]) for x in site_counts.values()]
    axes[0].boxplot(total_counts, labels=['ALL'])
    axes[0].set_ylabel('sites per genome')

    counts_by_type = [[x[svtype] for x in site_counts.values()] for svtype in SVTYPES_LIST]
    axes[1].boxplot(counts_by_type, labels=SVTYPES_LIST)
    axes[1].set_ylabel('sites per genome')
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(pf)


def plot_info_field_histograms(vcf_path: str, out_path: str):
    vcf = VariantFile(vcf_path)
    site_infos = {info: {svtype: [] for svtype in SVTYPES_LIST} for info in INFOS_LIST}
    for record in vcf.fetch():
        svtype = record.info['SVTYPE']
        for info in INFOS_LIST:
            site_infos[info][svtype].append(record.info[info])
    vcf.close()

    rows = len(SVTYPES_LIST)
    cols = len(INFOS_LIST)
    figure_width = 14
    figure_height = 4

    pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))
    unit_interval_bins = np.arange(0, 1.05, 0.05)
    phi_bins = np.arange(0, 1.575, 0.075)
    bins = {
        'PPE': unit_interval_bins,
        'PSR1': unit_interval_bins,
        'PSR2': unit_interval_bins,
        'PRD': unit_interval_bins,
        'EPE': unit_interval_bins,
        'ESR1': unit_interval_bins,
        'ESR2': unit_interval_bins,
        'PHI_PE': phi_bins,
        'PHI_SR1': phi_bins,
        'PHI_SR2': phi_bins
    }
    for i in range(len(SVTYPES_LIST)):
        for j in range(len(INFOS_LIST)):
            axes[i, j].hist(site_infos[INFOS_LIST[j]][SVTYPES_LIST[i]], bins=bins[INFOS_LIST[j]], density=True)
            axes[i, j].set_yticks([])
            if j == 0:
                axes[i, j].set_ylabel(SVTYPES_LIST[i])
            if i == 0:
                axes[i, j].set_title(INFOS_LIST[j])
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close(pf)


def plot_sites(vcf_path: str, mean_coverage_path: str, vids: list = None, out_dir: str = 'figures'):

    rows = 3
    cols = 5
    figure_width = 12
    figure_height = 6
    markersize = 3
    image_ext = ".png"

    def _plot_counts(axis, x: np.ndarray, y: np.ndarray, xlabel: str = None, ylabel: str = None, xmin: float = -0.05,
                     xmax: float = 3., ymin: float = -0.05, ymax: float = 1.05):
        l1 = axis.plot(x[gt0], y[gt0], 'o', markersize=markersize)[0]
        l2 = axis.plot(x[gt1], y[gt1], 'o', markersize=markersize)[0]
        l3 = axis.plot(x[gt2], y[gt2], 'o', markersize=markersize)[0]
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)
        axis.set_xlim([xmin, xmax])
        axis.set_ylim([ymin, ymax])
        return l1, l2, l3

    mean_cov_df = pd.read_csv(mean_coverage_path, sep='\t', header=None, index_col=0)

    vcf = VariantFile(vcf_path)
    samples_list = list(vcf.header.samples)
    mean_cov = mean_cov_df.loc[samples_list].values.squeeze(-1)
    if vids is not None:
        vids_set = set(vids)
    for record in vcf.fetch():
        vid = record.id
        if vids is not None and vid not in vids_set:
            continue
        # TODO: skip depth calls
        if record.info['ALGORITHMS'] == 'depth':
            continue

        gt = np.asarray([sum(record.samples[sample]['GT']) for sample in samples_list])
        pe = np.asarray([record.samples[sample]['PE'] for sample in samples_list]) / mean_cov
        sr1 = np.asarray([record.samples[sample]['SSR'] for sample in samples_list]) / mean_cov
        sr2 = np.asarray([record.samples[sample]['ESR'] for sample in samples_list]) / mean_cov
        pl = np.asarray([record.samples[sample]['PL'] for sample in samples_list])
        gq = np.asarray([record.samples[sample]['GQ'] for sample in samples_list])
        cnlp = np.asarray([record.samples[sample]['CNLP'] for sample in samples_list])

        pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

        title = "{:s} {:s}:{:d}-{:d} {:s} size={:d}".format(vid, record.chrom, record.start, record.stop,
                                                                record.info['SVTYPE'], record.info['SVLEN'])
        pf.suptitle(title)

        gt0 = gt == 0
        gt1 = gt == 1
        gt2 = gt >= 2

        count_lim = max(3., pe.max() + 0.1, sr1.max() + 0.1, sr2.max() + 0.1)
        l1, l2, l3 = _plot_counts(axes[0, 0], x=pe, y=pl, xlabel='Norm PE', ylabel='PL', xmax=count_lim)
        _plot_counts(axes[0, 1], x=sr1, y=pl, xlabel='Norm SR1', ylabel='PL', xmax=count_lim)
        _plot_counts(axes[0, 2], x=sr2, y=pl, xlabel='Norm SR2', ylabel='PL', xmax=count_lim)

        gq_lim = 105.
        _plot_counts(axes[1, 0], x=pe, y=gq, xlabel='Norm PE', ylabel='GQ', xmax=count_lim, ymax=gq_lim)
        _plot_counts(axes[1, 1], x=sr1, y=gq, xlabel='Norm SR1', ylabel='GQ', xmax=count_lim, ymax=gq_lim)
        _plot_counts(axes[1, 2], x=sr2, y=gq, xlabel='Norm SR2', ylabel='GQ', xmax=count_lim, ymax=gq_lim)

        cnlp_lim = 21.
        cnlp_lb = -1
        _plot_counts(axes[2, 0], x=cnlp[..., 0], y=pl, xlabel='LP-CN=0', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim)
        _plot_counts(axes[2, 1], x=cnlp[..., 1], y=pl, xlabel='LP-CN=1', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim)
        _plot_counts(axes[2, 2], x=cnlp[..., 2], y=pl, xlabel='LP-CN=2', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim)
        _plot_counts(axes[2, 3], x=cnlp[..., 3], y=pl, xlabel='LP-CN=3', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim)
        _plot_counts(axes[2, 4], x=cnlp[..., 4], y=pl, xlabel='LP-CN=4', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim)

        pf.delaxes(axes[0, 3])
        pf.delaxes(axes[0, 4])

        info_keys = ['PPE', 'PSR1', 'PSR2', 'PRD']
        axes[1, 3].bar(x=info_keys, height=[record.info[x] for x in info_keys])
        axes[1, 3].set_ylim([0, 1.05])

        info_keys = ['EPE', 'ESR1', 'ESR2']
        axes[1, 4].bar(x=info_keys, height=[record.info[x] for x in info_keys])

        labels = ['HOM_REF', 'HET', 'HOM_VAR']
        pf.legend(handles=[l1, l2, l3], labels=labels, loc='center', bbox_to_anchor=(0.65, 0.8))
        plt.subplots_adjust(left=0.05, right=0.95, wspace=0.4, hspace=0.4)

        figure_path = os.path.join(out_dir, vid + image_ext)
        plt.savefig(figure_path)
        plt.close(pf)
    vcf.close()


def plot_genotype_ternary_diagrams(vcf_path: str, out_path: str):

    figure_width = 16
    figure_height = 8
    gq_cutoffs = [0, 2, 30]

    def _draw_ternary(data, axis, label=""):
        scale = 100
        p = np.arange(0, 1., 0.01)
        q = 1 - p
        hwe = np.stack([p*p, 2*p*q, q*q], axis=-1)
        figure, tax = ternary.figure(ax=axis, scale=scale)
        tax.boundary(linewidth=2.0)
        tax.gridlines(multiple=20, color="gray")
        fontsize = 16.
        tax.set_title(title="%s (n=%d)" % (label, len(data)))
        if len(data) > 0:
            tax.scatter(np.asarray(data)*scale, marker='D', color='g', s=1, alpha=0.2)
        tax.plot(hwe*scale, linewidth=1.0, label="HWE", color='k')
        tax.get_axes().axis('off')
        tax.clear_matplotlib_ticks()
        tax.right_corner_label("REF", fontsize=fontsize)
        tax.top_corner_label("HET", fontsize=fontsize)
        tax.left_corner_label("HOM", fontsize=fontsize)


    vcf = VariantFile(vcf_path)
    genotypes = {svtype: [] for svtype in SVTYPES_LIST}
    gqs = {svtype: [] for svtype in SVTYPES_LIST}
    samples = vcf.header.samples
    for record in vcf.fetch():
        svtype = record.info['SVTYPE']
        genotypes[svtype].append([sum(record.samples[sample]['GT']) for sample in samples])
        gqs[svtype].append([record.samples[sample]['GQ'] for sample in samples])
    vcf.close()

    pf, axes = plt.subplots(len(gq_cutoffs), len(SVTYPES_LIST), figsize=(figure_width, figure_height))
    for i in range(len(SVTYPES_LIST)):
        svtype = SVTYPES_LIST[i]
        gt = np.asarray(genotypes[svtype])
        if gt.shape[0] > 0:
            gq = np.asarray(gqs[svtype])
            gt[gt > 2] = 2  # TODO: multiallelics probably should be filtered instead
            gt_counts = np.zeros((gt.shape[0], 3))
            for j in range(3):
                gt_copy = gt.copy()
                locs = gt_copy == j
                gt_copy[locs] = 1
                gt_copy[~locs] = 0
                gt_counts[:, j] = gt_copy.sum(axis=1)
            gt_freq = gt_counts / gt_counts.sum(axis=1, keepdims=True)
        for k in range(len(gq_cutoffs)):
            if gt.shape[0] > 0:
                gt_freq_k = gt_freq.copy()
                cutoff_locs = gq.max(axis=1) < gq_cutoffs[k]
                gt_freq_k[cutoff_locs, :] = 0
            else:
                gt_freq_k = np.zeros(np.shape([]))
            label = "{:s}, maxGQ>{:d}".format(svtype, gq_cutoffs[k])
            _draw_ternary(gt_freq_k, axis=axes[k, i], label=label)
    plt.tight_layout()
    plt.savefig(out_path)


def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', help='VCF path', required=True)
    parser.add_argument('--out-name', help='Output VCF base name', required=True)
    parser.add_argument('--out-dir', help='Output VCF directory', default='./plots')

    parser.add_argument('--coverage-file', help='Table of sample mean per-base coverage (required for site plots)')
    parser.add_argument('--image-ext', help='Image extension', default='png')

    parser.add_argument('--sites', nargs='+', help='If specified, plot sites with these variant IDs')
    parser.add_argument('--all-sites', action='store_true', help='Plot all sites')

    args = parser.parse_args()

    if args.sites is not None and args.all_sites:
        raise ValueError('Cannot specify both --sites and --all-sites')
    if (args.sites is not None or args.all_sites) and args.coverage_file is None:
        raise ValueError('Either --sites or --all-sites was used but --coverage-file is required in this case')

    if args.image_ext is not None and args.image_ext.startswith('.'):
        args.image_ext = args.image_ext[1:]

    return args


def main():
    args = parse_args()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    print("Creating sites per genome count plot...")
    sites_per_genome_summary_path = os.path.join(args.out_dir, args.out_name + ".sites_per_genome_counts." + args.image_ext)
    plot_sites_per_genome_summary(vcf_path=args.vcf, out_path=sites_per_genome_summary_path)

    print("Creating info field histograms...")
    info_field_histograms_path = os.path.join(args.out_dir, args.out_name + ".info_field_histograms." + args.image_ext)
    plot_info_field_histograms(vcf_path=args.vcf, out_path=info_field_histograms_path)

    print("Creating genotype ternary diagrams...")
    ternary_diagram_path = os.path.join(args.out_dir, args.out_name + ".genotype_ternary_diagrams." + args.image_ext)
    plot_genotype_ternary_diagrams(vcf_path=args.vcf, out_path=ternary_diagram_path)

    if args.sites is not None or args.all_sites:
        print("Creating per-site plots...")
        site_plots_dir = os.path.join(args.out_dir, SAMPLE_SITE_PLOT_DIR)
        if not os.path.exists(site_plots_dir):
            os.makedirs(site_plots_dir)
        plot_sites(vcf_path=args.vcf, mean_coverage_path=args.coverage_file, vids=args.sites, out_dir=site_plots_dir)


if __name__== "__main__":
    main()