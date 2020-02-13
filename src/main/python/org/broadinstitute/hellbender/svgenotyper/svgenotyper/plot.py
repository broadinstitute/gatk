from matplotlib import pyplot as plt
import numpy as np
import os
import pandas as pd
from pysam import VariantFile


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
    mean_cov = mean_cov_df.loc[samples_list].to_numpy().squeeze(-1)
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
        if not os.path.exists(os.path.dirname(figure_path)):
            os.makedirs(os.path.dirname(figure_path))
        plt.savefig(figure_path)
        plt.close(pf)
    vcf.close()
