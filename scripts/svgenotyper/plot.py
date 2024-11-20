import argparse
import os

import hail as hl
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import ternary
from scipy import interpolate


INFOS_LIST = ['PPE', 'PSR1', 'PSR2', 'EPE', 'ESR1', 'ESR2', 'PHI_PE', 'PHI_SR1', 'PHI_SR2']
SAMPLE_SITE_PLOT_DIR = 'sites'
BREAKEND_POSITIONS_PLOT_DIR = 'position_distributions'


def is_depth_only(record):
    return record.info.ALGORITHMS == 'depth'

def max_carrier_cnlp(record, cutoff: int):
    if 'CNLP' not in record.samples[0]:
        return False
    svtype = record.info['SVTYPE']
    min_cnlps = [min(record.samples[called_sample]['CNLP']) for called_sample in [s for s in record.samples if is_carrier(record.samples[s], svtype)]]
    if len(min_cnlps) == 0:
        return False
    return min(min_cnlps) <= cutoff

def median(x: list):
    if len(x) == 0:
        raise ValueError('Empty list')
    midpoint = len(x) / 2
    return sorted(x)[midpoint]


def is_carrier(record_sample, svtype):
    if svtype == 'DUP':
        return record_sample['CN'] > record_sample['NCN']
    return sum(record_sample['GT']) > 0


class GATKSVPlottingSuite(object):
    def __init__(self,
                 hail_table_path,
                 out_dir: str,
                 out_name: str,
                 image_ext: str,
                 contigs: list,
                 args):
        self.out_dir = out_dir
        self.out_name = out_name
        self.image_ext = image_ext
        self.contigs = contigs
        self.reference = hl.get_reference('GRCh38')

        mt = hl.read_matrix_table(hail_table_path)

        if args.depth_only:
            mt = mt.filter_rows(mt.info.ALGORITHMS == 'depth')
        if args.non_depth_only:
            if args.depth_only:
                raise ValueError("Both depth-only and non-depth-only filters are on")
            mt = mt.filter_rows(mt.info.ALGORITHMS != 'depth')
        if args.min_gq is not None:
            mt = mt.annotate_rows(min_gq=hl.agg.min(mt.GQ))
            mt = mt.filter_rows(mt.min_gq >= args.min_gq).drop('min_gq')
        if args.min_p_pesr is not None:
            mt = mt.filter_rows((mt.info.PPE >= args.min_p_pesr) | ((mt.info.PSR1 * mt.info.PSR2) >= args.min_p_pesr))
        if args.min_svlen is not None:
            mt = mt.filter_rows(mt.info.SVLEN >= args.min_svlen)
        if args.max_svlen is not None:
            mt = mt.filter_rows(mt.info.SVLEN <= args.max_svlen)
        if args.max_carrier_cnlp is not None:
            mt = mt.annotate_rows(passes_max_carrier_cnlp=hl.agg.all((~mt.GT.is_non_ref()) | (hl.expr.min(mt.CNLP) <= args.max_carrier_cnlp)))
            mt = mt.filter_rows(mt.passes_max_carrier_cnlp).drop('passes_max_carrier_cnlp')
        if args.max_cnlp is not None:
            mt = mt.annotate_rows(passes_max_cnlp=hl.agg.all(hl.expr.min(mt.CNLP) <= args.max_cnlp))
            mt = mt.filter_rows(mt.passes_max_cnlp).drop('passes_max_cnlp')
        self.mt = mt.cache()
        self.types = sorted(list(mt.aggregate_rows(hl.agg.collect_as_set(mt.info.SVTYPE))))
        self.samples = mt.s.collect()

    def fetch(self, contig: str = None):
        if contig is None:
            return self.mt
        else:
            return self.mt.filter_rows(self.mt.locus.contig == contig)

    def plot_sites_per_genome_summary(self):
        out_path = os.path.join(self.out_dir, self.out_name + ".sites_per_genome_counts." + self.image_ext)

        rows = 1
        cols = 2
        figure_width = max(4, len(self.types) * 2)
        figure_height = 4

        pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

        mt = self.mt
        mt = mt.annotate_cols(svtype_counts=hl.agg.counter((mt.info.SVTYPE, mt.GT.is_non_ref())))
        svtype_counts = mt.annotate_cols(svtype_counts=hl.agg.counter((mt.s, mt.info.SVTYPE, mt.GT.is_non_ref()))).svtype_counts.collect()
        counts_by_type = [[count_dict[key] for count_dict in svtype_counts for key in count_dict if key[1] == t and key[2]] for t in self.types]
        counts_by_sample = [sum([count_dict[key] for count_dict in svtype_counts for key in count_dict if key[0] == s and key[2]]) for s in self.samples]

        axes[0].boxplot(counts_by_sample, labels=['ALL'])
        axes[0].set_ylabel('sites per genome')

        axes[1].boxplot(counts_by_type, labels=self.types)
        axes[1].set_ylabel('sites per genome')
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

    def plot_contig_variant_position_histograms(self):

        def _histograms_plotter(out_path: str, contig: str, positions_contig: list, types=list):
            rows = len(types)
            cols = 1
            figure_width = 8
            figure_height = 3 * len(types)
            n_bins = 100

            pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

            for i in range(len(types)):
                if len(types) == 1:
                    ax = axes
                else:
                    ax = axes[i]
                ax.set_title("{} Breakend Distribution - {}".format(types[i], contig))
                ax.hist(positions_contig[i], bins=n_bins)
                ax.xaxis.set_ticklabels([])
                ax.set_ylabel('freq')
                if i == len(types) - 1:
                    ax.set_xlabel('position')

            plt.tight_layout()
            plt.savefig(out_path)
            plt.close(pf)

        mt = self.fetch()

        mt = mt.group_rows_by(type_contig=hl.expr.tuple((mt.info.SVTYPE, mt.locus.contig)))\
            .aggregate_rows(position=hl.agg.collect(mt.locus.position))\
            .result()
        data = mt.aggregate_rows(hl.agg.collect((mt.type_contig[0], mt.type_contig[1], mt.position)))
        positions_by_contig_dict = {c: {} for c in self.contigs}
        for d in data:
            if d[1] in positions_by_contig_dict and d[0] in self.types:
                positions_by_contig_dict[d[1]][d[0]] = d[2]
        positions_by_contig = [[positions_by_contig_dict[c][t] if t in positions_by_contig_dict[c] else [] for t in self.types] if c in positions_by_contig_dict else [] for c in self.contigs]

        for i in range(len(self.contigs)):
            variant_plots_dir = os.path.join(self.out_dir, BREAKEND_POSITIONS_PLOT_DIR)
            if not os.path.exists(variant_plots_dir):
                os.makedirs(variant_plots_dir)
            variant_pos_out_path = os.path.join(variant_plots_dir, self.out_name + ".breakend_positions." + self.contigs[i] + "." + self.image_ext)
            _histograms_plotter(out_path=variant_pos_out_path, contig=self.contigs[i], positions_contig=positions_by_contig[i], types=self.types)

    def plot_variant_counts(self):

        def _totals_plotter(out_path: str, totals: list, types: list):
            rows = 1
            cols = 1
            figure_width = max(3, len(types))
            figure_height = 4

            pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

            axes.set_title("Variant Counts")
            axes.bar(x=types, height=totals)
            axes.set_ylabel('freq')

            plt.tight_layout()
            plt.savefig(out_path)
            plt.close(pf)

        def _totals_by_contig_plotter(out_path: str, totals_by_contig: list, types: list):
            rows = len(types)
            cols = 1
            figure_width = 16
            figure_height = 3 * len(types)

            pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

            for i in range(len(types)):
                if len(types) == 1:
                    ax = axes
                else:
                    ax = axes[i]
                ax.set_title("Variant Counts - {}".format(types[i]))
                ax.bar(x=self.contigs, height=[totals_by_contig[(c, self.types[i])] if (c, self.types[i]) in totals_by_contig else 0 for c in self.contigs])
                ax.set_ylabel('freq')

            plt.tight_layout()
            plt.savefig(out_path)
            plt.close(pf)

        mt = self.fetch()
        totals = [mt.aggregate_rows(hl.agg.counter(mt.info.SVTYPE))[t] for t in self.types]
        totals_by_contig = mt.aggregate_rows(hl.agg.counter((mt.locus.contig, mt.info.SVTYPE)))

        base_out_path = os.path.join(self.out_dir, self.out_name)
        variant_counts_by_contig_out_path = base_out_path + ".variant_counts_by_contig." + self.image_ext
        _totals_by_contig_plotter(out_path=variant_counts_by_contig_out_path, totals_by_contig=totals_by_contig, types=self.types)
        variant_totals_out_path = base_out_path + ".variant_counts." + self.image_ext
        _totals_plotter(out_path=variant_totals_out_path, totals=totals, types=self.types)

    def plot_svlen_distributions(self):
        out_path = os.path.join(self.out_dir, self.out_name + ".svlen_distributions." + self.image_ext)
        mt = self.fetch()

        mt_rows = mt.rows()
        svlens = mt_rows.aggregate(hl.agg.group_by(mt_rows.info.SVTYPE, hl.agg.collect(mt_rows.info.SVLEN)))
        # Flatten if necessary
        svlens = {key: val if len(val) > 0 and not isinstance(val, list) else [x for y in val for x in y] for key,val in svlens.items()}

        rows = 1
        cols = 1
        figure_width = 4
        figure_height = 3
        max_len = max([max(svlens[key]) for key in svlens])
        stop = max(15, np.log(max_len / 50.))
        bins = 50 * np.logspace(start=0, stop=stop, num=128, base=2.0)

        pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))
        types = sorted(list(svlens.keys() - set(['BND'])))
        for type in types:
            hist, bin_edges = np.histogram(svlens[type], bins=bins)
            bin_centers = (bin_edges[0:-1] + bin_edges[1:]) / 2.0
            axes.plot(bin_centers, hist, '-')
        axes.set_xlabel('SVLEN')
        axes.set_ylabel('freq')
        axes.set_xscale('log')
        axes.set_yscale('log')
        pf.legend(labels=types, loc='right', bbox_to_anchor=(0.95, 0.8))
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

    def plot_sites(self, mean_coverage_path: str, vids: list = None, out_dir: str = 'figures'):

        def _plot_site_counts(axis, x: np.ndarray, y: np.ndarray, xlabel: str = None, ylabel: str = None, xmin: float = -0.05,
                              xmax: float = 3., ymin: float = -0.05, ymax: float = 1.05):
            l1 = axis.plot(x[gt0], y[gt0], 'o', markersize=markersize, color='g')[0]
            l2 = axis.plot(x[gt1], y[gt1], 'o', markersize=markersize, color='b')[0]
            l3 = axis.plot(x[gt2], y[gt2], 'o', markersize=markersize, color='r')[0]
            axis.set_xlabel(xlabel)
            axis.set_ylabel(ylabel)
            axis.set_xlim([xmin, xmax])
            axis.set_ylim([ymin, ymax])
            return l1, l2, l3

        rows = 3
        cols = 5
        figure_width = 12
        figure_height = 6
        markersize = 3
        image_ext = ".png"
        mean_cov_ploidy = 2.0

        mean_cov_df = pd.read_csv(mean_coverage_path, sep='\t', header=None, index_col=0)

        samples_list = list(self.header)
        mean_cov = mean_cov_df.loc[samples_list].values.squeeze(-1) / mean_cov_ploidy
        if vids is not None:
            vids_set = set(vids)
        for record in self.fetch():
            vid = record.id
            if vids is not None and vid not in vids_set:
                continue
            # TODO: skip depth calls
            if record.info['ALGORITHMS'] == 'depth':
                continue

            gt = np.asarray([sum(record.samples[sample]['GT']) for sample in samples_list])
            pe = np.asarray([record.samples[sample]['PE'] for sample in samples_list]) / mean_cov
            sr1 = np.asarray([record.samples[sample]['SR1'] for sample in samples_list]) / mean_cov
            sr2 = np.asarray([record.samples[sample]['SR2'] for sample in samples_list]) / mean_cov
            pl = np.asarray([record.samples[sample]['PL'] for sample in samples_list])
            gq = np.asarray([record.samples[sample]['GQ'] for sample in samples_list])
            if record.info['SVTYPE'] == 'DEL' or record.info['SVTYPE'] == 'DUP':
                cnlp = np.asarray([record.samples[sample]['CNLP'] for sample in samples_list])
            else:
                cnlp = None

            pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

            title = "{} {}:{}-{} {} size={}".format(vid, record.chrom, record.start, record.stop,
                                                    record.info['SVTYPE'], record.info['SVLEN'])
            pf.suptitle(title)

            gt0 = gt == 0
            gt1 = gt == 1
            gt2 = gt >= 2

            pl_max = 31
            pl[pl > pl_max] = pl_max
            pl_lim = pl_max + 1.
            count_lim = max(3., pe.max() + 0.1, sr1.max() + 0.1, sr2.max() + 0.1)
            l1, l2, l3 = _plot_site_counts(axes[0, 0], x=pe, y=pl, xlabel='Norm PE', ylabel='PL', xmax=count_lim, ymax=pl_lim)
            _plot_site_counts(axes[0, 1], x=sr1, y=pl, xlabel='Norm SR1', ylabel='PL', xmax=count_lim, ymax=pl_lim)
            _plot_site_counts(axes[0, 2], x=sr2, y=pl, xlabel='Norm SR2', ylabel='PL', xmax=count_lim, ymax=pl_lim)

            gq_max = 31
            gq[gq > gq_max] = gq_max
            gq_lim = gq_max + 1.
            _plot_site_counts(axes[1, 0], x=pe, y=gq, xlabel='Norm PE', ylabel='GQ', xmax=count_lim, ymax=gq_lim)
            _plot_site_counts(axes[1, 1], x=sr1, y=gq, xlabel='Norm SR1', ylabel='GQ', xmax=count_lim, ymax=gq_lim)
            _plot_site_counts(axes[1, 2], x=sr2, y=gq, xlabel='Norm SR2', ylabel='GQ', xmax=count_lim, ymax=gq_lim)

            if cnlp is not None:
                cnlp_max = 20
                cnlp[cnlp > cnlp_max] = cnlp_max
                cnlp_lim = cnlp_max + 1.
                cnlp_lb = -1
                _plot_site_counts(axes[2, 0], x=cnlp[..., 0], y=pl, xlabel='LP-CN=0', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim, ymax=pl_lim)
                _plot_site_counts(axes[2, 1], x=cnlp[..., 1], y=pl, xlabel='LP-CN=1', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim, ymax=pl_lim)
                _plot_site_counts(axes[2, 2], x=cnlp[..., 2], y=pl, xlabel='LP-CN=2', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim, ymax=pl_lim)
                _plot_site_counts(axes[2, 3], x=cnlp[..., 3], y=pl, xlabel='LP-CN=3', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim, ymax=pl_lim)
                _plot_site_counts(axes[2, 4], x=cnlp[..., 4], y=pl, xlabel='LP-CN=4', ylabel='PL', xmin=cnlp_lb, xmax=cnlp_lim, ymax=pl_lim)

            pf.delaxes(axes[0, 3])
            pf.delaxes(axes[0, 4])

            info_keys = ['PPE', 'PSR1', 'PSR2']
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

    def plot_sites_hist(self, mean_coverage_path: str, vids: list = []):

        def _plot_site_hist_counts(axis, x: np.ndarray, bins: np.ndarray, xlabel: str = None):
            l1 = axis.hist(np.asarray([x[gt0], x[gt1], x[gt2]]), bins=bins, stacked=True, edgecolor='k', linewidth=1)
            axis.set_xlabel(xlabel)
            axis.set_xlim([bins[0], bins[-1]])
            return l1

        rows = 3
        cols = 5
        figure_width = 12
        figure_height = 6
        image_ext = ".png"
        mean_cov_ploidy = 2.0

        mean_cov_df = pd.read_csv(mean_coverage_path, sep='\t', header=None, index_col=0)

        mt = self.fetch()
        samples_list = list(mt.s.collect())
        mean_cov = mean_cov_df.loc[samples_list].values.squeeze(-1) / mean_cov_ploidy


        records = mt.filter_rows(hl.expr.set(vids).contains(mt.rsid)).entries().collect()

        for vid in vids:
            entries = [r for r in records if r.rsid == vid]

            depth_only = is_depth_only(entries[0])

            gt = np.asarray([sum(e.GT) for e in entries])

            cnlp = [e.CNLP for e in entries if e.CNLP is not None]
            if len(cnlp) == 0:
                cnlp = None
            else:
                cnlp = np.asarray(cnlp)

            if not depth_only:
                pe = np.asarray([e.PE for e in entries]) / mean_cov
                sr1 = np.asarray([e.SR1 for e in entries]) / mean_cov
                sr2 = np.asarray([e.SR2 for e in entries]) / mean_cov

            pl = np.asarray([e.PL for e in entries])
            gq = np.asarray([e.GQ for e in entries])

            pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))

            title = "{} {}:{}-{} {} size={}".format(vid, entries[0].locus.contig, entries[0].locus.position, entries[0].info.END,
                                                    entries[0].info.SVTYPE, entries[0].info.SVLEN)
            pf.suptitle(title)

            gt0 = gt == 0
            gt1 = gt == 1
            gt2 = gt >= 2

            if not depth_only:
                pl_max = 31
                pl[pl > pl_max] = pl_max
                pl_lim = pl_max + 1.
                count_lim = max(3., pe.max() + 0.1, sr1.max() + 0.1, sr2.max() + 0.1)
                bins = np.arange(start=0, stop=count_lim, step=0.1)
                l1 = _plot_site_hist_counts(axes[0, 0], x=pe, bins=bins, xlabel='Norm PE')
                _plot_site_hist_counts(axes[0, 1], x=sr1, bins=bins, xlabel='Norm SR1')
                _plot_site_hist_counts(axes[0, 2], x=sr2, bins=bins, xlabel='Norm SR2')

            gq_max = 31
            gq[gq > gq_max] = gq_max
            gq_lim = gq_max + 1.
            bins = np.arange(start=0, stop=gq_lim, step=1)
            _plot_site_hist_counts(axes[1, 0], x=pl[..., 0], bins=bins, xlabel='PL0')
            _plot_site_hist_counts(axes[1, 1], x=pl[..., 1], bins=bins, xlabel='PL1')
            _plot_site_hist_counts(axes[1, 2], x=pl[..., 2], bins=bins, xlabel='PL2')
            _plot_site_hist_counts(axes[1, 3], x=gq, bins=bins, xlabel='GQ')

            if cnlp is not None:
                cnlp_max = 20
                cnlp[cnlp > cnlp_max] = cnlp_max
                cnlp_lim = cnlp_max + 1.
                cnlp_lb = -1
                bins = np.arange(start=0, stop=cnlp_lim, step=1)
                _plot_site_hist_counts(axes[2, 0], x=cnlp[..., 0], bins=bins, xlabel='LP-CN=0')
                _plot_site_hist_counts(axes[2, 1], x=cnlp[..., 1], bins=bins, xlabel='LP-CN=1')
                _plot_site_hist_counts(axes[2, 2], x=cnlp[..., 2], bins=bins, xlabel='LP-CN=2')
                _plot_site_hist_counts(axes[2, 3], x=cnlp[..., 3], bins=bins, xlabel='LP-CN=3')
                _plot_site_hist_counts(axes[2, 4], x=cnlp[..., 4], bins=bins, xlabel='LP-CN=4')

            pf.delaxes(axes[0, 3])

            if not depth_only:
                info_keys = ['PPE', 'PSR1', 'PSR2']
                axes[0, 4].bar(x=info_keys, height=[entries[0].info.PPE, entries[0].info.PSR1, entries[0].info.PSR2])
                axes[0, 4].set_ylim([0, 1.05])

                info_keys = ['EPE', 'ESR1', 'ESR2']
                axes[1, 4].bar(x=info_keys, height=[entries[0].info.EPE, entries[0].info.ESR1, entries[0].info.ESR2])

            labels = ['HOM_REF', 'HET', 'HOM_VAR']
            pf.legend(labels=labels, loc='center', bbox_to_anchor=(0.65, 0.8))
            plt.subplots_adjust(left=0.05, right=0.95, wspace=0.4, hspace=0.4)

            site_plots_dir = os.path.join(self.out_dir, SAMPLE_SITE_PLOT_DIR)
            if not os.path.exists(site_plots_dir):
                os.makedirs(site_plots_dir)
            figure_path = os.path.join(site_plots_dir, self.out_name + '.' + vid + image_ext)
            plt.savefig(figure_path)
            plt.close(pf)

    def plot_genotype_ternary_diagrams(self):

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
                tax.scatter(np.asarray(data)*scale, marker='D', color='g', s=4, alpha=0.2)
            tax.plot(hwe*scale, linewidth=1.0, label="HWE", color='k')
            tax.get_axes().axis('off')
            tax.clear_matplotlib_ticks()
            tax.right_corner_label("REF", fontsize=fontsize)
            tax.top_corner_label("HET", fontsize=fontsize)
            tax.left_corner_label("HOM", fontsize=fontsize)

        mt = self.fetch()
        mt = mt.filter_rows((mt.locus.contig != 'chrX') & (mt.locus.contig != 'chrY'))
        mt = mt.annotate_rows(gt_fractions=(hl.agg.fraction(mt.GT.is_hom_ref()), hl.agg.fraction(mt.GT.is_het()), hl.agg.fraction(mt.GT.is_hom_var())))

        out_path = os.path.join(self.out_dir, self.out_name + ".genotype_ternary_diagrams." + self.image_ext)
        figure_width = 4 * len(self.types)
        figure_height = 3

        pf, axes = plt.subplots(1, len(self.types), figsize=(figure_width, figure_height))
        for i in range(len(self.types)):
            gt_freq = mt.filter_rows(mt.info.SVTYPE == self.types[i]).gt_fractions.collect()
            if len(self.types) == 1:
                ax = axes
            else:
                ax = axes[i]
            _draw_ternary(gt_freq, axis=ax, label=self.types[i])
        plt.tight_layout()
        plt.savefig(out_path)

    def plot_info_field_histograms(self):
        out_path = os.path.join(self.out_dir, self.out_name + ".info_field_histograms." + self.image_ext)
        mt = self.fetch()

        mt_rows = mt.rows()
        site_infos = mt_rows.aggregate(hl.agg.group_by(mt_rows.info.SVTYPE, hl.agg.collect(tuple(mt_rows.info[i] for i in INFOS_LIST))))

        rows = len(self.types)
        cols = len(INFOS_LIST)
        figure_width = 14
        figure_height = max(2, len(self.types))

        pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))
        unit_interval_bins = np.arange(0, 1.05, 0.05)
        phi_bins = np.arange(0, 1.575, 0.075)
        bins = {
            'PPE': unit_interval_bins,
            'PSR1': unit_interval_bins,
            'PSR2': unit_interval_bins,
            'EPE': unit_interval_bins,
            'ESR1': unit_interval_bins,
            'ESR2': unit_interval_bins,
            'PHI_PE': phi_bins,
            'PHI_SR1': phi_bins,
            'PHI_SR2': phi_bins
        }
        for i in range(len(self.types)):
            for j in range(len(INFOS_LIST)):
                data = site_infos[self.types[i]]
                x = [d[j] for d in data if d[j] is not None]
                if len(self.types) == 1:
                    ax = axes[j]
                else:
                    ax = axes[i, j]
                ax.hist(x, bins=bins[INFOS_LIST[j]], density=True)
                ax.set_yticks([])
                if j == 0:
                    ax.set_ylabel(self.types[i])
                if i == 0:
                    ax.set_title(INFOS_LIST[j])
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

    def plot_allele_freq(self):
        out_path = os.path.join(self.out_dir, self.out_name + ".allele_frequency." + self.image_ext)
        mt = self.fetch()

        num_samples = mt.cols().count()
        # Autosomal only, TODO : hard-coded allosome
        mt = hl.variant_qc(mt).rows()
        mt = mt.filter((mt.locus.contig != 'chrX') & (mt.locus.contig != 'chrY'))
        vqc = mt.aggregate(hl.agg.group_by(mt.info.SVTYPE, hl.agg.collect(hl.struct(AF=mt.variant_qc.AF, N_NON_REF=mt.variant_qc.n_non_ref, P_HWE=mt.variant_qc.p_value_hwe, MIN_GQ=mt.variant_qc.gq_stats.min))))

        types = sorted(vqc.keys())

        rows = len(types)
        cols = 5
        figure_width = 20
        figure_height = max(5, 3 * rows)

        bins_lin_af = min(num_samples, 30)
        bins_log_start_af = - np.log(float(num_samples)) / np.log(10.)
        bins_log_af = np.logspace(start=bins_log_start_af, stop=0, num=min(num_samples, 30), base=10.0)

        bins_log_start_p_hwe = np.log(min([min([x.P_HWE for x in vqc[t]]) for t in vqc])) / np.log(10.)
        bins_log_p_hwe = np.logspace(start=max(-6, bins_log_start_p_hwe), stop=0, num=30, base=10.0)

        bins_min_gq = np.arange(31)

        alt_allele_index = 1  # TODO: does not work for multi-allelic, assumes alt index
        rwidth = 0.9

        pf, axes = plt.subplots(rows, cols, figsize=(figure_width, figure_height))
        for i in range(rows):
            allele_freq = [x.AF[alt_allele_index] for x in vqc[types[i]]]
            carrier_freq = [x.N_NON_REF / num_samples for x in vqc[types[i]]]
            rolling_af_fn = interpolate.interp1d(carrier_freq, allele_freq)
            rolling_af_x = np.arange(min(carrier_freq), max(carrier_freq), 0.01)
            rolling_af = rolling_af_fn(rolling_af_x)
            p_hwe = [x.P_HWE for x in vqc[types[i]]]
            min_gq = [x.MIN_GQ for x in vqc[types[i]]]

            if rows == 1:
                ax0 = axes[0]
                ax1 = axes[1]
                ax2 = axes[2]
                ax3 = axes[3]
                ax4 = axes[4]
            else:
                ax0 = axes[i, 0]
                ax1 = axes[i, 1]
                ax2 = axes[i, 2]
                ax3 = axes[i, 3]
                ax4 = axes[i, 4]
            ax0.hist(allele_freq, bins=bins_lin_af, rwidth=rwidth)
            ax0.set_ylabel('{} freq'.format(types[i]))

            ax1.hist(allele_freq, bins=bins_log_af, rwidth=rwidth)
            ax1.set_xscale('log')

            ax2.hist(p_hwe, bins=bins_log_p_hwe, rwidth=rwidth)
            ax2.set_xscale('log')

            ax3.hist(min_gq, bins=bins_min_gq, rwidth=rwidth)

            ax4.plot(carrier_freq, allele_freq, 'ok', markersize=2)
            ax4.plot([0, 1], [0, 1], '-k')
            ax4.plot([0, 1], [0, 0.5], '-k')
            ax4.plot(rolling_af_x, rolling_af, '-r')
            ax4.set_ylabel('{} AF'.format(types[i]))
            if i == rows - 1:
                ax0.set_xlabel('AF')
                ax1.set_xlabel('AF')
                ax2.set_xlabel('HWE p-val')
                ax3.set_xlabel('Min GQ')
                ax4.set_xlabel('Carrier Freq')
        plt.tight_layout()
        plt.savefig(out_path)
        plt.close(pf)

def parse_args():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('--hail-table', help='Hail table path', required=True)
    parser.add_argument('--out-name', help='Output base name', required=True)
    parser.add_argument('--contigs-list', help='List of contigs', required=True)
    parser.add_argument('--out-dir', help='Output directory', default='.')

    parser.add_argument('--depth-only', help='Depth-only variants filter', action='store_true')
    parser.add_argument('--non-depth-only', help='Non-depth-only variants filter', action='store_true')
    parser.add_argument('--min-gq', help='Plot variants above this min GQ (inclusive)', type=int)
    parser.add_argument('--min-p-pesr', help='Plot variants with PE and/or SR support above this probability (inclusive)', type=float)
    parser.add_argument('--min-svlen', help='Plot variants above this length (inclusive)', type=int)
    parser.add_argument('--max-svlen', help='Plot variants below this length (inclusive)', type=int)
    parser.add_argument('--max-carrier-cnlp', help='Plot variants with min CNLP below this value (inclusive) in carrier samples (also filters no-call sites and without CNLP format fields)', type=int)
    parser.add_argument('--max-cnlp', help='Plot variants with min CNLP below this value (inclusive) in all samples (also filters no-call sites and without CNLP format fields)', type=int)

    parser.add_argument('--coverage-file', help='Table of sample mean per-base coverage (required for site plots)')
    parser.add_argument('--image-ext', help='Image extension', default='png')

    parser.add_argument('--sites', nargs='+', help='If specified, plot sites with these variant IDs')

    args = parser.parse_args()

    if args.sites is not None and args.coverage_file is None:
        raise ValueError('Either --sites or --all-sites was used but --coverage-file is required in this case')

    if args.image_ext is not None and args.image_ext.startswith('.'):
        args.image_ext = args.image_ext[1:]

    return args


def main():
    args = parse_args()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    with open(args.contigs_list, 'r') as f:
        contigs = [x.strip() for x in f.readlines()]

    suite = GATKSVPlottingSuite(hail_table_path=args.hail_table, out_dir=args.out_dir, out_name = args.out_name,
                                image_ext=args.image_ext, contigs=contigs, args=args)

    if args.sites is not None:
        print("Creating per-site plots...")
        suite.plot_sites_hist(mean_coverage_path=args.coverage_file, vids=args.sites)
    else:
        print("Creating sites per genome count plot...")
        suite.plot_sites_per_genome_summary()

        print("Creating variant count plots...")
        suite.plot_variant_counts()

        print("Creating info field histograms...")
        suite.plot_info_field_histograms()

        print("Creating AF plot...")
        suite.plot_allele_freq()

        print("Creating genotype ternary diagrams...")
        suite.plot_genotype_ternary_diagrams()

        print("Creating SVLEN distribution plot...")
        suite.plot_svlen_distributions()

        print("Creating contig breakend position distribution plots...")
        suite.plot_contig_variant_position_histograms()


if __name__ == "__main__":
    main()
