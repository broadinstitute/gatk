"""
genome_wide_dropouts.py
-----------------------
Detect variant-call dropout regions across the entire genome for every
superpartition in a GVS Hail VDS.

Designed to run as a PySpark job on a Hail Dataproc cluster, submitted by
run_in_hail_cluster.py.  Arguments are supplied via --key value pairs derived
from a script-arguments JSON file.

Inputs
------
* A VDS (Hail Variant Dataset) at ``--vds-path``.
* A TSV file at ``--samples-path`` (GCS) produced by the BigQuery
  superpartition-sampling query.  Required columns:
      sample_name    (string)
      sample_id      (integer)
      superpartition (integer)

Output
------
* A TSV written to ``--output-path`` (GCS) with columns:
      contig, bin_start, bin_end, superpartition, n_variants,
      cross_sp_median, cross_sp_dropout_flag,
      within_sp_median, within_sp_dropout_flag

  cross_sp_median       — median n_variants at this bin position across all SPs
                          that have any variants there (the primary baseline).
  cross_sp_dropout_flag — 1 when n_variants < dropout-fraction * cross_sp_median.
                          Primary dropout signal: flags bins where a specific SP
                          has far fewer variants than other SPs at the same locus.
  within_sp_median      — median n_variants across all non-zero bins on the same
                          chromosome within the same SP (secondary baseline).
  within_sp_dropout_flag — 1 when n_variants < dropout-fraction * within_sp_median.
                          Secondary signal: flags bins sparse relative to that SP's
                          own chromosomal baseline.

Usage (direct)
--------------
python genome_wide_dropouts.py \
    --vds-path          gs://<bucket>/my.vds \
    --temp-path         gs://<bucket>/hail-temp \
    --samples-path      gs://<bucket>/superpartition_samples.tsv \
    --output-path       gs://<bucket>/dropout_bins.tsv \
    [--bin-size         50000] \
    [--dropout-fraction 0.5]
"""

import argparse

import hail as hl
import pandas as pd


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_BIN_SIZE = 50_000
DEFAULT_DROPOUT_FRACTION = 0.5  # bins below 50 % of the per-SP/contig median are flagged


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_sample_map(samples_path: str) -> tuple:
    """Read the superpartition TSV from a GCS or local path via hl.hadoop_open.

    sample_name is forced to str at read time because production GVS sample
    names are all-numeric and pandas would otherwise silently infer them as
    int64, causing hl.literal type-check failures downstream.
    """
    with hl.hadoop_open(samples_path, 'r') as f:
        df = pd.read_csv(f, sep='\t', dtype={'sample_name': str})

    required = {'sample_name', 'sample_id', 'superpartition'}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f'TSV is missing required columns: {missing}')

    sp_to_samples = (
        df.groupby('superpartition')['sample_name']
        .apply(list)
        .to_dict()
    )
    return sp_to_samples, df


def build_bin_counts(vd: hl.MatrixTable, bin_size: int) -> hl.Table:
    """
    For every (superpartition, contig, bin_start) triple, count the number of
    loci at which at least one sample in that superpartition has a defined LGT.

    Parameters
    ----------
    vd : MatrixTable
        variant_data from the VDS, already filtered to sampled columns and
        annotated with a ``superpartition`` column field (int32).
    bin_size : int
        Width of each genomic bin in base-pairs.

    Returns
    -------
    hl.Table
        Schema: { contig: str, bin_start: int, superpartition: int,
                  n_variants: int64 }
    """
    # --- 1. Collapse to one column per superpartition ----------------------
    # For each (locus, superpartition) entry, record whether *any* sampled
    # sample in that superpartition has a defined LGT at this locus.
    vd_by_sp = vd.group_cols_by(vd.superpartition).aggregate(
        any_defined=hl.agg.any(hl.is_defined(vd.LGT))
    )

    # Drop rows where no superpartition has any call (reduces shuffle size).
    vd_by_sp = vd_by_sp.filter_rows(hl.agg.any(vd_by_sp.any_defined))

    # --- 2. Flatten to entries, bin positions, count -----------------------
    entries = vd_by_sp.entries()
    entries = entries.filter(entries.any_defined)

    bin_counts = entries.group_by(
        contig=entries.locus.contig,
        bin_start=(entries.locus.position // bin_size) * bin_size,
        superpartition=entries.superpartition,
    ).aggregate(n_variants=hl.agg.count())

    return bin_counts


def build_bin_counts_for_contig(vd: hl.MatrixTable, contig: str, bin_size: int) -> hl.Table:
    """
    Same as build_bin_counts but restricted to a single contig.
    Called once per chromosome so each Spark job is chromosome-sized rather
    than whole-genome, which avoids executor OOM during the group_by shuffle.
    """
    interval = hl.parse_locus_interval(
        hl.eval(hl.format('%s:1-%d', contig,
                          hl.get_reference('GRCh38').lengths[contig])),
        reference_genome='GRCh38',
    )
    vd_chr = hl.filter_intervals(vd, [interval])
    return build_bin_counts(vd_chr, bin_size)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='Genome-wide variant dropout detector across all superpartitions',
    )
    parser.add_argument('--vds-path',          type=str, required=True,
                        help='Path to the Hail VDS (GCS)')
    parser.add_argument('--temp-path',         type=str, required=True,
                        help='Hail temporary directory (GCS path)')
    parser.add_argument('--samples-path',      type=str, required=True,
                        help='GCS path to the superpartition samples TSV '
                             '(columns: sample_name, sample_id, superpartition)')
    parser.add_argument('--output-path',       type=str, required=True,
                        help='Output TSV path (GCS)')
    parser.add_argument('--bin-size',          type=int, default=DEFAULT_BIN_SIZE,
                        help=f'Genomic bin size in bp (default: {DEFAULT_BIN_SIZE:,})')
    parser.add_argument('--dropout-fraction',  type=float, default=DEFAULT_DROPOUT_FRACTION,
                        help=f'Flag bins below this fraction of the per-SP/contig median '
                             f'(default: {DEFAULT_DROPOUT_FRACTION})')
    args = parser.parse_args()

    hl.init(tmp_dir=f'{args.temp_path}/hail_tmp_dropouts')

    # ------------------------------------------------------------------
    # 1. Load the superpartition sample map from GCS
    # ------------------------------------------------------------------
    print(f'Loading sample map from {args.samples_path} ...')
    sp_to_samples, sp_df = load_sample_map(args.samples_path)
    all_sampled = list(sp_df['sample_name'])
    print(f'  {len(sp_to_samples)} superpartitions, {len(all_sampled)} total sampled samples')

    # ------------------------------------------------------------------
    # 2. Load VDS and filter to sampled samples only using the idiomatic
    #    hl.vds.filter_samples API, which propagates the sample filter
    #    through both variant_data and reference_data more efficiently
    #    than a manual filter_cols on the raw MatrixTable.
    # ------------------------------------------------------------------
    print(f'Loading VDS from {args.vds_path} ...')
    vds = hl.vds.read_vds(args.vds_path)

    # hl.vds.filter_samples expects a plain Python list of sample names
    vds_filtered = hl.vds.filter_samples(vds, all_sampled, keep=True,
                                          remove_dead_alleles=False)
    vd = vds_filtered.variant_data
    print(f'  Filtered VDS to {vd.count_cols()} sampled columns')

    # ------------------------------------------------------------------
    # 3. Annotate columns with their superpartition number
    # ------------------------------------------------------------------
    sample_to_sp = {
        row['sample_name']: int(row['superpartition'])
        for _, row in sp_df.iterrows()
    }
    sp_map_literal = hl.literal(sample_to_sp, hl.tdict(hl.tstr, hl.tint32))
    vd = vd.annotate_cols(superpartition=sp_map_literal[vd.s])

    # ------------------------------------------------------------------
    # 4. Compute per-(superpartition, bin) variant counts ONE CHROMOSOME
    #    AT A TIME and checkpoint each result.
    #
    #    Processing the whole genome in one Spark job causes executors to
    #    OOM during the group_by shuffle (Stage 0 has 119 K partitions for
    #    the full AoU VDS).  Scattering by chromosome limits each job to
    #    ~1/25 of that data, which fits comfortably in executor memory.
    # ------------------------------------------------------------------
    CONTIGS = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

    chr_tables = []
    for contig in CONTIGS:
        chk_path = f'{args.temp_path}/bin_counts_{contig}.ht'
        if hl.hadoop_exists(chk_path + '/_SUCCESS'):
            print(f'  {contig}: checkpoint already exists, reusing ...')
            chr_tables.append(hl.read_table(chk_path))
            continue

        print(f'  {contig}: computing bin counts ...')
        ht = build_bin_counts_for_contig(vd, contig, args.bin_size)
        ht = ht.checkpoint(chk_path, overwrite=True)
        chr_tables.append(ht)

    print('Unioning per-chromosome tables ...')
    bin_counts_ht = hl.Table.union(*chr_tables)

    # ------------------------------------------------------------------
    # 5. Expand to a COMPLETE grid: every (contig, bin_start) that has
    #    data for ANY superpartition × every superpartition.
    #
    #    Strategy: group by (contig, bin_start) first (~110K unique bins)
    #    to collect each bin's SP counts into a single array row, then
    #    expand to 136 entries per bin with a local map+explode.  This
    #    avoids any cross-product shuffle; the only shuffle is the initial
    #    group_by which reduces ~20M rows → ~110K rows.
    # ------------------------------------------------------------------
    print('Building complete bin × superpartition grid (filling absent bins with 0)...')

    # Collect all SP numbers to the driver (just 136 ints — trivial)
    all_sps = sorted(bin_counts_ht.aggregate(
        hl.agg.collect_as_set(bin_counts_ht.superpartition)
    ))
    print(f'  {len(all_sps)} superpartitions, expanding grid...')
    all_sps_lit = hl.literal(all_sps, hl.tarray(hl.tint32))

    # Group: one row per (contig, bin_start) with all SP counts as an array
    bins_ht = (
        bin_counts_ht
        .group_by('contig', 'bin_start')
        .aggregate(sp_counts=hl.agg.collect(
            hl.struct(
                superpartition=bin_counts_ht.superpartition,
                n_variants=bin_counts_ht.n_variants,
            )
        ))
    )

    # Checkpoint the grouped table (~110K rows) before the explode
    bins_chk = f'{args.temp_path}/bins_grouped.ht'
    if hl.hadoop_exists(bins_chk + '/_SUCCESS'):
        print(f'  bins_grouped checkpoint already exists, reusing ...')
        bins_ht = hl.read_table(bins_chk)
    else:
        print(f'Checkpointing grouped bins to {bins_chk} ...')
        bins_ht = bins_ht.checkpoint(bins_chk, overwrite=True)

    # For each bin, generate one entry per SP; fill absent SPs with 0.
    # The map+explode is local (no shuffle) because all 136 SP counts are
    # already co-located on the same row after the group_by.
    bins_ht = bins_ht.annotate(
        entries=all_sps_lit.map(lambda sp: hl.struct(
            superpartition=sp,
            n_variants=hl.coalesce(
                hl.dict(bins_ht.sp_counts.map(
                    lambda x: (x.superpartition, x.n_variants)
                )).get(sp),
                hl.int64(0),
            )
        ))
    )
    bins_ht = bins_ht.drop('sp_counts')
    bins_ht = bins_ht.explode('entries')
    complete_ht = bins_ht.transmute(
        superpartition=bins_ht.entries.superpartition,
        n_variants=bins_ht.entries.n_variants,
    )
    complete_ht = complete_ht.key_by('contig', 'bin_start', 'superpartition')

    # Checkpoint the complete grid before the median join
    complete_chk = f'{args.temp_path}/complete_grid.ht'
    if hl.hadoop_exists(complete_chk + '/_SUCCESS'):
        print(f'  complete_grid checkpoint already exists, reusing ...')
        complete_ht = hl.read_table(complete_chk)
    else:
        print(f'Checkpointing complete grid to {complete_chk} ...')
        complete_ht = complete_ht.checkpoint(complete_chk, overwrite=True)

    # ------------------------------------------------------------------
    # 6. Compute dropout flags — two complementary signals
    #
    # within_sp_dropout_flag  (secondary)
    #   Compares each bin against the median of ALL bins on the same
    #   chromosome for the same superpartition.  Catches bins that are
    #   unusually sparse relative to that SP's own chr-level baseline.
    #
    # cross_sp_dropout_flag  (primary)
    #   Compares each bin against the median count at the SAME genomic
    #   position across ALL superpartitions.  A bin with n_variants=0
    #   for SP83 while every other SP has ~50 variants there is a clean
    #   data-dropout signal.  Genuine low-variant genomic regions
    #   (centromeres, etc.) are not flagged because they depress counts
    #   equally across all SPs.
    # ------------------------------------------------------------------
    print('Computing per-SP/contig medians and dropout flags in Hail ...')

    nonzero_ht = complete_ht.filter(complete_ht.n_variants > 0)

    # Within-SP median: group by (SP, contig)
    within_sp_medians_ht = (
        nonzero_ht
        .group_by(nonzero_ht.superpartition, nonzero_ht.contig)
        .aggregate(within_sp_median=hl.float64(hl.agg.approx_quantiles(
            hl.float64(nonzero_ht.n_variants), 0.5)))
    )

    # Cross-SP median: group by (contig, bin_start) — the median count at
    # each genomic position across all SPs that have any variants there.
    # Intentionally excludes zero-count SPs from the median so that dropout
    # bins don't pull the baseline down.
    cross_sp_medians_ht = (
        nonzero_ht
        .group_by(nonzero_ht.contig, nonzero_ht.bin_start)
        .aggregate(cross_sp_median=hl.float64(hl.agg.approx_quantiles(
            hl.float64(nonzero_ht.n_variants), 0.5)))
    )

    result_ht = complete_ht.annotate(
        bin_end=complete_ht.bin_start + args.bin_size,
        within_sp_median=hl.coalesce(
            within_sp_medians_ht[complete_ht.superpartition, complete_ht.contig].within_sp_median,
            hl.float64(0)
        ),
        cross_sp_median=hl.coalesce(
            cross_sp_medians_ht[complete_ht.contig, complete_ht.bin_start].cross_sp_median,
            hl.float64(0)
        ),
    )
    result_ht = result_ht.annotate(
        within_sp_dropout_flag=hl.int32(
            result_ht.n_variants < args.dropout_fraction * result_ht.within_sp_median
        ),
        cross_sp_dropout_flag=hl.int32(
            result_ht.n_variants < args.dropout_fraction * result_ht.cross_sp_median
        ),
    )
    # Clear the composite key so we can freely reorder fields in select
    # and then sort by the desired output order.
    result_ht = result_ht.key_by()
    result_ht = result_ht.select(
        'contig', 'bin_start', 'bin_end', 'superpartition',
        'n_variants',
        'cross_sp_median', 'cross_sp_dropout_flag',
        'within_sp_median', 'within_sp_dropout_flag',
    )
    result_ht = result_ht.order_by(
        result_ht.superpartition, result_ht.contig, result_ht.bin_start
    )

    # ------------------------------------------------------------------
    # 6. Write output TSV directly to GCS — no to_pandas() needed.
    # ------------------------------------------------------------------
    print(f'Writing results to {args.output_path} ...')
    result_ht.export(args.output_path)

    # Collect a small summary to the driver just for the printed report.
    print('Collecting summary statistics ...')
    n_total, n_cross_flagged, n_within_flagged = result_ht.aggregate((
        hl.agg.count(),
        hl.agg.count_where(result_ht.cross_sp_dropout_flag == 1),
        hl.agg.count_where(result_ht.within_sp_dropout_flag == 1),
    ))
    print(f'Done. Wrote {n_total:,} rows to {args.output_path}')
    print(f'  cross_sp_dropout_flag=1  : {n_cross_flagged:,} bins')
    print(f'  within_sp_dropout_flag=1 : {n_within_flagged:,} bins')

    if n_cross_flagged > 0:
        top_sps = (
            result_ht
            .filter(result_ht.cross_sp_dropout_flag == 1)
            .group_by('superpartition')
            .aggregate(n_dropout_bins=hl.agg.count())
            .order_by(hl.desc('n_dropout_bins'))
            .head(20)
            .collect()
        )
        print('\nTop superpartitions by cross-SP dropout bin count:')
        for row in top_sps:
            print(f'  SP {row.superpartition}: {row.n_dropout_bins} dropout bins')
