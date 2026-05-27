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
      contig, bin_start, bin_end, superpartition,
      n_variants, median_bin_count, dropout_flag
  where dropout_flag=1 when n_variants < --dropout-fraction * median.

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
    # 5. Compute per-(superpartition, contig) median inside Hail and flag
    #    dropout bins — avoids pulling the full table to the driver at all.
    # ------------------------------------------------------------------
    print('Computing per-SP/contig medians and dropout flags in Hail ...')

    # Median over non-zero bins only, so silent centromere/telomere regions
    # don't drag the baseline down.
    nonzero_ht = bin_counts_ht.filter(bin_counts_ht.n_variants > 0)
    medians_ht = (
        nonzero_ht
        .group_by(nonzero_ht.superpartition, nonzero_ht.contig)
        .aggregate(median_bin_count=hl.float64(hl.agg.approx_quantiles(
            hl.float64(nonzero_ht.n_variants), 0.5)))
    )

    result_ht = bin_counts_ht.annotate(
        bin_end=bin_counts_ht.bin_start + args.bin_size,
        median_bin_count=hl.coalesce(
            medians_ht[bin_counts_ht.superpartition, bin_counts_ht.contig].median_bin_count,
            hl.float64(0)
        )
    )
    result_ht = result_ht.annotate(
        dropout_flag=hl.int32(
            result_ht.n_variants < args.dropout_fraction * result_ht.median_bin_count
        )
    )
    result_ht = result_ht.select(
        'contig', 'bin_start', 'bin_end', 'superpartition',
        'n_variants', 'median_bin_count', 'dropout_flag'
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
    n_total, n_flagged = result_ht.aggregate(
        (hl.agg.count(), hl.agg.count_where(result_ht.dropout_flag == 1))
    )
    print(f'Done. Wrote {n_total:,} rows ({n_flagged:,} dropout bins flagged) to {args.output_path}')

    if n_flagged > 0:
        top_sps = (
            result_ht
            .filter(result_ht.dropout_flag == 1)
            .group_by('superpartition')
            .aggregate(n_dropout_bins=hl.agg.count())
            .order_by(hl.desc('n_dropout_bins'))
            .head(20)
            .collect()
        )
        print('\nTop superpartitions by dropout bin count:')
        for row in top_sps:
            print(f'  SP {row.superpartition}: {row.n_dropout_bins} dropout bins')
