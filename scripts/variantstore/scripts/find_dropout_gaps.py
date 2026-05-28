"""
find_dropout_gaps.py
--------------------
Mine the SPARSE variant_dropouts.tsv (produced before the zero-fill fix)
to identify regions where bins are completely absent for a given
superpartition.  The absence is the dropout signal.

A "gap" is a contiguous stretch of expected bins — spaced bin_size apart —
that are entirely missing from the file for a specific (superpartition,
contig).  Gaps that appear in nearly all superpartitions at the same
position are typically genuine low-variant genomic regions (centromeres,
telomeres, segmental duplications); gaps limited to a small number of SPs
are the likely data dropouts we care about.

Usage
-----
python find_dropout_gaps.py \
    --input   variant_dropouts.tsv \
    --output  dropout_gaps.tsv \
    [--bin-size       20000]   \
    [--min-gap-bp     100000]  \
    [--max-sp-fraction 0.2]

Outputs
-------
dropout_gaps.tsv — one row per (superpartition, contig, gap), columns:
    superpartition, contig, gap_start, gap_end, gap_size_bp,
    n_missing_bins, n_sps_with_same_gap, frac_sps_with_same_gap,
    likely_dropout
"""

import argparse
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# Chromosome sort key (numeric for chr1-chr22, then chrX, chrY, chrM)
# ---------------------------------------------------------------------------
_CONTIG_ORDER = {f'chr{i}': i for i in range(1, 23)}
_CONTIG_ORDER.update({'chrX': 23, 'chrY': 24, 'chrM': 25})

def _contig_key(c: str) -> int:
    return _CONTIG_ORDER.get(c, 99)


# ---------------------------------------------------------------------------
# Core gap-finding logic
# ---------------------------------------------------------------------------

def find_gaps(df: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    """
    For every (superpartition, contig), walk the sorted bin_start sequence
    and record every stretch where consecutive entries are more than
    bin_size apart.

    Returns a DataFrame with columns:
        superpartition, contig, gap_start, gap_end,
        gap_size_bp, n_missing_bins
    """
    records = []
    for (sp, contig), grp in df.groupby(['superpartition', 'contig'], sort=False):
        bins = grp['bin_start'].sort_values().values
        for prev, nxt in zip(bins[:-1], bins[1:]):
            diff = int(nxt) - int(prev)
            if diff > bin_size:
                gap_start = int(prev) + bin_size
                gap_end   = int(nxt)          # exclusive: first present bin after gap
                records.append(dict(
                    superpartition  = int(sp),
                    contig          = contig,
                    gap_start       = gap_start,
                    gap_end         = gap_end,
                    gap_size_bp     = gap_end - gap_start,
                    n_missing_bins  = (diff // bin_size) - 1,
                ))

    if not records:
        return pd.DataFrame(columns=[
            'superpartition', 'contig', 'gap_start', 'gap_end',
            'gap_size_bp', 'n_missing_bins',
        ])
    return pd.DataFrame(records)


def annotate_cross_sp(gaps: pd.DataFrame, n_sps: int) -> pd.DataFrame:
    """
    For each gap, count how many SPs have an identical gap at the same
    (contig, gap_start, gap_end).  Gaps shared by many SPs are likely
    genuine low-variant genomic regions, not data dropouts.
    """
    counts = (
        gaps
        .groupby(['contig', 'gap_start', 'gap_end'])
        .size()
        .rename('n_sps_with_same_gap')
        .reset_index()
    )
    gaps = gaps.merge(counts, on=['contig', 'gap_start', 'gap_end'], how='left')
    gaps['frac_sps_with_same_gap'] = gaps['n_sps_with_same_gap'] / n_sps
    return gaps


def flag_likely_dropouts(gaps: pd.DataFrame, max_sp_fraction: float) -> pd.DataFrame:
    """Mark gaps present in <= max_sp_fraction of all SPs as likely dropouts."""
    gaps['likely_dropout'] = gaps['frac_sps_with_same_gap'] <= max_sp_fraction
    return gaps


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Find dropout gaps in a sparse variant-dropout TSV'
    )
    parser.add_argument('--input',  required=True,
                        help='Path to sparse variant_dropouts.tsv')
    parser.add_argument('--output', required=True,
                        help='Output TSV of gap regions')
    parser.add_argument('--bin-size', type=int, default=20_000,
                        help='Bin size used when generating the input file (default: 20000)')
    parser.add_argument('--min-gap-bp', type=int, default=100_000,
                        help='Only report gaps of at least this many bp (default: 100000)')
    parser.add_argument('--max-sp-fraction', type=float, default=0.2,
                        help='Gaps present in more than this fraction of SPs are '
                             'labelled non-dropout; likely centromeres/telomeres '
                             '(default: 0.2)')
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Load
    # ------------------------------------------------------------------
    print(f'Reading {args.input} ...', flush=True)
    df = pd.read_csv(args.input, sep='\t')
    required = {'superpartition', 'contig', 'bin_start'}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f'ERROR: missing columns in input: {missing}')

    n_sps = df['superpartition'].nunique()
    print(f'  {len(df):,} rows  |  {n_sps} superpartitions  |  '
          f'{df["contig"].nunique()} contigs', flush=True)

    # ------------------------------------------------------------------
    # 2. Find gaps
    # ------------------------------------------------------------------
    print(f'Scanning for gaps > {args.bin_size:,} bp ...', flush=True)
    gaps = find_gaps(df, args.bin_size)

    if gaps.empty:
        print('No gaps found — file may already be fully expanded.')
        gaps.to_csv(args.output, sep='\t', index=False)
        return

    print(f'  {len(gaps):,} raw gaps found', flush=True)

    # ------------------------------------------------------------------
    # 3. Cross-SP annotation and dropout flagging
    # ------------------------------------------------------------------
    gaps = annotate_cross_sp(gaps, n_sps)
    gaps = flag_likely_dropouts(gaps, args.max_sp_fraction)

    # ------------------------------------------------------------------
    # 4. Filter to minimum size and sort
    # ------------------------------------------------------------------
    gaps = gaps[gaps['gap_size_bp'] >= args.min_gap_bp].copy()
    gaps['_contig_order'] = gaps['contig'].map(_contig_key)
    gaps = (
        gaps
        .sort_values(
            ['likely_dropout', 'gap_size_bp', '_contig_order', 'gap_start'],
            ascending=[False, False, True, True],
        )
        .drop(columns=['_contig_order'])
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------------
    # 5. Report
    # ------------------------------------------------------------------
    n_dropout = gaps['likely_dropout'].sum()
    n_genomic = len(gaps) - n_dropout
    print(f'\n  {len(gaps):,} gaps ≥ {args.min_gap_bp:,} bp after filter:')
    print(f'    {n_dropout:,} likely DROPOUT  '
          f'(≤ {args.max_sp_fraction:.0%} of SPs share same gap)')
    print(f'    {n_genomic:,} likely genomic  '
          f'(> {args.max_sp_fraction:.0%} of SPs share same gap)')

    print(f'\nTop likely dropout gaps:')
    header = (f"  {'SP':>4}  {'contig':6}  {'start':>12}  {'end':>12}  "
              f"{'size_bp':>10}  {'bins':>5}  {'SPs w/gap':>9}")
    print(header)
    print('  ' + '-' * (len(header) - 2))
    for _, r in gaps[gaps['likely_dropout']].head(30).iterrows():
        print(f"  {r['superpartition']:>4d}  {r['contig']:6s}  "
              f"{r['gap_start']:>12,}  {r['gap_end']:>12,}  "
              f"{r['gap_size_bp']:>10,}  {r['n_missing_bins']:>5,}  "
              f"{r['n_sps_with_same_gap']:>4}/{n_sps}")

    gaps.to_csv(args.output, sep='\t', index=False)
    print(f'\nWrote {len(gaps):,} rows to {args.output}')


if __name__ == '__main__':
    main()

