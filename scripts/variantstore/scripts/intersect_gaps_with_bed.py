"""
intersect_gaps_with_bed.py
--------------------------
Intersect the dropout gaps TSV (from find_dropout_gaps.py) with a BED file
of expected-coverage regions.  For each gap, compute how many bases fall
within BED-covered intervals and flag gaps where a significant fraction of
the gap lies in covered territory — those are genuine data-loss concerns.

The two key outputs per gap:
  covered_bp        — bases of the gap that overlap ≥1 BED interval
  covered_fraction  — covered_bp / gap_size_bp
  genuine_concern   — True when covered_fraction >= --min-covered-fraction

Example results:
  SP64 chr19:39980000-40640000  covered=660000/660000 (100%)  → genuine concern
  SP81 chr9:43400000-60520000   covered=10000/17120000  (0.1%) → not a concern

Usage
-----
python intersect_gaps_with_bed.py \
    --gaps    dropout_gaps.tsv \
    --bed     coverage.bed \
    --output  annotated_gaps.tsv \
    [--min-covered-fraction 0.1] \
    [--likely-dropout-only]       # only process rows where likely_dropout=True

BED file format assumed: tab-separated, columns contig / start (0-based) / end.
Additional columns (name, score, etc.) are ignored.  The script auto-detects
whether the BED file has a header row.
"""

import argparse
import bisect
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# BED loading and indexing
# ---------------------------------------------------------------------------

def load_bed(bed_path: str) -> pd.DataFrame:
    """
    Load a BED file into a DataFrame with columns: contig, start, end.
    Handles files with or without a header and with any number of extra
    columns.  BED coordinates are 0-based half-open [start, end).
    """
    with open(bed_path) as f:
        first = f.readline()

    # Detect header: if the first field of the first line is not a known
    # chromosome prefix assume it is a header row.
    has_header = not (
        first.startswith('chr') or
        first.split('\t')[0].isdigit()
    )
    df = pd.read_csv(
        bed_path,
        sep='\t',
        header=0 if has_header else None,
        dtype=str,
    )
    # Keep only the first three columns regardless of how many exist
    df = df.iloc[:, :3].copy()
    df.columns = ['contig', 'start', 'end']
    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)
    return df


def build_bed_index(bed_df: pd.DataFrame) -> dict:
    """
    Build a per-contig sorted list of (start, end) tuples for fast overlap
    queries.  Overlapping/adjacent intervals are merged first so that the
    covered_bases calculation doesn't double-count.
    """
    index = {}
    for contig, grp in bed_df.groupby('contig'):
        # Sort and merge overlapping intervals
        intervals = sorted(zip(grp['start'].values, grp['end'].values))
        merged = []
        for s, e in intervals:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append([s, e])
        index[contig] = merged
    return index


# ---------------------------------------------------------------------------
# Overlap calculation
# ---------------------------------------------------------------------------

def covered_bases(gap_start: int, gap_end: int, intervals: list) -> int:
    """
    Return the number of bases in [gap_start, gap_end) that are covered by
    at least one interval in ``intervals`` (a sorted list of [start, end]
    pairs as produced by build_bed_index).

    Uses bisect to skip intervals that end before the gap starts.
    """
    if not intervals:
        return 0

    # Find the index of the first interval whose END > gap_start
    # (any interval ending at or before gap_start cannot overlap the gap)
    ends = [e for _, e in intervals]
    lo = bisect.bisect_right(ends, gap_start)

    covered = 0
    for i in range(lo, len(intervals)):
        s, e = intervals[i]
        if s >= gap_end:
            break  # all remaining intervals start after the gap ends
        overlap = min(gap_end, e) - max(gap_start, s)
        if overlap > 0:
            covered += overlap
    return covered


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Intersect dropout gaps with a BED coverage file'
    )
    parser.add_argument('--gaps',   required=True,
                        help='dropout_gaps.tsv from find_dropout_gaps.py')
    parser.add_argument('--bed',    required=True,
                        help='BED file of expected-coverage regions')
    parser.add_argument('--output', required=True,
                        help='Output annotated TSV')
    parser.add_argument('--min-covered-fraction', type=float, default=0.1,
                        help='Flag gap as genuine_concern when covered_fraction '
                             '>= this value (default: 0.1 = 10%%)')
    parser.add_argument('--likely-dropout-only', action='store_true',
                        help='Only process rows where likely_dropout == True')
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # 1. Load inputs
    # ------------------------------------------------------------------
    print(f'Loading gaps from {args.gaps} ...', flush=True)
    gaps = pd.read_csv(args.gaps, sep='\t')
    required = {'superpartition', 'contig', 'gap_start', 'gap_end', 'gap_size_bp'}
    missing = required - set(gaps.columns)
    if missing:
        sys.exit(f'ERROR: gaps file is missing columns: {missing}')

    if args.likely_dropout_only and 'likely_dropout' in gaps.columns:
        n_before = len(gaps)
        gaps = gaps[gaps['likely_dropout'] == True].copy()
        print(f'  Filtered to likely_dropout=True: {len(gaps):,} / {n_before:,} rows')
    else:
        print(f'  {len(gaps):,} gaps loaded')

    print(f'Loading BED from {args.bed} ...', flush=True)
    bed_df = load_bed(args.bed)
    print(f'  {len(bed_df):,} BED intervals on {bed_df["contig"].nunique()} contigs')

    # ------------------------------------------------------------------
    # 2. Build BED index
    # ------------------------------------------------------------------
    bed_index = build_bed_index(bed_df)
    total_covered_bp = sum(e - s for intervals in bed_index.values()
                           for s, e in intervals)
    print(f'  {total_covered_bp:,} total covered bases after merging overlaps')

    # ------------------------------------------------------------------
    # 3. Compute overlap for each gap
    # ------------------------------------------------------------------
    print('Computing gap × BED overlaps ...', flush=True)

    covered_bps      = []
    covered_fracs    = []
    genuine_concerns = []

    for _, row in gaps.iterrows():
        contig     = row['contig']
        gap_start  = int(row['gap_start'])
        gap_end    = int(row['gap_end'])
        gap_size   = int(row['gap_size_bp'])

        intervals  = bed_index.get(contig, [])
        cov        = covered_bases(gap_start, gap_end, intervals)
        frac       = cov / gap_size if gap_size > 0 else 0.0

        covered_bps.append(cov)
        covered_fracs.append(frac)
        genuine_concerns.append(frac >= args.min_covered_fraction)

    gaps['covered_bp']       = covered_bps
    gaps['covered_fraction'] = covered_fracs
    gaps['genuine_concern']  = genuine_concerns

    # ------------------------------------------------------------------
    # 4. Sort: genuine concerns first, then by covered_fraction descending
    # ------------------------------------------------------------------
    gaps['_contig_order'] = gaps['contig'].map(
        {f'chr{i}': i for i in range(1, 23)} | {'chrX': 23, 'chrY': 24, 'chrM': 25}
    ).fillna(99)
    gaps = (
        gaps
        .sort_values(
            ['genuine_concern', 'covered_fraction', 'gap_size_bp'],
            ascending=[False, False, False],
        )
        .drop(columns=['_contig_order'])
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------------
    # 5. Report
    # ------------------------------------------------------------------
    n_genuine  = gaps['genuine_concern'].sum()
    n_expected = len(gaps) - n_genuine
    print(f'\n  {len(gaps):,} gaps processed:')
    print(f'    {n_genuine:,}  genuine concerns '
          f'(≥ {args.min_covered_fraction:.0%} of gap in BED-covered territory)')
    print(f'    {n_expected:,}  expected / in uncovered regions')

    if n_genuine > 0:
        print(f'\nGenuine concerns (covered_fraction ≥ {args.min_covered_fraction:.0%}):')
        header = (f"  {'SP':>4}  {'contig':6}  {'start':>12}  {'end':>12}  "
                  f"{'size_bp':>10}  {'covered_bp':>10}  {'covered%':>8}")
        print(header)
        print('  ' + '-' * (len(header) - 2))
        for _, r in gaps[gaps['genuine_concern']].iterrows():
            print(f"  {r['superpartition']:>4d}  {r['contig']:6s}  "
                  f"{r['gap_start']:>12,}  {r['gap_end']:>12,}  "
                  f"{r['gap_size_bp']:>10,}  {r['covered_bp']:>10,}  "
                  f"{r['covered_fraction']:>7.1%}")

    gaps.to_csv(args.output, sep='\t', index=False)
    print(f'\nWrote {len(gaps):,} annotated rows to {args.output}')


if __name__ == '__main__':
    main()

