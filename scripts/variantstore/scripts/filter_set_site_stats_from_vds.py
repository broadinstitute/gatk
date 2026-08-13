"""
Compute filter-set site statistics from a GVS Hail VDS.

If run with `--skip-non-ref-filter` this should be efficient enough to run in Terra Hail cluster
with two workers even for AoU-scale VDSes, no massive autoscaling cluster required.

This produces a table equivalent to the following BigQuery query run against
the filter_set_sites table:

    SELECT DISTINCT filters, COUNT(*)
    FROM `<project>.<dataset>.filter_set_sites`
    WHERE filter_set_name = '<filter_set_name>'
    GROUP BY 1
    ORDER BY 1

but with two important improvements:
  1. Only sites actually carried (as a non-reference call) by at least one
     sample in the VDS are counted.  This correctly handles any sites that
     were made monomorphic by sample withdrawals.
  2. The sample universe is exactly the VDS contents, so control samples and
     withdrawn samples are automatically excluded.

Usage (standalone, e.g. on a Dataproc cluster):
    python filter_set_site_stats_from_vds.py \
        --vds-path  gs://bucket/path/to/foxtrot_v9_r2.vds \
        --output-path gs://bucket/path/to/filter_site_stats.tsv

"""

import argparse
import sys

import hail as hl

# If the number of distinct filter combinations ever exceeds this threshold
# something unexpected has likely happened.  The count is printed to stderr
# so the caller is alerted even when stdout is redirected.
_FILTER_COMBO_WARN_THRESHOLD = 200


def compute_filter_set_site_stats(
        vds_path: str,
        output_path: str | None = None,
        require_non_ref: bool = True,
) -> hl.Table:
    """
    Read a GVS VDS and return a Hail Table with filter-set site statistics.

    When require_non_ref=True (default), only sites carried by at least one
    non-ref sample are counted — a safety net for sites made monomorphic by
    withdrawals but not yet pruned from the VDS.
    """
    vds = hl.vds.read_vds(vds_path)
    mt = vds.variant_data

    # -----------------------------------------------------------------
    # Optionally restrict to sites carried by ≥1 non-ref sample.
    # This is a safety net for sites that became monomorphic after
    # sample withdrawals but were not yet pruned from the VDS.
    # -----------------------------------------------------------------
    if require_non_ref:
        mt = mt.filter_rows(hl.agg.any(mt.LGT.is_non_ref()))

    stats = filter_stats_from_rows_table(mt.rows())

    if output_path is not None:
        stats.export(output_path)
    return stats


def filter_stats_from_rows_table(rows: hl.Table) -> hl.Table:
    """
    Core aggregation: rows Table with a ``filters: set<str>`` field →
    Table of (filters_string, site_count) pairs.

    Separated from VDS I/O so it can be unit-tested without a real VDS on disk.
    """
    # -----------------------------------------------------------------
    # Convert the filters set<str> to a canonical, human-readable string
    # that mirrors the BigQuery output format:
    #   - empty set  → "PASS"
    #   - non-empty  → alphabetically sorted, comma-space-separated values
    #                  e.g. "EXCESS_ALLELES, ExcessHet"
    # -----------------------------------------------------------------
    rows = rows.annotate(
        filter_str=hl.if_else(
            hl.len(rows.filters) == 0,
            'PASS',
            hl.delimit(hl.sorted(hl.array(rows.filters)), ', '),
        )
    )

    # -----------------------------------------------------------------
    # Aggregate: count sites per unique filter combination.
    # -----------------------------------------------------------------
    return (
        rows
        .group_by(rows.filter_str)
        .aggregate(n_sites=hl.agg.count())
        .rename({'filter_str': 'filters'})
        .order_by('filters')
    )


def _print_stats(stats: hl.Table) -> None:
    """Print all rows to stdout. Uses collect() over show(n) to avoid silent truncation."""
    rows = stats.collect()
    n = len(rows)

    if n > _FILTER_COMBO_WARN_THRESHOLD:
        print(
            f"WARNING: {n} distinct filter combinations found, which exceeds "
            f"the expected threshold of {_FILTER_COMBO_WARN_THRESHOLD}. "
            "Inspect the output carefully.",
            file=sys.stderr,
        )

    col_width = max(max((len(row.filters) for row in rows), default=0), len('filters'))
    header = f"{'filters':<{col_width}}  {'n_sites':>12}"
    print(header)
    print('-' * len(header))
    for row in rows:
        print(f"{row.filters:<{col_width}}  {row.n_sites:>12,}")
    print(f"\n{n} filter combination(s) total.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            'Compute filter-set site statistics from a GVS Hail VDS. '
            'Produces a TSV equivalent to the BigQuery filter_set_sites query, '
            'but restricted to sites actually carried by the samples in the VDS.'
        ),
    )
    parser.add_argument(
        '--vds-path',
        type=str,
        required=True,
        help='Path to the input VDS (GCS or local).',
    )
    parser.add_argument(
        '--output-path',
        type=str,
        required=False,
        default=None,
        help='Path to write the output TSV. Omit to print results to stdout only.',
    )
    parser.add_argument(
        '--temp-path',
        type=str,
        required=False,
        default=None,
        help='Optional path to a Hail temporary directory.',
    )
    parser.add_argument(
        '--skip-non-ref-filter',
        action='store_true',
        default=False,
        help=(
            'Skip the per-row filter that requires at least one non-ref call. '
            'Use only when the VDS is known to have no monomorphic rows.'
        ),
    )

    args = parser.parse_args()

    hl.init(
        idempotent=True,
        tmp_dir=args.temp_path,
    )
    hl.default_reference('GRCh38')

    stats = compute_filter_set_site_stats(
        vds_path=args.vds_path,
        output_path=args.output_path,
        require_non_ref=not args.skip_non_ref_filter,
    )

    # Only call collect/print when no output file was written.  stats is a lazy
    # Hail Table expression, so any additional action after export() would
    # re-execute the entire query from scratch.  When an output_path was
    # provided the results are already on disk; just echo the path.
    if args.output_path is None:
        _print_stats(stats)
    else:
        print(f"Results written to: {args.output_path}")

