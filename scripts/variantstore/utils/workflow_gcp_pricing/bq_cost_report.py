#!/usr/bin/env python3
"""
BigQuery cost reporter.

Queries a cost_observability table and reports per-stage BigQuery costs,
broken down by Query Scanned and Storage API Scanned charges.

Usage:
    python bq_cost_report.py <project.dataset> [<project.dataset> ...]

Example:
    python bq_cost_report.py aou-genomics-curation-prod.parquet_20k_scale_test_2
"""

import argparse
import re
import sys
from typing import Optional

from google.cloud import bigquery

COST_OBSERVABILITY_TABLE = "cost_observability"

# ── Pricing ──────────────────────────────────────────────────────────────────
QUERY_SCANNED_RATE_PER_TIB = 6.25   # USD / TiB
STORAGE_API_RATE_PER_TIB   = 1.10   # USD / TiB
GIB_PER_TIB                = 1024.0

# ── Stage definitions ─────────────────────────────────────────────────────────
# Each tuple: (display_name, step, event_key_filter, cost_type)
#   event_key_filter=None  → no event_key WHERE clause (treated as Query Scanned)
#   cost_type              → "query" | "storage_api"
STAGES = [
    (
        "Populate Alt Allele",
        "CreateAltAlleles",
        "BigQuery Query Scanned",
        "query",
    ),
    (
        "Create Filter Set – Query Scanned",
        "GvsCreateFilterSet",
        "BigQuery Query Scanned",
        "query",
    ),
    (
        "Create Filter Set – Storage API",
        "GvsCreateFilterSet",
        "Storage API Scanned",
        "storage_api",
    ),
    (
        "Prepare Ranges",
        "GvsPrepareRanges",
        None,           # no event_key filter → all bytes; treated as Query Scanned
        "query",
    ),
    (
        "Extract – Storage API",
        "GvsExtractCallset",
        "Storage API Scanned",
        "storage_api",
    ),
]

# Strictly allow project.dataset identifiers and prevent malformed table names.
DATASET_PATTERN = re.compile(r"^[A-Za-z0-9:.\-]+\.[A-Za-z0-9_]+$")


# ── Helpers ───────────────────────────────────────────────────────────────────

def build_sql(table: str, step: str, event_key: Optional[str]) -> str:
    sql = f"""
SELECT
  ROUND(SUM(event_bytes) / (1024 * 1024 * 1024), 2) AS gib
FROM
  `{table}`
WHERE
  step = '{step}'"""
    if event_key:
        sql += f"\n  AND event_key = '{event_key}'"
    return sql.strip()


def cost_for_gib(gib: float, cost_type: str) -> float:
    tib = gib / GIB_PER_TIB
    rate = QUERY_SCANNED_RATE_PER_TIB if cost_type == "query" else STORAGE_API_RATE_PER_TIB
    return tib * rate


def rate_label(cost_type: str) -> tuple[str, float]:
    if cost_type == "query":
        return "Query Scanned", QUERY_SCANNED_RATE_PER_TIB
    return "Storage API Scanned", STORAGE_API_RATE_PER_TIB


def validate_dataset_identifier(dataset: str) -> None:
    if not DATASET_PATTERN.fullmatch(dataset):
        raise ValueError(
            f"Invalid dataset '{dataset}'. Expected format: project.dataset "
            "with only letters/numbers/._:- in project and letters/numbers/_ in dataset."
        )


# ── Per-table report ──────────────────────────────────────────────────────────

def report_table(client: bigquery.Client, table: str) -> float:
    """Run all stage queries against *table*, print results, return total cost."""

    print(f"\n{'=' * 72}")
    print(f"  Table: {table}")
    print(f"{'=' * 72}")

    total_cost = 0.0

    for display_name, step, event_key, cost_type in STAGES:
        sql = build_sql(table, step, event_key)
        try:
            rows = list(client.query(sql).result())
            gib = float(rows[0].gib) if rows and rows[0].gib is not None else 0.0
        except Exception as exc:
            raise RuntimeError(f"{display_name}: query failed for table '{table}': {exc}") from exc

        tib = gib / GIB_PER_TIB
        cost = cost_for_gib(gib, cost_type)
        total_cost += cost
        type_name, rate = rate_label(cost_type)

        print(f"\n  {display_name}")
        print(f"    Step  : {step}")
        print(f"    Type  : {type_name}  (${rate:.2f}/TiB)")
        print(f"    Size  : {gib:>12,.2f} GiB  →  {tib:>9,.4f} TiB")
        print(f"    Cost  : ${cost:>10,.4f}")

    print(f"\n  {'─' * 50}")
    print(f"  Table total:  ${total_cost:,.4f}")
    return total_cost


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report per-stage BigQuery costs from a cost_observability table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "datasets",
        nargs="+",
        metavar="project.dataset",
        help="BigQuery project.dataset path(s) to query (table is always cost_observability).",
    )
    parser.add_argument(
        "--project",
        metavar="PROJECT_ID",
        help="GCP project ID to use for the BigQuery client (billing/job project).",
    )
    args = parser.parse_args()

    try:
        for dataset in args.datasets:
            validate_dataset_identifier(dataset)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(2)

    client = bigquery.Client(project=args.project)

    tables = [f"{ds}.{COST_OBSERVABILITY_TABLE}" for ds in args.datasets]

    grand_total = 0.0
    try:
        for table in tables:
            grand_total += report_table(client, table)
    except Exception as exc:
        print(f"\nError: {exc}", file=sys.stderr)
        print("Aborted before totals could be completed.", file=sys.stderr)
        sys.exit(1)

    if len(tables) > 1:
        print(f"\n{'=' * 72}")
        print(f"  GRAND TOTAL ({len(tables)} tables):  ${grand_total:,.4f}")

    print(f"\n{'=' * 72}\n")


if __name__ == "__main__":
    main()
