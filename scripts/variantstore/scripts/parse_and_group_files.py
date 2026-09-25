#!/usr/bin/env python3
"""
Discovers Parquet files in GCS, groups them by target BigQuery table,
and filters out files that have already been loaded.
"""

import argparse
import json
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    from google.cloud import bigquery
except ImportError:
    bigquery = None  # type: ignore  # will fail at runtime if BigQuery calls are made without the package installed


log = logging.getLogger(__name__)


_PREFIX_RE = re.compile(r'^[a-zA-Z][A-Za-z0-9_]+$')


def _validate_table_prefixes(superpartitioned_table_prefixes, regular_table_prefixes):
    """
    Raise ValueError if any prefix does not match '^[a-zA-Z][A-Za-z0-9_]+$'.

    Prefixes are embedded in SQL queries and regex patterns, so invalid values
    could produce malformed queries or silently incorrect behavior.
    """
    invalid = [
        p for p in (superpartitioned_table_prefixes + regular_table_prefixes)
        if not _PREFIX_RE.match(p)
    ]
    if invalid:
        raise ValueError(
            f"Table prefix(es) contain invalid characters (must match "
            f"'^[a-zA-Z][A-Za-z0-9_]+$'): {invalid}"
        )


def parse_table_and_sample_id_from_file_path(file_path, superpartitioned_table_prefixes=None, regular_table_prefixes=None):
    """
    Extract table name and sample_id from superpartitioned or regular (non-superpartitioned) GCS paths.

    Example:
        "gs://bucket/ref_ranges/ref_ranges_001_1_input_vcf_0_ERS4367795.vcf.gz.parquet" -> ("ref_ranges_001", 1)
        "gs://bucket/vet/vet_123_4567_input_vcf_0_ERS4367795.vcf.gz.parquet" -> ("vet_123", 4567)
        "gs://bucket/sample_chromosome_ploidy/sample_chromosome_ploidy_1_filename.parquet" -> ("sample_chromosome_ploidy", 1)
    """
    if regular_table_prefixes is None:
        regular_table_prefixes = ["sample_chromosome_ploidy"]
    if superpartitioned_table_prefixes is None:
        superpartitioned_table_prefixes = ["vet", "ref_ranges"]

    _validate_table_prefixes(superpartitioned_table_prefixes, regular_table_prefixes)

    for prefix in superpartitioned_table_prefixes:
        # Parse prefix + superpartition, sample id out of filename.
        pattern = f'/({prefix}_[0-9]+)_([0-9]+)_[^/]+$'
        match = re.search(pattern, file_path)
        if match:
            table = match.group(1)
            sample_id = int(match.group(2))
            return table, sample_id

    for prefix in regular_table_prefixes:
        # Parse prefix, sample id out of filename.
        pattern = rf'/({prefix})_([0-9]+)[^/]+$'
        match = re.search(pattern, file_path)
        if match:
            table = match.group(1)
            sample_id = int(match.group(2))
            return table, sample_id

    return None, None


def get_already_loaded_tables_and_sample_ids(project_id, dataset_name,
                                              superpartitioned_table_prefixes=None,
                                              regular_table_prefixes=None):
    """
    Query BigQuery directly to find which (table_name, sample_id) pairs are already loaded.

    For superpartitioned tables (vet_%, ref_ranges_%) this inspects
    INFORMATION_SCHEMA.PARTITIONS for non-empty partitions.
    For regular tables (sample_chromosome_ploidy) this queries the table itself.

    Returns a set of (table_name, sample_id) tuples.
    """
    if superpartitioned_table_prefixes is None:
        superpartitioned_table_prefixes = ["vet", "ref_ranges"]
    if regular_table_prefixes is None:
        regular_table_prefixes = ["sample_chromosome_ploidy"]

    _validate_table_prefixes(superpartitioned_table_prefixes, regular_table_prefixes)

    sub_queries = []

    if superpartitioned_table_prefixes:
        # Build a regex that matches any of the superpartitioned table names, e.g.
        # "^vet_[0-9]+$|^ref_ranges_[0-9]+$"
        superpartitioned_regex = "|".join(
            f"^{prefix}_[0-9]+$" for prefix in superpartitioned_table_prefixes
        )
        sub_queries.append(f"""
            SELECT parti.table_name AS table_name, SAFE_CAST(partition_id AS INT64) AS sample_id
            FROM `{project_id}.{dataset_name}.INFORMATION_SCHEMA.PARTITIONS` parti
            WHERE
                REGEXP_CONTAINS(parti.table_name, '{superpartitioned_regex}') AND
                parti.total_logical_bytes > 0 AND
                NOT STARTS_WITH(partition_id, '__') AND
                SAFE_CAST(partition_id AS INT64) IS NOT NULL
        """)

    for prefix in regular_table_prefixes:
        sub_queries.append(f"""
            SELECT DISTINCT '{prefix}' AS table_name, t.sample_id AS sample_id
            FROM `{project_id}.{dataset_name}.{prefix}` t
        """)

    if not sub_queries:
        return set()

    query = " UNION ALL ".join(sub_queries) + " ORDER BY table_name, sample_id"

    try:
        client = bigquery.Client(project=project_id)
        results = client.query(query)
        loaded = {(row.table_name, row.sample_id) for row in results}
        log.info(f"Found {len(loaded)} already-loaded (table_name, sample_id) pairs in BigQuery")
        return loaded
    except Exception as e:
        log.error(f"ERROR: Could not query BigQuery for loaded data: {e}")
        raise


def main():
    parser = argparse.ArgumentParser(
        description="Group Parquet files by BigQuery table and filter loaded files"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="File containing list of Parquet file paths (one per line)"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write table-specific FOFNs"
    )
    parser.add_argument(
        "--project-id",
        required=True,
        help="BigQuery project ID"
    )
    parser.add_argument(
        "--dataset",
        required=True,
        help="BigQuery dataset name"
    )
    parser.add_argument(
        "--regular-table-prefixes",
        nargs="+",
        default=["sample_chromosome_ploidy"],
        help="Table prefixes to scan for which correspond to regular, non-superpartitioned tables (e.g., sample_chromosome_ploidy)"
    )
    parser.add_argument(
        "--superpartitioned-table-prefixes",
        nargs="+",
        default=["vet", "ref_ranges"],
        help="Table prefixes which correspond to superpartitioned tables (e.g., vet ref_ranges)"
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Determine which (table_name, sample_id) pairs are already loaded in BigQuery.
    # This is the authoritative source of truth — no parquet_load_status table needed.
    already_loaded = get_already_loaded_tables_and_sample_ids(
        args.project_id, args.dataset,
        superpartitioned_table_prefixes=args.superpartitioned_table_prefixes,
        regular_table_prefixes=args.regular_table_prefixes,
    )

    # Read all parquet files
    with open(args.input) as f:
        all_files = [line.strip() for line in f if line.strip()]

    log.info(f"Found {len(all_files)} total Parquet files")

    # Group by table and filter out files whose (table_name, sample_id) is already loaded
    table_files = defaultdict(list)
    skipped_count = 0
    unmatched_table_name_count = 0

    for file_path in all_files:
        # Extract table name and sample id
        table_name, sample_id = parse_table_and_sample_id_from_file_path(
            file_path, args.superpartitioned_table_prefixes, args.regular_table_prefixes
        )
        if table_name is None:
            unmatched_table_name_count += 1
            log.error(f"Could not determine table for: {file_path}")
            continue

        # Skip if data for this (table_name, sample_id) is already present in BigQuery
        if (table_name, sample_id) in already_loaded:
            skipped_count += 1
            continue

        table_files[table_name].append(file_path)

    if unmatched_table_name_count > 0:
        raise ValueError(f"Error(s) examining Parquet files to load, see messages above for details.")

    log.info(f"Skipped {skipped_count} files for already-loaded (table_name, sample_id) pairs")
    log.info(f"Could not match {unmatched_table_name_count} files to tables")
    log.info(f"Grouped {sum(len(files) for files in table_files.values())} files into {len(table_files)} tables")

    # Write FOFNs for each table
    table_names = []
    fofn_paths = []

    for table_name in sorted(table_files.keys()):
        files = table_files[table_name]
        fofn_path = output_dir / f"{table_name}.fofn"

        with open(fofn_path, 'w') as f:
            for file_path in sorted(files):
                f.write(f"{file_path}\n")

        table_names.append(table_name)
        fofn_paths.append(str(fofn_path))
        log.info(f"  {table_name}: {len(files)} files -> {fofn_path}")

    # Write summary outputs
    with open(output_dir / "table_names.txt", 'w') as f:
        for name in sorted(table_files.keys()):
            f.write(f"{name}\n")

    with open(output_dir / "fofn_paths.txt", 'w') as f:
        for path in sorted(fofn_paths):
            f.write(f"{path}\n")

    # Write statistics
    stats = {
        "total_files": len(all_files),
        "already_loaded": skipped_count,
        "unmatched": unmatched_table_name_count,
        "to_load": sum(len(files) for files in table_files.values()),
        "table_count": len(table_files),
        "tables": {name: len(files) for name, files in table_files.items()}
    }

    with open(output_dir / "stats.json", 'w') as f:
        json.dump(stats, f, indent=2)

    log.info(f"\nSummary written to {output_dir / 'stats.json'}")
    return 0


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format='%(levelname)s - %(message)s')
    sys.exit(main())
