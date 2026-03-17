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

from google.cloud import bigquery

logging.basicConfig(stream=sys.stderr, level=logging.WARNING, format='%(levelname)s - %(message)s')
log = logging.getLogger(__name__)


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

    return None


def get_loaded_tables_and_sample_ids_from_information_schema(project_id, dataset_name):
    client = bigquery.Client(project=project_id)

    try:
        # Left outer join vet and ref_ranges partition info to the Parquet load status table to only return rows for
        # which there appears to be a loaded partition but no entry in the Parquet load status table. Similar logic
        # applies for ploidy but without looking at partitions as the ploidy table is unpartitioned.
        query = f"""

            SELECT parti.table_name AS table_name, CAST(partition_id AS INT64) AS sample_id
            FROM
                `{project_id}.{dataset_name}.INFORMATION_SCHEMA.PARTITIONS` parti
            LEFT OUTER JOIN
                `{project_id}.{dataset_name}.parquet_load_status` load_status
            ON
                parti.table_name = load_status.table_name AND
                CAST(partition_id AS INT64) = load_status.sample_id
            WHERE
                REGEXP_CONTAINS(parti.table_name, "^ref_ranges_[0-9]+$|^vet_[0-9]+$") AND
                parti.total_logical_bytes > 0 AND
                load_status.table_name IS NULL

            UNION ALL

            SELECT "sample_chromosome_ploidy" AS table_name, ploidy.sample_id AS sample_id
            FROM
                `{project_id}.{dataset_name}.sample_chromosome_ploidy` ploidy
            LEFT OUTER JOIN
                `{project_id}.{dataset_name}.parquet_load_status` load_status
            ON
                ploidy.sample_id = load_status.sample_id AND
                load_status.table_name = "sample_chromosome_ploidy"
            WHERE
                ploidy.chromosome = 1 * 1000 * 1000 * 1000 * 1000 AND
                load_status.sample_id IS NULL

            ORDER BY table_name, sample_id
        """
        results = client.query(query)
        return {(table_name, sample_id) for (table_name, sample_id) in results}
    except Exception as e:
        log.error(f"ERROR: Could not query partitions table: {e}")
        raise e


def get_loaded_files_from_tracking_table(project_id, dataset_name):
    """Query tracking table to get set of already-loaded file paths."""
    client = bigquery.Client(project=project_id)
    tracking_table = f"{project_id}.{dataset_name}.parquet_load_status"

    try:
        query = f"""
            SELECT file_path
            FROM `{tracking_table}`
        """
        results = client.query(query)
        loaded_files = {row.file_path for row in results}
        log.info(f"Found {len(loaded_files)} already-loaded files in tracking table")
        return loaded_files
    except Exception as e:
        log.error(f"Could not query Parquet load status tracking table: {e}")
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

    # Get already-loaded files per the tracking table.
    tracking_table_loaded_files = get_loaded_files_from_tracking_table(args.project_id, args.dataset)

    # Get already-loaded table + sample_id combinations per INFORMATION_SCHEMA.
    information_schema_loaded_tables_sample_ids = get_loaded_tables_and_sample_ids_from_information_schema(args.project_id, args.dataset)

    # Read all parquet files
    with open(args.input) as f:
        all_files = [line.strip() for line in f if line.strip()]

    log.info(f"Found {len(all_files)} total Parquet files")

    # Group by table and filter out loaded files
    table_files = defaultdict(list)
    skipped_count = 0
    unmatched_table_name_count = 0
    already_loaded_data_count = 0

    for file_path in all_files:
        # Skip if already loaded
        if file_path in tracking_table_loaded_files:
            skipped_count += 1
            continue

        # Extract table name and sample id
        table_name, sample_id = parse_table_and_sample_id_from_file_path(file_path, args.superpartitioned_table_prefixes, args.regular_table_prefixes)
        if table_name:
            table_files[table_name].append(file_path)
        else:
            unmatched_table_name_count += 1
            log.error(f"Could not determine table for: {file_path}")

        if (table_name, sample_id) in information_schema_loaded_tables_sample_ids:
            log.error(f"No entry in Parquet load status table for {file_path}, but sample_id {sample_id} appears to already have data in table {table_name}.")
            already_loaded_data_count += 1

    if unmatched_table_name_count > 0 or already_loaded_data_count > 0:
        raise ValueError(f"Error(s) examining Parquet files to load, see messages above for details.")

    log.info(f"Skipped {skipped_count} already-loaded files")
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
    sys.exit(main())
