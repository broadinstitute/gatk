#!/usr/bin/env python3
"""
Discovers Parquet files in GCS, groups them by target BigQuery table,
and filters out files that have already been loaded.
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

from google.cloud import bigquery


def parse_table_from_path(file_path, superpartitioned_table_prefixes, regular_table_prefixes):
    """
    Extract table name from GCS path.
    
    Example:
        gs://bucket/vet/001/vet_001_sample.parquet -> vet_001
        gs://bucket/ref_ranges/042/ref_042_sample.parquet -> ref_042
        gs://bucket/sample_chromosome_ploidy/sample_chromosome_ploidy.parquet -> sample_chromosome_ploidy
    """
    for prefix in superpartitioned_table_prefixes:
        # Match pattern: {prefix}/{digits}/
        pattern = rf'{prefix}/(\d+)/'
        match = re.search(pattern, file_path)
        if match:
            table_number = match.group(1)
            # ref_ranges maps to ref_XXX tables
            # Why would you remap ref_ranges to ref??
            # table_prefix = "ref" if prefix == "ref_ranges" else prefix
            # return f"{table_prefix}_{table_number}"
            return f"{prefix}_{table_number}"
    return None


def get_loaded_files(project_id, dataset_name):
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
        print(f"Found {len(loaded_files)} already-loaded files in tracking table")
        return loaded_files
    except Exception as e:
        print(f"Warning: Could not query tracking table: {e}")
        print("Assuming no files have been loaded yet")
        return set()


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
    
    # Get already-loaded files
    loaded_files = get_loaded_files(args.project_id, args.dataset)
    
    # Read all parquet files
    with open(args.input) as f:
        all_files = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(all_files)} total Parquet files")
    
    # Group by table and filter out loaded files
    table_files = defaultdict(list)
    skipped_count = 0
    unmatched_count = 0
    
    for file_path in all_files:
        # Skip if already loaded
        if file_path in loaded_files:
            skipped_count += 1
            continue
        
        # Extract table name
        table_name = parse_table_from_path(file_path, args.superpartitioned_table_prefixes, args.regular_table_prefixes)
        if table_name:
            table_files[table_name].append(file_path)
        else:
            unmatched_count += 1
            print(f"Warning: Could not determine table for: {file_path}")
    
    print(f"Skipped {skipped_count} already-loaded files")
    print(f"Could not match {unmatched_count} files to tables")
    print(f"Grouped {sum(len(files) for files in table_files.values())} files into {len(table_files)} tables")
    
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
        print(f"  {table_name}: {len(files)} files -> {fofn_path}")
    
    # Write summary outputs
    with open(output_dir / "table_names.txt", 'w') as f:
        for name in table_names:
            f.write(f"{name}\n")
    
    with open(output_dir / "fofn_paths.txt", 'w') as f:
        for path in fofn_paths:
            f.write(f"{path}\n")
    
    # Write statistics
    stats = {
        "total_files": len(all_files),
        "already_loaded": skipped_count,
        "unmatched": unmatched_count,
        "to_load": sum(len(files) for files in table_files.values()),
        "table_count": len(table_files),
        "tables": {name: len(files) for name, files in table_files.items()}
    }
    
    with open(output_dir / "stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"\nSummary written to {output_dir / 'stats.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
