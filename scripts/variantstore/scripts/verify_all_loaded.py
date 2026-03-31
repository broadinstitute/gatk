#!/usr/bin/env python3
"""
Verify that all Parquet files in GCS have been loaded to BigQuery.

Idempotency is determined by inspecting BigQuery directly:
- For vet_% and ref_ranges_% tables: checks INFORMATION_SCHEMA.PARTITIONS for non-empty partitions.
- For sample_chromosome_ploidy: checks for the sample_id in the table itself.

This approach does not require or consult a parquet_load_status tracking table.
"""

import argparse
import json
import logging
import os
import sys

from parse_and_group_files import (
    parse_table_and_sample_id_from_file_path,
    get_already_loaded_tables_and_sample_ids,
)

log = logging.getLogger(__name__)


def verify_all_loaded(project_id, dataset_name, gcs_files_list, output_dir,
                      superpartitioned_table_prefixes=None, regular_table_prefixes=None):
    """
    Compare GCS-derived (table_name, sample_id) pairs against what is actually
    present in BigQuery to find loads that are missing.

    Args:
        project_id: BigQuery project ID
        dataset_name: BigQuery dataset name
        gcs_files_list: Path to file listing all GCS Parquet URIs
        output_dir: Directory to write results
        superpartitioned_table_prefixes: Prefixes for superpartitioned tables (default: ["vet", "ref_ranges"])
        regular_table_prefixes: Prefixes for regular tables (default: ["sample_chromosome_ploidy"])

    Returns:
        Dictionary with verification results
    """
    if superpartitioned_table_prefixes is None:
        superpartitioned_table_prefixes = ["vet", "ref_ranges"]
    if regular_table_prefixes is None:
        regular_table_prefixes = ["sample_chromosome_ploidy"]

    os.makedirs(output_dir, exist_ok=True)

    # Read GCS file list
    with open(gcs_files_list) as f:
        gcs_files = [line.strip() for line in f if line.strip()]

    log.info(f"Found {len(gcs_files)} files in GCS")

    # Parse each GCS file path to determine its (table_name, sample_id).
    # Files whose path cannot be parsed are counted as unmatched and reported.
    gcs_pairs_to_files = {}   # (table_name, sample_id) -> list of file paths
    unmatched_files = []

    for file_path in gcs_files:
        table_name, sample_id = parse_table_and_sample_id_from_file_path(
            file_path, superpartitioned_table_prefixes, regular_table_prefixes
        )
        if table_name is None:
            unmatched_files.append(file_path)
            continue
        key = (table_name, sample_id)
        gcs_pairs_to_files.setdefault(key, []).append(file_path)

    if unmatched_files:
        log.warning(f"Could not parse {len(unmatched_files)} file path(s); they will be excluded from verification:")
        for p in unmatched_files[:20]:
            log.warning(f"  {p}")
        if len(unmatched_files) > 20:
            log.warning(f"  ... and {len(unmatched_files) - 20} more")

    total_files = len(gcs_files)
    all_gcs_pairs = set(gcs_pairs_to_files.keys())
    log.info(f"Parsed {len(all_gcs_pairs)} unique (table_name, sample_id) pairs from GCS file list")

    # Query BigQuery for all already-loaded (table_name, sample_id) pairs.
    log.info("Querying BigQuery for loaded data...")
    loaded_pairs = get_already_loaded_tables_and_sample_ids(
        project_id, dataset_name,
        superpartitioned_table_prefixes=superpartitioned_table_prefixes,
        regular_table_prefixes=regular_table_prefixes,
    )

    # Determine which GCS (table_name, sample_id) pairs are not yet in BigQuery.
    missing_pairs = all_gcs_pairs - loaded_pairs

    # Count files as loaded/missing based on their pair's presence in BigQuery.
    loaded_files_count = sum(
        len(files) for pair, files in gcs_pairs_to_files.items() if pair not in missing_pairs
    )
    missing_files_count = total_files - loaded_files_count - len(unmatched_files)

    log.info("Results:")
    log.info(f"  Total files in GCS:       {total_files}")
    log.info(f"  Unmatched (unparseable):  {len(unmatched_files)}")
    log.info(f"  Files with loaded pairs:  {loaded_files_count}")
    log.info(f"  Files with missing pairs: {missing_files_count}")
    log.info(f"  Missing (table, sample_id) pairs: {len(missing_pairs)}")

    all_loaded = (len(missing_pairs) == 0 and len(unmatched_files) == 0)

    # Write list of missing file paths if there are any
    missing_files_list_path = None
    if missing_pairs:
        missing_files_list_path = f"{output_dir}/missing_files.txt"
        missing_file_paths = sorted(
            path
            for pair in sorted(missing_pairs)
            for path in gcs_pairs_to_files.get(pair, [])
        )
        with open(missing_files_list_path, 'w') as f:
            for path in missing_file_paths:
                f.write(f"{path}\n")
        log.error(f"Wrote {len(missing_file_paths)} missing file path(s) to {missing_files_list_path}")
        log.error(f"Missing (table_name, sample_id) pairs:")
        for pair in sorted(missing_pairs)[:20]:
            log.error(f"  {pair}")
        if len(missing_pairs) > 20:
            log.error(f"  ... and {len(missing_pairs) - 20} more")

    results_dict = {
        "all_loaded": all_loaded,
        "total_files": total_files,
        "loaded_files": loaded_files_count,
        "missing_files": missing_files_count,
        "missing_files_list": missing_files_list_path,
        "unmatched_files": len(unmatched_files),
    }

    results_file = f"{output_dir}/verification_results.json"
    with open(results_file, 'w') as f:
        json.dump(results_dict, f, indent=2)

    log.info(f"Results written to {results_file}")
    return results_dict


def main():
    parser = argparse.ArgumentParser(
        description="Verify all Parquet files have been loaded to BigQuery"
    )
    parser.add_argument("--project-id", required=True, help="BigQuery project ID")
    parser.add_argument("--dataset-name", required=True, help="BigQuery dataset name")
    parser.add_argument(
        "--gcs-files-list",
        required=True,
        help="File containing list of all GCS Parquet URIs"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write verification results"
    )
    parser.add_argument(
        "--superpartitioned-table-prefixes",
        nargs="+",
        default=["vet", "ref_ranges"],
        help="Table prefixes for superpartitioned tables (default: vet ref_ranges)"
    )
    parser.add_argument(
        "--regular-table-prefixes",
        nargs="+",
        default=["sample_chromosome_ploidy"],
        help="Table prefixes for regular (non-superpartitioned) tables (default: sample_chromosome_ploidy)"
    )

    args = parser.parse_args()

    results = verify_all_loaded(
        project_id=args.project_id,
        dataset_name=args.dataset_name,
        gcs_files_list=args.gcs_files_list,
        output_dir=args.output_dir,
        superpartitioned_table_prefixes=args.superpartitioned_table_prefixes,
        regular_table_prefixes=args.regular_table_prefixes,
    )

    if results["all_loaded"]:
        log.info("✓ SUCCESS: All files have been loaded!")
    else:
        missing_count = results.get("missing_files", 0) or 0
        unmatched_count = results.get("unmatched_files", 0) or 0

        reasons = []
        if missing_count:
            reasons.append(f"{missing_count} file(s) not yet loaded")
        if unmatched_count:
            reasons.append(
                f"{unmatched_count} file(s) could not be parsed or matched to a table/sample_id"
            )

        if reasons:
            log.error("✗ INCOMPLETE: " + "; ".join(reasons))
        else:
            # Fallback in case all_loaded is False but no counts are available.
            log.error("✗ INCOMPLETE: Verification failed for unknown reasons")

        if results.get("missing_files_list"):
            log.error(f"  See missing files list: {results['missing_files_list']}")
        if results.get("unmatched_files_list"):
            log.error(f"  See unmatched files list: {results['unmatched_files_list']}")

    if not results["all_loaded"]:
        sys.exit(1)

    return 0


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format='%(levelname)s - %(message)s')
    sys.exit(main())
