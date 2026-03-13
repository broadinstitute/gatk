#!/usr/bin/env python3
"""
Verify that all Parquet files in GCS have been loaded to BigQuery.
Compares GCS inventory against the tracking table.
"""

import argparse
import json
import sys

from google.cloud import bigquery


def verify_all_loaded(project_id, dataset_name, gcs_files_list, output_dir):
    """
    Compare GCS files against tracking table to find missing loads.
    
    Args:
        project_id: BigQuery project ID
        dataset_name: BigQuery dataset name
        gcs_files_list: Path to file listing all GCS Parquet URIs
        output_dir: Directory to write results
    
    Returns:
        Dictionary with verification results
    """
    client = bigquery.Client(project=project_id)
    tracking_table = f"{project_id}.{dataset_name}.parquet_load_status"
    staging_table = f"{project_id}.{dataset_name}.__parquet_verification_staging"
    
    # Read GCS file list
    with open(gcs_files_list) as f:
        gcs_files = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(gcs_files)} files in GCS")
    
    # Create staging table with GCS file list
    # For large lists, we stage in BigQuery rather than using query parameters
    print(f"Creating staging table {staging_table}...")
    
    schema = [bigquery.SchemaField("file_path", "STRING")]
    
    # Create table with 1-hour expiration
    table = bigquery.Table(staging_table, schema=schema)
    table.expires = None  # Will set via OPTIONS in query
    
    # Use LOAD to populate from file or inline data
    # For simplicity, we'll use a CREATE TABLE AS SELECT with UNNEST
    rows_to_insert = [{"file_path": path} for path in gcs_files]
    
    # Split into chunks to avoid request size limits
    chunk_size = 10000
    
    # Create table with first chunk
    create_query = f"""
    CREATE OR REPLACE TABLE `{staging_table}`
    OPTIONS(
      expiration_timestamp = TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 1 HOUR)
    ) AS
    SELECT file_path
    FROM UNNEST(@file_paths) AS file_path
    """
    
    first_chunk = gcs_files[:chunk_size]
    job_config = bigquery.QueryJobConfig(
        query_parameters=[
            bigquery.ArrayQueryParameter("file_paths", "STRING", first_chunk)
        ]
    )
    
    create_job = client.query(create_query, job_config=job_config)
    create_job.result()
    
    print(f"  Created with {len(first_chunk)} rows")
    
    # Insert remaining chunks
    for i in range(chunk_size, len(gcs_files), chunk_size):
        chunk = gcs_files[i : i + chunk_size]
        
        insert_query = f"""
        INSERT INTO `{staging_table}` (file_path)
        SELECT file_path
        FROM UNNEST(@file_paths) AS file_path
        """
        
        job_config = bigquery.QueryJobConfig(
            query_parameters=[
                bigquery.ArrayQueryParameter("file_paths", "STRING", chunk)
            ]
        )
        
        insert_job = client.query(insert_query, job_config=job_config)
        insert_job.result()
        
        print(f"  Inserted {len(chunk)} rows (total: {min(i + chunk_size, len(gcs_files))})")
    
    # Compare against tracking table
    print(f"\nComparing against tracking table...")
    
    comparison_query = f"""
    WITH gcs_files AS (
      SELECT file_path FROM `{staging_table}`
    ),
    loaded_files AS (
      SELECT file_path FROM `{tracking_table}`
    )
    SELECT 
      COUNT(*) as total_files,
      COUNT(loaded_files.file_path) as loaded_files,
      COUNT(*) - COUNT(loaded_files.file_path) as missing_files
    FROM gcs_files
    LEFT JOIN loaded_files USING (file_path)
    """
    
    results = client.query(comparison_query)
    row = list(results)[0]
    
    total_files = row.total_files
    loaded_files = row.loaded_files
    missing_files = row.missing_files
    
    print(f"\nResults:")
    print(f"  Total files in GCS: {total_files}")
    print(f"  Files loaded: {loaded_files}")
    print(f"  Files missing: {missing_files}")
    
    all_loaded = (missing_files == 0)
    
    # Get list of missing files if any
    missing_files_list = None
    if not all_loaded:
        print(f"\nQuerying for missing file paths...")
        
        missing_query = f"""
        WITH gcs_files AS (
          SELECT file_path FROM `{staging_table}`
        ),
        loaded_files AS (
          SELECT file_path FROM `{tracking_table}`
        )
        SELECT gcs_files.file_path
        FROM gcs_files
        LEFT JOIN loaded_files USING (file_path)
        WHERE loaded_files.file_path IS NULL
        ORDER BY gcs_files.file_path
        """
        
        missing_results = client.query(missing_query)
        missing_paths = [row.file_path for row in missing_results]
        
        missing_files_list = f"{output_dir}/missing_files.txt"
        with open(missing_files_list, 'w') as f:
            for path in missing_paths:
                f.write(f"{path}\n")
        
        print(f"  Wrote {len(missing_paths)} missing file paths to {missing_files_list}")
    
    # Clean up staging table
    print(f"\nCleaning up staging table...")
    client.delete_table(staging_table, not_found_ok=True)
    
    # Write results
    results_dict = {
        "all_loaded": all_loaded,
        "total_files": total_files,
        "loaded_files": loaded_files,
        "missing_files": missing_files,
        "missing_files_list": missing_files_list,
    }
    
    results_file = f"{output_dir}/verification_results.json"
    with open(results_file, 'w') as f:
        json.dump(results_dict, f, indent=2)
    
    print(f"\nResults written to {results_file}")
    
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
    
    args = parser.parse_args()
    
    results = verify_all_loaded(
        project_id=args.project_id,
        dataset_name=args.dataset_name,
        gcs_files_list=args.gcs_files_list,
        output_dir=args.output_dir,
    )
    
    # Print summary
    print("\n" + "="*60)
    if results["all_loaded"]:
        print("✓ SUCCESS: All files have been loaded!")
    else:
        print(f"✗ INCOMPLETE: {results['missing_files']} files not yet loaded")
        print(f"  See: {results['missing_files_list']}")
    print("="*60)
    
    # Exit with error if not all loaded
    if not results["all_loaded"]:
        sys.exit(1)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
