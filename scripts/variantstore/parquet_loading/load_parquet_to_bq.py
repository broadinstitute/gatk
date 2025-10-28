#!/usr/bin/env python3
"""
Load Parquet files from GCS to BigQuery with fault tolerance and idempotency.
Uses deterministic job IDs and job state persistence to handle preemptions.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

from google.cloud import bigquery
from google.api_core import exceptions


def load_table_from_parquet_files(
    project_id,
    dataset_name,
    table_name,
    files_fofn,
    schema_path,
    pending_jobs_path,
    location=None,
    batch_size=10000,
):
    """
    Load Parquet files into a BigQuery table with fault tolerance.
    
    Args:
        project_id: BigQuery project ID
        dataset_name: BigQuery dataset name
        table_name: Target table name (e.g., 'vet_001')
        files_fofn: Path to file containing GCS URIs to load
        schema_path: Path to JSON schema file
        pending_jobs_path: Path to persist in-flight job state
        location: BigQuery location (auto-detected if None)
        batch_size: Number of files per load job
    
    Returns:
        Dictionary with load statistics
    """
    client = bigquery.Client(project=project_id)
    table_id = f"{project_id}.{dataset_name}.{table_name}"
    tracking_table_id = f"{project_id}.{dataset_name}.parquet_load_status"
    
    # Load schema
    with open(schema_path) as f:
        schema_json = json.load(f)
    schema = [bigquery.SchemaField.from_api_repr(field) for field in schema_json]
    
    # Fail fast if table was not pre-created
    try:
        table = client.get_table(table_id)
        location = location or table.location
        print(f"Target table: {table_id} in {location}")
    except exceptions.NotFound:
        print(f"ERROR: Table {table_id} does not exist. Tables must be pre-created.")
        sys.exit(1)
    
    # Load pending jobs state
    pending_jobs = _load_pending_jobs(pending_jobs_path)
    
    # Read files to load
    with open(files_fofn) as f:
        files = [line.strip() for line in f if line.strip()]
    
    if not files:
        print("No files to load")
        return {"files_loaded": 0, "rows_loaded": 0, "batches_processed": 0}
    
    print(f"Processing {len(files)} files in batches of {batch_size}")
    
    # Configure load job
    job_config = bigquery.LoadJobConfig(
        source_format=bigquery.SourceFormat.PARQUET,
        schema=schema,
        write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
    )
    
    total_rows = 0
    successful_batches = 0
    failed_batches = 0
    successful_loads = []
    
    # Load files in batches
    for i in range(0, len(files), batch_size):
        batch = files[i : i + batch_size]
        batch_num = i // batch_size + 1
        
        print(f"\nBatch {batch_num}/{(len(files) + batch_size - 1) // batch_size}: {len(batch)} files")
        
        try:
            job_id = _make_job_id(table_name, batch)
            
            if job_id in pending_jobs:
                # Job was previously submitted; check its status
                print(f"  Checking status of existing job {job_id}")
                try:
                    load_job = client.get_job(job_id=job_id, location=location)
                    print(f"  Resumed job in state: {load_job.state}")
                except exceptions.NotFound:
                    # Job never started or was cleaned up
                    print(f"  Job not found, will resubmit")
                    load_job = client.load_table_from_uri(
                        batch,
                        table_id,
                        job_config=job_config,
                        job_id=job_id,
                        location=location,
                    )
            else:
                # New job - submit it
                print(f"  Submitting new job {job_id}")
                load_job = client.load_table_from_uri(
                    batch,
                    table_id,
                    job_config=job_config,
                    job_id=job_id,
                    location=location,
                )
                pending_jobs[job_id] = {
                    "files": batch,
                    "location": location,
                    "table_name": table_name,
                }
                _store_pending_jobs(pending_jobs_path, pending_jobs)
            
            # Wait for completion
            load_job.result()
            
            # Check for errors
            if load_job.errors:
                print(f"  ERROR: Job completed with errors:")
                for error in load_job.errors:
                    print(f"    {error}")
                pending_jobs.pop(job_id, None)
                _store_pending_jobs(pending_jobs_path, pending_jobs)
                failed_batches += 1
                continue
            
            # Success!
            rows_loaded = load_job.output_rows or 0
            print(f"  ✓ Loaded {rows_loaded:,} rows from {len(batch)} files")
            total_rows += rows_loaded
            successful_batches += 1
            
            # Record successful loads in tracking table
            # For single-file batches, we know exact row count; for multi-file, set to None
            rows_per_file = rows_loaded if len(batch) == 1 else None
            for file_path in batch:
                successful_loads.append({
                    "file_path": file_path,
                    "table_name": table_name,
                    "file_size_bytes": None,  # Could populate via storage API
                    "load_job_id": load_job.job_id,
                    "rows_loaded": rows_per_file,
                })
            
            # Remove from pending jobs
            pending_jobs.pop(job_id, None)
            _store_pending_jobs(pending_jobs_path, pending_jobs)
            
        except Exception as exc:
            print(f"  ERROR: Exception loading batch: {exc}")
            failed_batches += 1
            # Continue with next batch rather than failing entire task
            continue
    
    # Insert tracking records
    if successful_loads:
        print(f"\nRecording {len(successful_loads)} files in tracking table...")
        insert_tracking_records(client, tracking_table_id, successful_loads)
    
    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY for {table_name}:")
    print(f"  Total files: {len(files)}")
    print(f"  Successful batches: {successful_batches}")
    print(f"  Failed batches: {failed_batches}")
    print(f"  Total rows loaded: {total_rows:,}")
    print(f"{'='*60}")
    
    return {
        "files_loaded": len(successful_loads),
        "rows_loaded": total_rows,
        "batches_processed": successful_batches,
        "batches_failed": failed_batches,
        "completion_status": "SUCCESS" if failed_batches == 0 else "PARTIAL",
    }


def insert_tracking_records(client, tracking_table_id, records):
    """
    Insert load tracking records, skipping duplicates.
    Batches records to avoid exceeding BigQuery request limits.
    """
    chunk_size = 1000  # Keep well below 16MB limit
    
    for i in range(0, len(records), chunk_size):
        chunk = records[i : i + chunk_size]
        
        merge_query = f"""
        MERGE `{tracking_table_id}` T
        USING UNNEST(@records) S
        ON T.file_path = S.file_path
        WHEN NOT MATCHED THEN
          INSERT (file_path, table_name, file_size_bytes, load_timestamp, load_job_id, rows_loaded)
          VALUES (S.file_path, S.table_name, S.file_size_bytes, CURRENT_TIMESTAMP(), S.load_job_id, S.rows_loaded)
        """
        
        job_config = bigquery.QueryJobConfig(
            query_parameters=[
                bigquery.ArrayQueryParameter(
                    "records",
                    "STRUCT",
                    [
                        {
                            "file_path": r["file_path"],
                            "table_name": r["table_name"],
                            "file_size_bytes": r["file_size_bytes"],
                            "load_job_id": r["load_job_id"],
                            "rows_loaded": r["rows_loaded"],
                        }
                        for r in chunk
                    ],
                )
            ]
        )
        
        query_job = client.query(merge_query, job_config=job_config)
        query_job.result()


def _make_job_id(table_name, batch):
    """Create deterministic job ID from table name and batch contents."""
    digest = hashlib.sha1("\n".join(sorted(batch)).encode("utf-8")).hexdigest()[:16]
    return f"load_{table_name}_{digest}"


def _load_pending_jobs(pending_jobs_path):
    """Load pending jobs state from disk."""
    path = Path(pending_jobs_path)
    if path.exists():
        try:
            return json.loads(path.read_text())
        except Exception as e:
            print(f"Warning: Could not load pending jobs: {e}")
            return {}
    return {}


def _store_pending_jobs(pending_jobs_path, pending_jobs):
    """Store pending jobs state to disk."""
    path = Path(pending_jobs_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(pending_jobs, indent=2))


def main():
    parser = argparse.ArgumentParser(
        description="Load Parquet files into BigQuery with fault tolerance"
    )
    parser.add_argument("--project-id", required=True, help="BigQuery project ID")
    parser.add_argument("--dataset-name", required=True, help="BigQuery dataset name")
    parser.add_argument("--table-name", required=True, help="Target table name")
    parser.add_argument("--files-fofn", required=True, help="File listing Parquet URIs to load")
    parser.add_argument("--schema-path", required=True, help="Path to JSON schema file")
    parser.add_argument(
        "--pending-jobs-path",
        required=True,
        help="Path to persist pending job state"
    )
    parser.add_argument("--location", help="BigQuery location (auto-detected if not specified)")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        help="Number of files per load job (default: 10000)"
    )
    parser.add_argument("--output-stats", help="Path to write statistics JSON")
    
    args = parser.parse_args()
    
    stats = load_table_from_parquet_files(
        project_id=args.project_id,
        dataset_name=args.dataset_name,
        table_name=args.table_name,
        files_fofn=args.files_fofn,
        schema_path=args.schema_path,
        pending_jobs_path=args.pending_jobs_path,
        location=args.location,
        batch_size=args.batch_size,
    )
    
    if args.output_stats:
        with open(args.output_stats, 'w') as f:
            json.dump(stats, f, indent=2)
    
    # Exit with error if there were failures
    if stats["completion_status"] != "SUCCESS":
        sys.exit(1)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
