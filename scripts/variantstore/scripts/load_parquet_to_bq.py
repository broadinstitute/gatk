#!/usr/bin/env python3
"""
Load Parquet files from GCS to BigQuery with fault tolerance and idempotency.
Uses deterministic job IDs to handle preemptions; no local state file is needed.
A job that fails permanently is superseded by a fresh generation of its job ID so
that a task retry is not stuck re-fetching the same dead job (see
_submit_resumable_load_job).
"""

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path


def load_table_from_parquet_files(
    project_id,
    dataset_name,
    table_name,
    files_fofn,
    schema_path,
    location=None,
    batch_size=10000,
):
    """
    Load Parquet files into a BigQuery table with fault tolerance.
    
    Args:
        project_id: BigQuery project ID
        dataset_name: BigQuery dataset name
        table_name: Target table name (e.g., 'vet_001'). If None, extracted from files_fofn filename.
        files_fofn: Path to file containing GCS URIs to load
        schema_path: Path to JSON schema file (optional, can be None for schema autodetection)
        location: BigQuery location (auto-detected if None)
        batch_size: Number of files per load job
    
    Returns:
        Dictionary with load statistics
    """
    from google.cloud import bigquery
    from google.api_core import exceptions

    client = bigquery.Client(project=project_id)
    
    # Extract table name from FOFN filename if not provided
    if table_name is None:
        fofn_path = Path(files_fofn)
        # Remove .fofn extension to get table name
        table_name = fofn_path.stem
        print(f"Extracted table name from FOFN: {table_name}")
    
    table_id = f"{project_id}.{dataset_name}.{table_name}"

    # Load schema if provided, otherwise use autodetection
    schema = None
    if schema_path:
        with open(schema_path) as f:
            schema_json = json.load(f)
        schema = [bigquery.SchemaField.from_api_repr(field) for field in schema_json]
    else:
        print("No schema provided; using autodetection from Parquet files")
    
    # Fail fast if table was not pre-created
    try:
        table = client.get_table(table_id)
        location = location or table.location
        print(f"Target table: {table_id} in {location}")
    except exceptions.NotFound:
        print(f"ERROR: Table {table_id} does not exist. Tables must be pre-created.")
        sys.exit(1)
    
    # Read files to load
    with open(files_fofn) as f:
        files = [line.strip() for line in f if line.strip()]
    
    if not files:
        print("No files to load")
        return {
            "files_loaded": 0,
            "rows_loaded": 0,
            "batches_processed": 0,
            "batches_failed": 0,
            "completion_status": "SUCCESS",
        }
    
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
    successful_files_count = 0
    
    # Load files in batches
    for i in range(0, len(files), batch_size):
        batch = files[i : i + batch_size]
        batch_num = i // batch_size + 1
        
        print(f"\nBatch {batch_num}/{(len(files) + batch_size - 1) // batch_size}: {len(batch)} files")

        job_id = None  # ensure it is always defined for the except handler
        try:
            # Submit the job, resuming an in-flight or already-succeeded job on a
            # deterministic job-ID collision and superseding a permanently-failed one.
            load_job, job_id = _submit_resumable_load_job(
                client, project_id, dataset_name, table_name, batch,
                table_id, job_config, location,
            )

            # Wait for completion with retry logic
            load_job = _wait_for_job_with_retry(load_job)
            
            # Check for errors
            if load_job.errors:
                error_msg = "; ".join(
                    e.get("message", str(e)) for e in load_job.errors
                )
                print(f"  ERROR: Job completed with errors:")
                for error in load_job.errors:
                    print(f"    {error}")
                _log_batch_result(batch, table_name, load_job.job_id, error=error_msg)
                failed_batches += 1
                continue
            
            # Success!
            rows_loaded = load_job.output_rows or 0
            print(f"  ✓ Loaded {rows_loaded:,} rows from {len(batch)} files")
            total_rows += rows_loaded
            successful_batches += 1
            successful_files_count += len(batch)
            _log_batch_result(batch, table_name, load_job.job_id)
            
        except Exception as exc:
            print(f"  ERROR: Exception loading batch: {exc}")
            _log_batch_result(batch, table_name, job_id, error=str(exc))
            failed_batches += 1
            # Continue with next batch rather than failing entire task
            continue
    
    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY for {table_name}:")
    print(f"  Total files: {len(files)}")
    print(f"  Successful batches: {successful_batches}")
    print(f"  Failed batches: {failed_batches}")
    print(f"  Total rows loaded: {total_rows:,}")
    print(f"{'='*60}")
    
    return {
        "files_loaded": successful_files_count,
        "rows_loaded": total_rows,
        "batches_processed": successful_batches,
        "batches_failed": failed_batches,
        "completion_status": "SUCCESS" if failed_batches == 0 else "PARTIAL",
    }



def _make_job_id(project_id, dataset_name, table_name, batch, generation=0):
    """Create deterministic job ID from project, dataset, table name and batch contents.

    The optional generation lets a permanently-failed job be superseded by a fresh
    job ID on a later task attempt (see _submit_resumable_load_job).  Generation 0
    (the default) is unsuffixed so IDs minted before this parameter existed still
    resolve identically; generation N>=1 appends an "_rN" suffix.
    """
    digest = hashlib.sha1(
        "\n".join([project_id, dataset_name] + sorted(batch)).encode("utf-8")
    ).hexdigest()[:16]
    job_id = f"load_{table_name}_{digest}"
    if generation:
        job_id = f"{job_id}_r{generation}"
    return job_id




def _extract_sample_id_from_file_path(file_path, table_name):
    """
    Extract the sample_id integer from a GCS Parquet file path given the resolved table name.

    For superpartitioned tables the filename encodes both the superpartition and the sample_id:
        gs://bucket/vet/vet_001_42_input_vcf_0_SAMPLE.parquet  -> "42"   (table_name="vet_001")

    For regular tables the filename encodes only the sample_id:
        gs://bucket/sample_chromosome_ploidy/sample_chromosome_ploidy_42_SAMPLE.parquet -> "42"

    Returns the sample_id as a string, or None if the pattern is not found.
    """
    pattern = rf'/{re.escape(table_name)}_([0-9]+)(?:[_.]|$)'
    match = re.search(pattern, file_path)
    return match.group(1) if match else None


def _log_batch_result(files, table_name, job_id, error=None):
    """
    Print a structured log block for a completed (successful or failed) batch.

    Format:
        BATCH WRITE {SUCCEEDED|FAILED} for job <job_id>, table <table_name>.
        [ERROR: <error message>]          <- only present when the batch failed
        <sample_id>\t<file_path>          <- one line per file
    """
    status = "FAILED" if error is not None else "SUCCEEDED"
    job_id_str = str(job_id) if job_id is not None else ""
    print(f"BATCH WRITE {status} for job {job_id_str}, table {table_name}.")
    if error is not None:
        print(f"ERROR: {error}")
    for file_path in files:
        sample_id = _extract_sample_id_from_file_path(file_path, table_name)
        sample_id_str = sample_id if sample_id is not None else ""
        print(f"{sample_id_str}\t{file_path}")


def _is_permanently_failed_job(load_job):
    """
    True if load_job is an existing job that finished in DONE state with an
    error_result.

    A BigQuery load job is atomic: one that terminates with an error_result loaded
    zero rows, so it is safe to supersede with a fresh job ID rather than re-fetching
    it forever.  A freshly submitted job is PENDING/RUNNING with no error_result and
    is therefore never considered permanently failed here.
    """
    return load_job.state == "DONE" and load_job.error_result is not None


def _submit_resumable_load_job(
    client, project_id, dataset_name, table_name, batch,
    table_id, job_config, location, max_generations=100,
):
    """
    Submit a load job for a batch, resuming on a deterministic job-ID collision and
    superseding a pre-existing permanently-failed job (one we collide with) with a
    fresh generation of the ID.

    A BigQuery job ID cannot be reused, so once a job reaches DONE with an
    error_result (a malformed Parquet file, a schema mismatch) its ID is burned:
    re-submitting the same ID only re-fetches the dead job, whose result() re-raises
    the original error.  Because the base ID is a pure function of the batch
    contents, no task retry could otherwise make the batch succeed.  Walking a
    generation suffix mints a genuinely new ID on the next attempt while still
    resuming a job that is merely PENDING/RUNNING or already succeeded -- which must
    never be duplicated, since the load is WRITE_APPEND.

    Returns:
        (load_job, job_id) for the generation that is in flight, already succeeded,
        or freshly submitted.

    Raises:
        Exception: if the job submission fails for a non-Conflict reason, or if every
        generation up to max_generations is an existing permanently-failed job.
    """
    for generation in range(max_generations):
        job_id = _make_job_id(project_id, dataset_name, table_name, batch, generation)
        print(f"  Submitting job {job_id}")
        load_job, resumed_existing = _submit_load_job_with_retry(
            client, batch, table_id, job_config, job_id, location
        )
        # Only supersede a job we *collided with*: a pre-existing job that reached DONE
        # with an error_result loaded zero rows and its ID is burned.  A job we just
        # submitted is returned as-is even if it fast-failed (e.g. a schema rejection
        # already terminal in the insert response), so the caller records the batch
        # failure and a later task attempt supersedes it.  Superseding a fresh failure
        # here would burn every generation in a single attempt, since each new
        # generation fast-fails identically -- leaving all later retries exhausted.
        if resumed_existing and _is_permanently_failed_job(load_job):
            print(
                f"  Job {job_id} already exists and failed permanently "
                f"({load_job.error_result}); superseding it with a fresh job ID"
            )
            continue
        return load_job, job_id

    raise Exception(
        f"Exhausted {max_generations} job-ID generations for table {table_name}; "
        f"every generation is an existing permanently-failed job"
    )


def _submit_load_job_with_retry(client, batch, table_id, job_config, job_id, location, max_retries=3):
    """
    Submit a BigQuery load job with exponential backoff for quota/rate limit errors.
    
    Args:
        client: BigQuery client
        batch: List of GCS URIs to load
        table_id: Target table ID
        job_config: LoadJobConfig
        job_id: Deterministic job ID
        location: BigQuery location
        max_retries: Maximum number of retry attempts (default: 3)

    Returns:
        (load_job, resumed_existing).  resumed_existing is True only when the job ID
        already existed and the job was fetched via get_job() (a Conflict); it is
        False for a job we freshly submitted.  Callers use this to distinguish a
        pre-existing job -- which may be safely superseded if it failed permanently --
        from one we just created, which must be returned as-is even if it is already
        terminal (see _submit_resumable_load_job).

    Raises:
        Exception: If all retries are exhausted or a non-retryable error occurs
    """
    from google.api_core import exceptions

    retry_delays = [30, 60, 120]  # Exponential backoff: 30s, 60s, 120s

    for attempt in range(max_retries + 1):
        try:
            load_job = client.load_table_from_uri(
                batch,
                table_id,
                job_config=job_config,
                job_id=job_id,
                location=location,
            )
            return load_job, False

        except (exceptions.TooManyRequests, 
                exceptions.ServiceUnavailable,
                exceptions.InternalServerError) as e:
            # These are retryable quota/availability errors
            if attempt < max_retries:
                delay = retry_delays[attempt]
                print(f"  Quota/rate limit error: {e}")
                print(f"  Retrying in {delay} seconds (attempt {attempt + 1}/{max_retries})...")
                time.sleep(delay)
            else:
                print(f"  ERROR: Max retries ({max_retries}) exhausted for quota errors")
                raise
        
        except exceptions.Conflict as e:
            # Job ID already exists - this is actually OK.  Fetch the existing job
            # rather than treating it as an error, and flag it as pre-existing so the
            # caller may supersede it if it turns out to have failed permanently.
            print(f"  Job {job_id} already exists, fetching existing job")
            return client.get_job(job_id=job_id, location=location), True
        
        except Exception as e:
            # Non-retryable error, fail immediately
            print(f"  Non-retryable error submitting job: {e}")
            raise
    
    # Should never reach here, but just in case
    raise Exception("Unexpected error in retry logic")


def _wait_for_job_with_retry(load_job, max_retries=3):
    """
    Wait for a BigQuery load job to complete with retry logic for transient errors.
    
    Args:
        load_job: LoadJob object
        max_retries: Maximum number of retry attempts (default: 3)
    
    Returns:
        Completed LoadJob object
    
    Raises:
        Exception: If all retries are exhausted or a non-retryable error occurs
    """
    from google.api_core import exceptions

    retry_delays = [30, 60, 120]  # Exponential backoff
    
    for attempt in range(max_retries + 1):
        try:
            # Wait for job completion
            result = load_job.result()
            return load_job
            
        except (exceptions.TooManyRequests,
                exceptions.ServiceUnavailable, 
                exceptions.InternalServerError) as e:
            # Transient error while waiting - refresh job state and retry
            if attempt < max_retries:
                delay = retry_delays[attempt]
                print(f"  Transient error while waiting for job: {e}")
                print(f"  Retrying in {delay} seconds (attempt {attempt + 1}/{max_retries})...")
                time.sleep(delay)
                # Refresh the job object
                load_job.reload()
            else:
                print(f"  ERROR: Max retries ({max_retries}) exhausted while waiting for job")
                raise
        
        except Exception as e:
            # Non-retryable error (e.g., actual job failure)
            print(f"  Job failed with non-retryable error: {e}")
            raise
    
    raise Exception("Unexpected error in wait retry logic")


def _execute_query_with_retry(client, query, job_config=None, max_retries=3):
    """
    Execute a BigQuery query with exponential backoff for quota/rate limit errors.
    
    Args:
        client: BigQuery client
        query: SQL query string
        job_config: QueryJobConfig (optional)
        max_retries: Maximum number of retry attempts (default: 3)
    
    Returns:
        Query results
    
    Raises:
        Exception: If all retries are exhausted or a non-retryable error occurs
    """
    from google.api_core import exceptions

    retry_delays = [30, 60, 120]  # Exponential backoff
    
    for attempt in range(max_retries + 1):
        try:
            query_job = client.query(query, job_config=job_config)
            result = query_job.result()
            return result
            
        except (exceptions.TooManyRequests,
                exceptions.ServiceUnavailable,
                exceptions.InternalServerError) as e:
            # Retryable quota/availability errors
            if attempt < max_retries:
                delay = retry_delays[attempt]
                print(f"  Query quota/rate limit error: {e}")
                print(f"  Retrying in {delay} seconds (attempt {attempt + 1}/{max_retries})...")
                time.sleep(delay)
            else:
                print(f"  ERROR: Max retries ({max_retries}) exhausted for query")
                raise
        
        except Exception as e:
            # Non-retryable error
            print(f"  Non-retryable query error: {e}")
            raise
    
    raise Exception("Unexpected error in query retry logic")



def main():
    parser = argparse.ArgumentParser(
        description="Load Parquet files into BigQuery with fault tolerance"
    )
    parser.add_argument("--project-id", required=True, help="BigQuery project ID")
    parser.add_argument("--dataset-name", required=True, help="BigQuery dataset name")
    parser.add_argument("--table-name", help="Target table name (if not provided, extracted from FOFN filename)")
    parser.add_argument("--files-fofn", required=True, help="File listing Parquet URIs to load")
    parser.add_argument("--schema-path", help="Path to JSON schema file (optional, autodetect if not provided)")
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
