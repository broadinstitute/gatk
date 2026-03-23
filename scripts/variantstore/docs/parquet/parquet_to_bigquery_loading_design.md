# Design Document: Robust Parquet to BigQuery Loading

## Overview

This document outlines the design for a WDL-based system to load hundreds of thousands of Parquet files from Google Cloud Storage (GCS) into BigQuery tables in a robust, idempotent manner that can handle interruptions and partial completions.

## Problem Statement

After `CreateVariantIngestFiles` writes Parquet files to GCS in a structured directory hierarchy, we need to:

1. Load all Parquet files into their corresponding BigQuery tables
2. Ensure each file is loaded exactly once (idempotency)
3. Handle cases where tables may already contain data from previous runs
4. Be resilient to task interruptions and preemptions
5. Provide definitive completion status when all files are loaded

## GCS Directory Structure

```
${OUTPUT_GCS_DIR}/
├── vet/
│   ├── 000/
│   │   ├── vet_000_sample1.parquet
│   │   ├── vet_000_sample2.parquet
│   │   └── ...
│   ├── 001/
│   │   ├── vet_001_sample1.parquet
│   │   └── ...
│   └── ...
└── ref_ranges/
  ├── 000/
  │   ├── ref_000_sample1.parquet
  │   └── ...
  └── ...
```

**Mapping Rule**: `${OUTPUT_GCS_DIR}/{table_prefix}/{table_number}/` → BigQuery table `{table_prefix}_{table_number}` (e.g., `vet/001/…` → `vet_001`, `ref_ranges/042/…` → `ref_042`)

**Constraints**:
- Vet Parquet files are ≈191 MB each; ref-range Parquet files are ≈160 MB each
- Maximum 4,000 files per table partition (system-imposed)
- Total files could be hundreds of thousands
- Files may be written incrementally across multiple ingestion runs
- BigQuery dataset and destination tables (including schemas) are provisioned explicitly before this workflow runs

## Design Principles

1. **Idempotency**: Safe to run multiple times; already-loaded files are skipped
2. **Fault Tolerance**: Interruptions don't corrupt state; can resume from where it left off
3. **Atomicity per File**: Each file load is tracked individually
4. **Verifiability**: Clear way to determine if all files have been loaded
5. **Parallelization**: Load multiple tables/files concurrently for performance

## Proposed Solution

### Architecture Overview

The solution consists of three main WDL tasks:

1. **DiscoverParquetFiles**: Inventory all Parquet files and group by target table
2. **LoadParquetFilesToBQ**: Load files for a single table with tracking
3. **VerifyAllLoaded**: Confirm all files were successfully loaded

### Idempotency Mechanism

Idempotency is determined by querying BigQuery directly — no separate tracking table is created or maintained.

**For superpartitioned tables** (`vet_%`, `ref_ranges_%`), `INFORMATION_SCHEMA.PARTITIONS` is queried for non-empty partitions whose `partition_id` parses as a valid integer `sample_id`:

```sql
SELECT table_name, SAFE_CAST(partition_id AS INT64) AS sample_id
FROM `{project}.{dataset}.INFORMATION_SCHEMA.PARTITIONS`
WHERE
    REGEXP_CONTAINS(table_name, "^vet_[0-9]+$|^ref_ranges_[0-9]+$")
    AND total_logical_bytes > 0
    AND partition_id NOT LIKE "__%"
    AND SAFE_CAST(partition_id AS INT64) IS NOT NULL
```

**For regular tables** (`sample_chromosome_ploidy`), the table itself is queried for distinct `sample_id` values:

```sql
SELECT DISTINCT "sample_chromosome_ploidy" AS table_name, sample_id
FROM `{project}.{dataset}.sample_chromosome_ploidy`
```

Any `(table_name, sample_id)` pair already present in BigQuery is skipped without attempting a reload. This is more reliable than a tracking table because:
- It reflects the actual committed state of BigQuery data (not a secondary record of intent)
- It guards against the case where a tracking table claims data is loaded that actually isn't
- It eliminates high-frequency MERGE DML writes that were causing quota pressure at scale

### Detailed Task Design

#### Task 1: DiscoverParquetFiles

**Purpose**: Scan GCS to inventory all Parquet files and group them by target BigQuery table.

**Inputs**:
- `output_gcs_dir`: Base GCS path
- `project_id`: BigQuery project
- `dataset_name`: BigQuery dataset
- `table_prefixes`: Array of prefixes to scan (e.g., ["vet", "ref_ranges"])

**Process**:
1. Use `gcloud storage ls --recursive` to list all `.parquet` files
2. Parse file paths to extract table numbers and sample IDs
3. Group files by table (e.g., `vet_000`, `vet_001`, etc.)
4. Query BigQuery (`INFORMATION_SCHEMA.PARTITIONS` and `sample_chromosome_ploidy`) to identify which `(table_name, sample_id)` pairs are already loaded
5. Create file-of-file-names (FOFN) for each table containing only files whose pair is not yet loaded

**Outputs**:
- Array of table names (e.g., `["vet_000", "vet_001", ..., "ref_ranges_000", ...]`)
- Array of FOFNs, one per table, containing GCS paths of files to load
- Summary statistics (total files, already loaded, remaining)

**Implementation Notes**:
```bash
# List all objects once, then filter for parquet paths
gcloud storage ls --recursive gs://${OUTPUT_GCS_DIR}/ > all_objects.txt
grep '\.parquet$' all_objects.txt > all_files.txt

# Parse and group by table; already-loaded (table_name, sample_id) pairs are
# detected via INFORMATION_SCHEMA.PARTITIONS and sample_chromosome_ploidy queries
python3 parse_and_group_files.py \
  --input all_files.txt \
  --output-dir grouped_files/ \
  --project-id ${PROJECT_ID} \
  --dataset ${DATASET} \
  --superpartitioned-table-prefixes vet ref_ranges \
  --regular-table-prefixes sample_chromosome_ploidy
```

**Volatility**: Mark as `volatile: true` since we always want fresh GCS and BigQuery state.

---

#### Task 2: LoadParquetFilesToBQ

**Purpose**: Load all unloaded Parquet files for a single BigQuery table.

**Inputs**:
- `project_id`: BigQuery project
- `dataset_name`: BigQuery dataset
- `table_name`: Target table (e.g., `vet_000`)
- `files_to_load`: File containing GCS paths to load (FOFN)
- `schema_path`: Path to the canonical table schema JSON (already provisioned during table creation)
- `pending_jobs_path`: Path to a JSON file used to persist in-flight BigQuery load jobs (stored on the Cromwell VM or in GCS)
- `billing_project_id`: Optional billing override

**Notes**:
- With 4,000 files per table at ~191 MB (vet) or ~160 MB (ref ranges), the worst-case payload per table is ~764 GB, well within BigQuery's 15 TB per-load limit.
- We still permit optional batching to keep retries fast and to respect project-level load quotas.
- Tables and schemas must already exist; the task should fail fast if `table_id` is missing rather than attempting to create it.
- Every batch is launched with a deterministic BigQuery `job_id` derived from the table name and batch contents. The `pending_jobs_path` tracks submitted jobs so that, after a preemption, the task can call `client.get_job(job_id)` to determine whether the work finished before resubmitting.

**Process**:

```python
#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path

from google.cloud import bigquery


def load_table_from_parquet_files(
  project_id,
  dataset_name,
  table_name,
  files_fofn,
  schema_path,
  pending_jobs_path,
  location=None,
):
  client = bigquery.Client(project=project_id)
  table_id = f"{project_id}.{dataset_name}.{table_name}"
  schema = client.schema_from_json(schema_path)

  # Fail fast if table was not pre-created
  table = client.get_table(table_id)
  location = location or table.location

  pending_jobs = _load_pending_jobs(pending_jobs_path)

  # Read files to load
  with open(files_fofn) as f:
    files = [line.strip() for line in f if line.strip()]

  if not files:
    print("No files to load")
    return 0

  # Configure load job using explicit schema
  job_config = bigquery.LoadJobConfig(
    source_format=bigquery.SourceFormat.PARQUET,
    schema=schema,
    write_disposition=bigquery.WriteDisposition.WRITE_APPEND,
  )

  total_rows = 0
  successful_loads = []

  # Load files in batches (BigQuery supports up to 10,000 URIs per job)
  batch_size = 10000  # Adjust (e.g., 500-1000) if throttling or faster retries are desired
  for i in range(0, len(files), batch_size):
    batch = files[i : i + batch_size]

    print(f"Loading batch {i // batch_size + 1}: {len(batch)} files")

    try:
      job_id = _make_job_id(table_name, batch)

      if job_id in pending_jobs:
        # Job was previously submitted; refresh its status
        load_job = client.get_job(job_id=job_id, location=location)
        print(f"Resuming existing load job {job_id} in state {load_job.state}")
      else:
        load_job = client.load_table_from_uri(
          batch,
          table_id,
          job_config=job_config,
          job_id=job_id,
          location=location,
        )
        pending_jobs[job_id] = {"files": batch, "location": location}
        _store_pending_jobs(pending_jobs_path, pending_jobs)

      load_job.result()  # Wait for completion

      if load_job.errors:
        print(f"Job {job_id} completed with errors: {load_job.errors}")
        pending_jobs.pop(job_id, None)
        _store_pending_jobs(pending_jobs_path, pending_jobs)
        continue

      print(f"Loaded {load_job.output_rows} rows from {len(batch)} files")
      total_rows += load_job.output_rows

      # Record successful loads in tracking table
      rows_per_file = load_job.output_rows if len(batch) == 1 else None
      for file_path in batch:
        successful_loads.append(
          {
            "file_path": file_path,
            "table_name": table_name,
            "file_size_bytes": None,  # Optional: populate via GCS metadata lookup
            "load_job_id": load_job.job_id,
            "rows_loaded": rows_per_file,
          }
        )

      pending_jobs.pop(job_id, None)
      _store_pending_jobs(pending_jobs_path, pending_jobs)

    except Exception as exc:  # pragma: no cover - design sample code
      print(f"Error loading batch: {exc}")
      # Continue with next batch rather than failing entire task
      continue

  # Insert tracking records
  if successful_loads:
    insert_tracking_records(client, tracking_table_id, successful_loads)

  print(f"Total rows loaded: {total_rows}")
  return len(successful_loads)


def insert_tracking_records(client, tracking_table_id, records):
  """Insert load tracking records, skipping duplicates."""

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
      bigquery.ArrayQueryParameter("records", "STRUCT", records)
    ]
  )

  query_job = client.query(merge_query, job_config=job_config)
  query_job.result()


def _make_job_id(table_name, batch):
  digest = hashlib.sha1("\n".join(sorted(batch)).encode("utf-8")).hexdigest()[:16]
  return f"load_{table_name}_{digest}"


def _load_pending_jobs(pending_jobs_path):
  path = Path(pending_jobs_path)
  if path.exists():
    return json.loads(path.read_text())
  return {}


def _store_pending_jobs(pending_jobs_path, pending_jobs):
  Path(pending_jobs_path).write_text(json.dumps(pending_jobs))
```

> **Batching note**: Keep each `records` payload comfortably below BigQuery's 16 MB request limit (e.g., chunks of ≤1,000 records). For very large batches, write tracking records to a temporary table and issue a single `MERGE` from that table instead of UNNESTing parameters.

> **Metadata capture**: Populate `file_size_bytes` via a `storage.Client().get_blob()` lookup when the additional API call is acceptable. When batches contain multiple files, `rows_loaded` is left null; derive aggregate row counts by querying BigQuery tables directly.

**Outputs**:
- `files_loaded_count`: Number of files successfully loaded
- `rows_loaded`: Total rows loaded
- `completion_status`: "SUCCESS" or "PARTIAL"

**Error Handling**:
- If a batch fails, log it but continue with remaining batches
- Return partial success status if some batches fail
- Failed files remain untracked and will be retried in next run
- Deterministic job IDs and the persisted `pending_jobs_path` allow resumed tasks to query `client.get_job()` before deciding whether to retry or skip a batch, eliminating duplicate writes after preemptions

**Parallelization**: This task is scattered across all tables.

**Preemptible**: Set to 3-5 to save costs; task is fully resumable.

---

#### Task 3: VerifyAllLoaded

**Purpose**: Confirm that all discovered files have been successfully loaded.

**Inputs**:
- `project_id`: BigQuery project
- `dataset_name`: BigQuery dataset
- `output_gcs_dir`: Base GCS path
- `table_prefixes`: Array of prefixes scanned
- `load_task_outputs`: Array from LoadParquetFilesToBQ tasks (for dependency)

**Process**:

Parse each GCS file path to extract its `(table_name, sample_id)` pair, then query BigQuery to find which pairs are already loaded. Files whose pair is not found in BigQuery are reported as missing.

```sql
-- For superpartitioned tables (vet/ref_ranges): check INFORMATION_SCHEMA.PARTITIONS
SELECT table_name, SAFE_CAST(partition_id AS INT64) AS sample_id
FROM `{project}.{dataset}.INFORMATION_SCHEMA.PARTITIONS`
WHERE
    REGEXP_CONTAINS(table_name, "^vet_[0-9]+$|^ref_ranges_[0-9]+$")
    AND total_logical_bytes > 0
    AND partition_id NOT LIKE "__%"
    AND SAFE_CAST(partition_id AS INT64) IS NOT NULL

UNION ALL

-- For regular tables (sample_chromosome_ploidy): query the table directly
SELECT DISTINCT "sample_chromosome_ploidy" AS table_name, sample_id
FROM `{project}.{dataset}.sample_chromosome_ploidy`;
```

GCS files whose `(table_name, sample_id)` pair is absent from the above result set are reported as missing.

**Outputs**:
- `all_loaded`: Boolean (true if missing_files == 0)
- `total_files`: Total count
- `loaded_files`: Loaded count
- `missing_files`: Unloaded count
- `missing_files_list`: File with paths of any missing files

**Failure Condition**: Fail the task if `missing_files > 0` after all load tasks complete.

---

### WDL Workflow Structure

```wdl
workflow LoadParquetsToBigQuery {
  input {
    String output_gcs_dir
    String project_id
    String dataset_name
    Array[String] table_prefixes = ["vet", "ref_ranges"]
    String? billing_project_id
    Int max_parallel_loads = 50
  }
  
  # Ensure tracking table exists
  call CreateTrackingTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name
  }
  
  # Discover all files and group by table
  call DiscoverParquetFiles {
    input:
      output_gcs_dir = output_gcs_dir,
      project_id = project_id,
      dataset_name = dataset_name,
      table_prefixes = table_prefixes,
      tracking_table_ready = CreateTrackingTable.done
  }
  
  # Load each table in parallel
  scatter (i in range(length(DiscoverParquetFiles.table_names))) {
    call LoadParquetFilesToBQ {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        table_name = DiscoverParquetFiles.table_names[i],
        files_to_load = DiscoverParquetFiles.file_fofns[i],
        billing_project_id = billing_project_id
    }
  }
  
  # Verify everything loaded successfully
  call VerifyAllLoaded {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      output_gcs_dir = output_gcs_dir,
      table_prefixes = table_prefixes,
      load_outputs = LoadParquetFilesToBQ.completion_status
  }
  
  output {
    Boolean all_loaded = VerifyAllLoaded.all_loaded
    Int total_files = VerifyAllLoaded.total_files
    Int loaded_files = VerifyAllLoaded.loaded_files
    Array[String] tables_loaded = DiscoverParquetFiles.table_names
    File? missing_files = VerifyAllLoaded.missing_files_list
  }
}
```

---

## Alternative: Batch Files Instead of Batch Tables

If the number of files per table is very large (approaching 4,000), loading all files for a table in a single task might take too long or risk preemption.

**Alternative Approach**: Batch files within each table into smaller groups (e.g., 500 files per load task).

```wdl
scatter (i in range(length(DiscoverParquetFiles.table_names))) {
  # Further split large tables into batches
  call CreateFileBatches {
    input:
      files_fofn = DiscoverParquetFiles.file_fofns[i],
      batch_size = 500
  }
  
  scatter (batch_fofn in CreateFileBatches.batch_fofns) {
    call LoadParquetFilesToBQ {
      input:
        table_name = DiscoverParquetFiles.table_names[i],
        files_to_load = batch_fofn,
        # ...
    }
  }
}
```

This provides more granular parallelization and better resilience to preemptions.

---

## Integration with GvsImportGenomes.wdl

### Option 1: Separate Workflow (Recommended)

Keep `LoadParquetsToBigQuery` as a separate workflow that runs after `GvsImportGenomes` completes.

**Advantages**:
- Clean separation of concerns
- Can re-run loading independently if needed
- Easier to test and debug
- Can handle multiple `GvsImportGenomes` runs aggregating to same tables

**Usage**:
```bash
# Run import (writes Parquet files)
cromwell run GvsImportGenomes.wdl -i import_inputs.json

# Run load (loads Parquet files to BQ)
cromwell run LoadParquetsToBigQuery.wdl -i load_inputs.json
```

### Option 2: Integrated Workflow

Add the loading tasks directly into `GvsImportGenomes.wdl` after the `LoadData` scatter completes.

**Advantages**:
- Single workflow execution
- Guaranteed ordering

**Disadvantages**:
- Tighter coupling
- Harder to re-run just the loading portion
- Workflow becomes more complex

**Implementation**: Add after the `LoadData` scatter:

```wdl
# After LoadData scatter...
call CreateTrackingTable { ... }
call DiscoverParquetFiles {
  input:
    load_data_done = LoadData.done,  # Dependency
    # ...
}
# ... rest of loading tasks
```

---

## Handling Edge Cases

### 1. Files Added During Load

**Scenario**: `GvsImportGenomes` is still writing files while loading is in progress.

**Solution**: 
- `DiscoverParquetFiles` is marked `volatile: true` and only captures files at scan time
- Run `VerifyAllLoaded` after all data writing is complete
- If needed, run the entire loading workflow again (idempotent)

### 2. Duplicate File Names

**Scenario**: Same filename written multiple times (e.g., re-ingesting a sample).

**Solution**:
- File paths are unique (include sample name and table partition)
- If truly duplicate, tracking table prevents reloading (file_path is primary key)
- BigQuery's `WRITE_APPEND` adds rows; use `WRITE_TRUNCATE` if table should be replaced

### 3. Task Preemption During Load

**Scenario**: `LoadParquetFilesToBQ` is preempted mid-execution.

**Solution**:
- Files already loaded are recorded in tracking table
- On retry, `DiscoverParquetFiles` excludes already-loaded files
- Task continues from where it left off

### 4. Schema Evolution

**Scenario**: Parquet files have different schemas over time.

**Solution**:
- Store a canonical schema JSON for each table and pass it explicitly to every load job
- BigQuery will raise schema mismatch errors immediately; capture and surface them so producers can correct upstream data
- If intentional schema changes are required, update the table definition and schema JSON first, then relaunch the loader

### 5. Very Large Number of Tables

**Scenario**: Hundreds of table partitions (vet_000 through vet_999, etc.).

**Solution**:
- Cromwell handles large scatters well
- Consider using sub-workflows if scatter exceeds thousands of elements
- Monitor BigQuery quota for concurrent load jobs

---

## Performance Considerations

### BigQuery Load Job Limits

- **Concurrent load jobs per project/region**: Soft limit ~100 (plan to throttle or request quota increases)
- **Concurrent load jobs per table**: Subject to the same project-level cap
- **Files per load job**: 10,000 (handled by batching)

### Optimization Strategies

1. **Batch Files per Load Job**: Load 500-10,000 files per BigQuery load job
2. **Parallel Table Loads**: Load different tables concurrently (controlled by WDL scatter) while respecting the ~100 job cap; throttle or stagger launches if quota errors appear
3. **Preemptible VMs**: Use for all tasks since they're resumable (save ~70% cost)
4. **Disk Space**: Minimal; only need to store FOFNs and tracking queries

### Cost Estimation

- **BigQuery Load**: Free (loads from GCS are free)
- **BigQuery Storage**: Standard rates for data and tracking table
- **Compute**: Minimal; mostly `gcloud` commands and Python scripts
- **GCS Operations**: List operations (minimal cost)

**Example**: Loading 100,000 files across 25 tables:
- Discovery: ~1 minute, $0.01
- Load tasks: ~25 concurrent tasks × 10 minutes = $0.50 (with preemptibles)
- Verification: ~1 minute, $0.01
- **Total**: ~$0.52 + BigQuery storage

---

## Testing Strategy

### Unit Tests

1. **Test File Grouping Logic**: Verify correct table name and sample ID extraction from paths
2. **Test BigQuery Query Construction**: Ensure already-loaded detection queries are correct
3. **Test Batch Creation**: Verify files are split into appropriate batch sizes

### Integration Tests

1. **Small Scale Test**: 10 files across 2 tables
2. **Medium Scale Test**: 1,000 files across 10 tables
3. **Preemption Test**: Kill a load task mid-execution and verify recovery
4. **Duplicate Run Test**: Run workflow twice and verify no duplicate loads

### Validation Queries

```sql
-- Overall progress per table (via partition metadata)
SELECT
  table_name,
  COUNT(*) as partitions_loaded,
  SUM(total_logical_bytes) as total_bytes
FROM `{dataset}.INFORMATION_SCHEMA.PARTITIONS`
WHERE REGEXP_CONTAINS(table_name, "^vet_[0-9]+$|^ref_ranges_[0-9]+$")
  AND total_logical_bytes > 0
  AND partition_id NOT LIKE "__%"
GROUP BY table_name
ORDER BY table_name;

-- Ploidy samples loaded
SELECT COUNT(DISTINCT sample_id) as samples_with_ploidy
FROM `{dataset}.sample_chromosome_ploidy`;

-- Find (table_name, sample_id) pairs present in GCS inventory but not yet loaded
-- (Run parse_and_group_files.py or verify_all_loaded.py for automated detection)
```

---

## Rollout Plan

### Phase 1: Prototype (Week 1)
- Implement tracking table schema
- Create `DiscoverParquetFiles` task
- Create basic `LoadParquetFilesToBQ` task (single table)
- Test with 10-100 files

### Phase 2: Scale Testing (Week 2)
- Add batching logic
- Implement `VerifyAllLoaded` task
- Test with 10,000 files
- Test preemption recovery

### Phase 3: Production Integration (Week 3)
- Integrate with `GvsImportGenomes.wdl`
- Full end-to-end testing with production data
- Document operational procedures

### Phase 4: Monitoring & Optimization (Week 4)
- Add monitoring and alerting
- Optimize batch sizes based on actual performance
- Create dashboards for load status

---

## Operational Procedures

### Normal Operation

1. Run `GvsImportGenomes.wdl` to generate Parquet files
2. Run `LoadParquetsToBigQuery.wdl` to load files into BigQuery
3. Verify `all_loaded = true` in workflow outputs
4. Check BigQuery tables for expected row counts

### Handling Failures

**If LoadParquetFilesToBQ fails:**
1. Check error message in task stderr
2. Verify BigQuery quotas and permissions
3. Re-run workflow (idempotent)

**If VerifyAllLoaded reports missing files:**
1. Check `missing_files_list` output
2. Investigate why files weren't loaded (errors in task outputs)
3. Re-run workflow to attempt loading missing files

**If file count doesn't match expectations:**
1. Query tracking table to see what was loaded
2. Compare against GCS listing
3. Investigate discrepancies in `GvsImportGenomes` logs

### Monitoring Queries

```sql
-- Overall load status: partitions loaded per table
SELECT
  table_name,
  COUNT(*) AS partitions_loaded,
  SUM(total_logical_bytes) AS total_bytes,
  MIN(last_modified_time) AS first_modified,
  MAX(last_modified_time) AS last_modified
FROM `{project}.{dataset}.INFORMATION_SCHEMA.PARTITIONS`
WHERE REGEXP_CONTAINS(table_name, "^vet_[0-9]+$|^ref_ranges_[0-9]+$")
  AND total_logical_bytes > 0
  AND partition_id NOT LIKE "__%"
GROUP BY table_name
ORDER BY table_name;

-- Ploidy samples loaded
SELECT COUNT(DISTINCT sample_id) AS samples_with_ploidy
FROM `{project}.{dataset}.sample_chromosome_ploidy`;
```

---

## Future Enhancements

1. **Incremental Loading**: Webhook or pub/sub trigger when new files are written
2. **Automated Retry**: Exponential backoff for failed load jobs
3. **Schema Validation**: Pre-validate Parquet schemas before loading
4. **Compression**: Optimize Parquet compression settings for BigQuery
5. **Partitioning**: Use BigQuery table partitioning for better query performance
6. **Streaming**: Investigate BigQuery streaming inserts for lower latency

---

## Summary

This design provides a robust, idempotent solution for loading Parquet files from GCS to BigQuery with the following key features:

✅ **Idempotent**: Safe to run multiple times  
✅ **Fault-Tolerant**: Survives preemptions and failures  
✅ **Scalable**: Handles hundreds of thousands of files  
✅ **Verifiable**: Clear completion status  
✅ **Maintainable**: Well-separated concerns with clear task boundaries  
✅ **Cost-Effective**: Uses preemptible VMs and free BigQuery loads  

The tracking table approach provides the best balance of reliability, queryability, and operational simplicity.
