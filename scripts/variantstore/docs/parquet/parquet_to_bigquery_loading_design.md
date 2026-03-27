# Design Document: Robust Parquet to BigQuery Loading

## Overview

After `CreateVariantIngestFiles` writes Parquet files to GCS, a set of WDL tasks
loads them into BigQuery tables robustly and idempotently. The system handles
preemptions, partial completions, and workflow restarts without duplicating work.

## GCS Directory Structure

Files are organised under `${OUTPUT_GCS_DIR}/{table_prefix}/{table_number}/` and
map directly to BigQuery tables of the form `{table_prefix}_{table_number}` (e.g.
`vet/001/…` → `vet_001`, `ref_ranges/042/…` → `ref_ranges_042`).

## Design Principles

- **Idempotent**: Safe to run multiple times; already-loaded files are skipped.
- **Fault-Tolerant**: Survives preemptions; no work is lost or duplicated.
- **Scalable**: Handles hundreds of thousands of files across many tables.
- **Verifiable**: A final verification step confirms every file was loaded.

## Architecture

Three WDL tasks run in sequence:

### 1. DiscoverParquetFiles

Lists all Parquet files in GCS and groups them by target BigQuery table. Before
producing its output FOFNs it queries BigQuery to identify which
`(table_name, sample_id)` pairs are already loaded and excludes those files.

Idempotency is determined by querying BigQuery directly — no separate tracking
table is created or maintained:

- **Superpartitioned tables** (`vet_%`, `ref_ranges_%`): `INFORMATION_SCHEMA.PARTITIONS`
  is checked for non-empty partitions whose `partition_id` is a valid integer
  `sample_id`.
- **Regular tables** (`sample_chromosome_ploidy`): the table itself is queried
  for distinct `sample_id` values.

This reflects the actual committed state of BigQuery data rather than a secondary
record of intent, and avoids the high-frequency MERGE DML writes that caused
quota pressure at scale. The task is marked `volatile: true` so it always
re-scans GCS and BigQuery.

**Outputs**: one FOFN per table containing only the files still to be loaded,
plus summary statistics.

### 2. LoadParquetFilesToBQ

Scattered across all tables, each instance loads the files listed in its FOFN
into its target BigQuery table, processing them in batches of up to 10,000 files
per BigQuery load job.

**Fault tolerance** is provided by deterministic job IDs: each batch's job ID is
a SHA-1 hash of the table name and sorted batch contents. If a VM is preempted
and retried, the same job IDs are regenerated. BigQuery returns a Conflict
exception for any already-submitted job, and the script fetches and waits on the
existing job rather than resubmitting. No local state file is needed or written.

Quota and transient errors (e.g. `TooManyRequests`, `ServiceUnavailable`) are
retried with exponential backoff (30 s → 60 s → 120 s). A failing batch is
logged and skipped so that the remaining batches continue; the task exits with a
partial-success status if any batch failed.

**Preemptible**: set to 3–5 attempts to save cost; fully resumable.

### 3. VerifyParquetLoading

After all load tasks complete, this task re-runs the same BigQuery partition
query and reports any GCS files whose `(table_name, sample_id)` pair is still
absent. The task fails if any files are missing, triggering a workflow error with
a list of unloaded paths.

## Integration with GvsImportGenomes.wdl

The three tasks above are called from `GvsImportGenomes.wdl` after the
`LoadData` scatter completes. They can also be re-run independently if loading
needs to be retried without repeating data creation.

## Edge Cases

| Scenario | Handling |
|---|---|
| VM preempted mid-load | Deterministic job IDs allow the retry VM to resume from BigQuery rather than resubmit |
| Workflow restarted | `DiscoverParquetFiles` skips already-loaded `(table_name, sample_id)` pairs |
| Batch quota error | Exponential backoff retries up to 3 times; batch is skipped if exhausted |
| Missing table | Task fails fast rather than attempting to create the table |
| Schema mismatch | BigQuery raises an error immediately; fix the upstream Parquet producer |

