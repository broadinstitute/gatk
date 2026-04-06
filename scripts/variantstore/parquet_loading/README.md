# Parquet-Based BigQuery Loading

## Overview

GVS previously ingested data into BigQuery row-by-row via the BigQuery Storage Write API, inside the `LoadData` task.  This approach was replaced with a two-phase design:

1. **`LoadData`** writes variant and reference data to local Parquet files and copies them to a temporary GCS path.
2. A downstream set of WDL tasks (described below) loads those Parquet files into BigQuery in bulk, using BigQuery's free GCS-to-BQ load jobs rather than the metered Write API.

The loading system is idempotent and survives VM preemptions without duplicating data.

## WDL Tasks (`GvsImportGenomes.wdl`)

Three tasks run in sequence after the `LoadData` scatter completes:

| Task | Purpose |
|---|---|
| `DiscoverParquetFiles` | Lists all Parquet files in GCS, filters out already-loaded `(table_name, sample_id)` pairs by querying BigQuery directly, and emits one FOFN per target table. |
| `LoadParquetFilesToBQ` | Scattered across all tables; loads each FOFN into its BigQuery table in batches of up to 10,000 files per load job. |
| `VerifyParquetLoading` | Re-queries BigQuery after all loads complete and fails the workflow if any `(table_name, sample_id)` pair is still missing. |

## Python Scripts (`scripts/variantstore/scripts/`)

All scripts are packaged into the variants Docker image at `/app/`.

| Script | Description |
|---|---|
| `parse_and_group_files.py` | Discovers Parquet files under a GCS prefix, maps each file to its target BigQuery table and `sample_id`, queries `INFORMATION_SCHEMA.PARTITIONS` (and `sample_chromosome_ploidy` directly) to identify already-loaded pairs, and writes per-table FOFNs for the remaining files. |
| `load_parquet_to_bq.py` | Loads one FOFN into BigQuery. Generates a deterministic SHA-1 job ID per batch so that a retry VM after preemption re-uses the in-flight BigQuery job rather than re-submitting it. Retries quota and transient errors with exponential backoff (30 s → 60 s → 120 s, up to 3 retries). |
| `verify_all_loaded.py` | Confirms that every `(table_name, sample_id)` pair derived from the original GCS file list is present in BigQuery; reports any missing pairs and exits non-zero if any are found. |

## Key Design Points

**Idempotency** is based on what BigQuery actually contains, not a secondary tracking table. `parse_and_group_files.py` queries `INFORMATION_SCHEMA.PARTITIONS` for non-empty integer partitions in the `vet_%` and `ref_ranges_%` tables, and queries `sample_chromosome_ploidy` directly for distinct `sample_id` values. Files belonging to already-loaded pairs are silently skipped on every run.

**Preemption recovery** relies on deterministic BigQuery job IDs (a SHA-1 hash of the project, dataset, table name, and sorted batch file list). When a retry VM submits the same batch, BigQuery returns a `Conflict` exception; `load_parquet_to_bq.py` catches this and fetches the pre-existing job to wait on, so no data is written twice and no local state file is required.

## Further Reading

- [BigQuery: Loading Parquet data from Cloud Storage](https://cloud.google.com/bigquery/docs/loading-data-cloud-storage-parquet)
- [BigQuery quotas and limits — load jobs](https://cloud.google.com/bigquery/quotas#load_jobs)
