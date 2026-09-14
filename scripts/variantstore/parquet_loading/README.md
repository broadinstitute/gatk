# Parquet-Based BigQuery Loading

## Overview

GVS previously ingested data into BigQuery row-by-row via the BigQuery Storage Write API, inside the `LoadData` task. This approach was replaced with a two-phase design:

1. **`LoadData`** writes variant and reference data to local Parquet files and copies them to a temporary GCS path.
2. A downstream set of WDL tasks (described below) loads those Parquet files into BigQuery in bulk, using BigQuery's free GCS-to-BQ load jobs rather than the metered Write API.

The loading system is idempotent and survives VM preemptions without duplicating data.

## WDL Tasks (`GvsImportGenomes.wdl`)

Four tasks coordinate loading, verification, and cleanup after the `LoadData` scatter completes:

| Task                   | Purpose                                                                                                                                                                                                                                                                                                 |
|------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `DiscoverParquetFiles` | Lists all Parquet files in GCS, filters out already-loaded `(table_name, sample_id)` pairs by querying BigQuery directly, and emits one FOFN per target table.                                                                                                                                          |
| `LoadParquetFilesToBQ` | Scattered across all tables; loads each FOFN into its BigQuery table in batches of up to 10,000 files per load job.                                                                                                                                                                                     |
| `VerifyParquetLoading` | Re-queries BigQuery after all loads complete, runs exact structural checks (`completeness`, `cardinality`, `cross_family`) to gate `all_loaded`, and evaluates cohort duplication/truncation screens to gate `safe_to_delete_parquet`. Fails the workflow if any expected data is missing or corrupted. |
| `DeleteParquetFiles`   | Deletes temporary GCS Parquet files only if `safe_to_delete_parquet` is true.                                                                                                                                                                                                                           |

## Python Scripts (`scripts/variantstore/scripts/`)

All scripts are packaged into the variants Docker image at `/app/`.

| Script                        | Description                                                                                                                                                                                                                                                                                                                                  |
|-------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `parse_and_group_files.py`    | Discovers Parquet files under a GCS prefix, maps each file to its target BigQuery table and `sample_id`, queries `INFORMATION_SCHEMA.PARTITIONS` (and `sample_chromosome_ploidy` directly) to identify already-loaded pairs, and writes per-table FOFNs for the remaining files.                                                             |
| `load_parquet_to_bq.py`       | Loads one FOFN into BigQuery. Generates a deterministic SHA-1 job ID per batch so that a retry VM after preemption re-uses the in-flight BigQuery job rather than re-submitting it. Retries quota and transient errors with exponential backoff (30 s → 60 s → 120 s, up to 3 retries).                                                      |
| `verify_all_loaded.py`        | CLI coordinator that confirms every `(table_name, sample_id)` pair is loaded, executes structural verification via `verify_structural_checks.py`, logs detailed diagnostic summaries, writes `verification_results.json` and `missing_files.txt`, and exits non-zero if `all_loaded` is false.                                               |
| `verify_structural_checks.py` | Independent structural verification engine (VS-1989). Evaluates BigQuery partition and table counts for family completeness (no empty partitions), ploidy cardinality (when an exact override is specified), cross-family consistency across co-produced tables, and cohort-relative duplication and truncation screens on `vet` row counts. |

## Verification and Deletion Gates (VS-1989)

Verification enforces a two-tier safety contract separating workflow success from data deletion authorization:

### 1. `all_loaded` (Workflow Gate)
Determines whether BigQuery ingestion succeeded. Governed exclusively by **exact** checks:
- **File Presence**: Every `(table_name, sample_id)` pair discovered in GCS must be confirmed present in BigQuery.
- **Family Completeness**: Every expected sample partition in superpartitioned tables (`vet_%`, `ref_ranges_%`) and every expected sample row in regular tables (`sample_chromosome_ploidy`, `vcf_header_lines_scratch`) must contain $> 0$ rows.
- **Ploidy Cardinality Consistency**: Since `SamplePloidyCreator` emits at most one row per chromosome per sample, `COUNT(*) == COUNT(DISTINCT chromosome)` is verified per sample to detect duplicated row ingestion exactly while allowing legitimate contig variations (such as 23 vs 24 or 25 with chrM). When an explicit exact count is configured via `--expected-ploidy-rows-per-sample` (e.g. `24` for WGS), every expected sample must also match that exact count.
- **Cross-Family Consistency**: A sample present in any co-produced family (`vet`, `ref_ranges`, `sample_chromosome_ploidy`, and conditionally `vcf_header_lines_scratch` when loaded alongside data) must be present in all configured co-produced families.

If any exact check fails, `all_loaded` is `false`, and `VerifyParquetLoading` exits non-zero, immediately aborting the workflow.

### 2. `safe_to_delete_parquet` (Deletion Gate)
A conservative gate governing whether `DeleteParquetFiles` is permitted to delete the source Parquet files from GCS:
- Requires `all_loaded == true`.
- Evaluates the cohort duplication screen (default ceiling: $1.6\times$ median / lower baseline for $N=2$) and truncation screen (default floor: $1/1.6 = 0.625\times$ median / upper baseline for $N=2$) on `vet` row counts.
- Singletons ($N=1$) are flagged conservatively by default since no cohort consensus baseline exists.
- Flagged samples retain the source Parquet files in GCS for investigation (`safe_to_delete_parquet = false`), but do not fail `all_loaded`.
- If flagged row counts are verified as legitimate (e.g. high-coverage samples or expected singletons), deletion can be explicitly authorized via `allow_flagged_vet_loads = true` (`--allow-flagged-vet-loads`).

### Diagnostic Persistence
On verification failure (`rc != 0`), `VerifyParquetLoading` copies `verification_results.json` to `verification_diagnostics_gcs_dir` (if specified) before exiting. This copy is gated on failure so successful retries cannot overwrite the durable failure diagnostics.

## Key Design Points

**Idempotency** is based on what BigQuery actually contains, not a secondary tracking table. `parse_and_group_files.py` queries `INFORMATION_SCHEMA.PARTITIONS` for non-empty integer partitions in the `vet_%` and `ref_ranges_%` tables, and queries `sample_chromosome_ploidy` directly for distinct `sample_id` values. Files belonging to already-loaded pairs are silently skipped on every run.

**Preemption recovery** relies on deterministic BigQuery job IDs (a SHA-1 hash of the project, dataset, table name, and sorted batch file list). When a retry VM submits the same batch, BigQuery returns a `Conflict` exception; `load_parquet_to_bq.py` catches this and fetches the pre-existing job to wait on, so no data is written twice and no local state file is required.

## Further Reading

- [BigQuery: Loading Parquet data from Cloud Storage](https://cloud.google.com/bigquery/docs/loading-data-cloud-storage-parquet)
- [BigQuery quotas and limits — load jobs](https://cloud.google.com/bigquery/quotas#load_jobs)
