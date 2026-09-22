# Parquet-Based BigQuery Loading

## Overview

GVS previously ingested data into BigQuery row-by-row via the BigQuery Storage Write API, inside the `LoadData` task. This approach was replaced with a two-phase design:

1. **`LoadData`** writes variant and reference data to local Parquet files and copies them to a temporary GCS path.
2. A downstream set of WDL tasks (described below) loads those Parquet files into BigQuery in bulk, using BigQuery's free GCS-to-BQ load jobs rather than the metered Write API.

The loading system is idempotent and survives VM preemptions without duplicating data.

## WDL Tasks (`GvsImportGenomes.wdl`)

Five tasks coordinate loading, verification, and cleanup after the `LoadData` scatter completes:

| Task                            | Purpose                                                                                                                                                                                                                                                                                                          |
|---------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `DiscoverParquetFiles`          | Lists all Parquet files in GCS, filters out already-loaded `(table_name, sample_id)` pairs by querying BigQuery directly, and emits one FOFN per target table.                                                                                                                                                   |
| `LoadParquetFilesToBQ`          | Scattered across all tables; loads each FOFN into its BigQuery table in batches of up to 10,000 files per load job.                                                                                                                                                                                              |
| `VerifyParquetLoading`          | Re-queries BigQuery after all loads complete, runs exact structural checks (`completeness`, `cardinality`, `cross_family`) to gate `all_loaded`, and evaluates cohort duplication/truncation screens to produce the per-sample quarantine list. Fails the workflow if any expected data is missing or corrupted. |
| `QuarantineFlaggedParquetFiles` | Moves the Parquet of screen-flagged samples into a `quarantine/` subdirectory, renamed so it survives the bulk delete. Runs unconditionally and is a no-op when nothing was flagged.                                                                                                                             |
| `DeleteParquetFiles`            | Deletes temporary GCS Parquet files only if `safe_to_delete_parquet` is true, and only after the quarantine step has succeeded.                                                                                                                                                                                  |

## Python Scripts (`scripts/variantstore/scripts/`)

All scripts are packaged into the variants Docker image at `/app/`.

| Script                        | Description                                                                                                                                                                                                                                                                                                                                                               |
|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `parse_and_group_files.py`    | Discovers Parquet files under a GCS prefix, maps each file to its target BigQuery table and `sample_id`, queries `INFORMATION_SCHEMA.PARTITIONS` (and `sample_chromosome_ploidy` directly) to identify already-loaded pairs, and writes per-table FOFNs for the remaining files.                                                                                          |
| `load_parquet_to_bq.py`       | Loads one FOFN into BigQuery. Generates a deterministic SHA-1 job ID per batch so that a retry VM after preemption re-uses the in-flight BigQuery job rather than re-submitting it. Retries quota and transient errors with exponential backoff (30 s → 60 s → 120 s, up to 3 retries).                                                                                   |
| `verify_all_loaded.py`        | CLI coordinator that confirms every `(table_name, sample_id)` pair is loaded, executes structural verification via `verify_structural_checks.py`, logs detailed diagnostic summaries, writes `verification_results.json`, `missing_files.txt` and `quarantine_files.txt`, and exits non-zero if `all_loaded` is false.                                                    |
| `verify_structural_checks.py` | Independent structural verification engine (VS-1989). Evaluates BigQuery partition and table counts for family completeness (no empty partitions), ploidy cardinality (against an exact override when given, otherwise the callset mode), cross-family consistency across co-produced tables, and cohort-relative duplication and truncation screens on `vet` row counts. |

## Verification and Deletion Gates (VS-1989)

Verification enforces a two-tier safety contract separating workflow success from data deletion authorization:

### 1. `all_loaded` (Workflow Gate)
Determines whether BigQuery ingestion succeeded. Governed exclusively by **exact** checks:
- **File Presence**: Every `(table_name, sample_id)` pair discovered in GCS must be confirmed present in BigQuery.
- **Family Completeness**: Every expected sample partition in superpartitioned tables (`vet_%`, `ref_ranges_%`) and every expected sample row in regular tables (`sample_chromosome_ploidy`, `vcf_header_lines_scratch`) must contain $> 0$ rows.
- **Ploidy Cardinality Consistency**: Detects both duplicated loads and below-mode partial loads. Since `SamplePloidyCreator` emits at most one row per chromosome per sample, `COUNT(*) == COUNT(DISTINCT chromosome)` is verified per sample to detect chromosome collisions. In mode-inferred mode (the default), samples are screened against a valid contig floor ($\ge \text{mode} - 2$, catching truncated loads missing autosomes while permitting legitimate karyotype variations such as 23 vs 24 contigs) and a modal duplication ceiling ($< 1.5\times$ baseline). When an explicit exact count is configured via `--expected-ploidy-rows-per-sample` (e.g. `24` for WGS), every expected sample must match that exact count.
- **Cross-Family Consistency**: A sample present in any co-produced family (`vet`, `ref_ranges`, `sample_chromosome_ploidy`, and conditionally `vcf_header_lines_scratch` when loaded alongside data) must be present in all configured co-produced families.

If any exact check fails, `all_loaded` is `false`, and `VerifyParquetLoading` exits non-zero, immediately aborting the workflow.

### 2. The heuristic screens and the per-sample quarantine
The duplication and truncation screens are heuristics, so a flag holds back *that sample's* Parquet rather than blocking deletion for the callset:
- The duplication screen flags `vet` row counts above a ceiling (`parquet_vet_duplication_threshold`, default $1.6\times$ median, or the lower of the two values for $N=2$); the truncation screen flags those below a floor (`parquet_vet_truncation_threshold`, default $1/1.6 = 0.625\times$).
- **The two thresholds are separate inputs that happen to share a default.** Only the high side has been calibrated: the Foxtrot pass measured $1.6\times$ flagging 2 `vet` samples out of 540,545, and nothing has been measured below the median. A variant-count distribution has no reason to be symmetric in ratio space — its upper tail is bounded by biology, while its lower tail absorbs low-coverage samples, smaller callable fractions, and more aggressive GQ dropping — so mirroring $1.6$ into the floor is a provisional default, not a derived one. Keeping the knobs apart is what lets the low side be retuned, or switched off with `parquet_vet_truncation_threshold = 0` (`--vet-truncation-threshold 0`), without disturbing the calibrated high side. A disabled screen is logged and reported as disabled rather than as clean. This is also why gate 4 below aborts only on a duplication flag.
- Singletons ($N=1$) are flagged conservatively, since no cohort consensus baseline exists.
- A flag never affects `all_loaded`, and never blocks the delete for the callset; whether it aborts the run for attention is a separate decision, made by gate 4 below. `QuarantineFlaggedParquetFiles` moves every Parquet file belonging to a flagged sample — across **all** families, not just the flagging one, since a sample is only re-ingestable from a complete set of its files — into a `quarantine/` subdirectory of the Parquet output directory, appending `.quarantined` to each name. The rest of the callset's Parquet is then deleted as usual.
- If flagged row counts are verified as legitimate (e.g. high-coverage samples or expected singletons), the screens can be waived with `allow_flagged_vet_loads = true` (`--allow-flagged-vet-loads`), in which case nothing is quarantined and flagged Parquet is deleted with the rest.

Two independent properties keep quarantined files alive, and both matter:
- **The rename.** `DeleteParquetFiles`' default strategy globs `*.parquet` across the whole output directory, so a quarantined object survives it only because its name no longer ends in `.parquet`. The task asserts that the configured suffix does not itself end in `.parquet`.
- **The location.** `DeleteParquetFiles`' alternate strategy deletes the four table directories, and the bucket lifecycle rule `ConfigureParquetLifecycle` installs matches those same four prefixes with a 14-day `Delete` action. `quarantine/` is in neither list, so quarantined files are exempt from both. **Nothing removes them**: clean the directory up by hand once the samples have been reviewed. A `README.txt` written alongside them explains how to restore one.

Quarantining per sample rather than withholding the run's delete is what keeps a heuristic out of a callset-wide decision. On Foxtrot the calibration pass flagged 2 `vet` samples out of 540,545; under a run-wide gate that would have retained every sample's Parquet — and only until the 14-day lifecycle rule deleted it anyway, silently, in a green run.

### 3. `safe_to_delete_parquet` (Deletion Gate)
The gate governing whether `DeleteParquetFiles` runs at all:
- Requires `all_loaded == true`.
- Requires every flagged sample to have at least one resolvable Parquet path, so the quarantine step has something to move. A flagged sample with no such path blocks deletion, because running it would destroy exactly the files the screens asked to keep.
- Screen flags alone do **not** block it; see the quarantine above.
- `DeleteParquetFiles` additionally consumes `QuarantineFlaggedParquetFiles.done`, so a failed quarantine takes the delete down with it. That ordering is structural rather than a boolean, because a failed quarantine followed by a bulk delete is the one unsafe sequence.

### 4. `parquet_fail_on_quarantine` (Notification Gate)
Quarantining makes a flagged run *succeed*, which leaves the opposite problem: a run that needs a human look is green, and the only trace is a workflow output nobody reads. `parquet_fail_on_quarantine` (default `true`, forwarded from `GvsJointVariantCalling` and `GvsBulkIngestGenomes`) closes that gap by aborting such a run via `Utils.TerminateWorkflow`.
- The abort is **purely a notification**. It runs after `QuarantineFlaggedParquetFiles`, so the samples are loaded and their Parquet is already safe: nothing is destroyed and nothing is rolled back. Re-running with `parquet_fail_on_quarantine = false` (or `parquet_allow_flagged_vet_loads = true`, once the flag has been reviewed and dismissed) finishes from a fully loaded dataset.
- It is conditioned on `QuarantineFlaggedParquetFiles.quarantined_files`, not on the equivalent `VerifyParquetLoading` count. The data dependency is what places the abort after the move; terminating off the verification output would let Cromwell kill the workflow while the flagged files were still under the prefixes the 14-day lifecycle rule matches, neither quarantined nor deleted.
- `DeleteParquetFiles` depends on the same task, so it and the abort are siblings and Cromwell may start the delete first. That race is benign: the delete skips quarantined files by both suffix and location, and deleting the unflagged samples' Parquet is wanted either way.
- **Only duplication flags abort.** A truncation-only flag still quarantines and still reports, but the truncation threshold has only been measured on its high side, and an uncalibrated heuristic should not be able to fail a completed 500k-sample ingest. `VerifyParquetLoading` reports the duplication subset separately (`quarantine_duplication_sample_count`) for exactly this reason.

### Observability on a successful run
When the abort is off, or when only the truncation screen fired, the run is green and the flags and quarantine counts are published as workflow outputs (`parquet_vet_duplication_flagged`, `parquet_vet_truncation_flagged`, `parquet_quarantined_samples`, `parquet_quarantined_duplication_samples`, `parquet_quarantined_files`, `parquet_quarantine_files_list`), from `GvsImportGenomes` up through `GvsBulkIngestGenomes` and `GvsJointVariantCalling`. A non-zero `parquet_quarantined_samples` on a green run is the signal to go and look.

`GvsQuickstartVcfIntegration` asserts that count is zero. The quickstart callset is ordinary, so a flag there is a screen regression rather than bad input data — and since a truncation-only flag does not abort, the assertion is what keeps that case from passing CI silently.

### Diagnostic Persistence
On verification failure (`rc != 0`), `VerifyParquetLoading` copies `verification_results.json` to `verification_diagnostics_gcs_dir` (if specified) before exiting. This copy is gated on failure so successful retries cannot overwrite the durable failure diagnostics.

A green run that quarantined anything gets the same treatment, under distinct names (`verification_results.quarantine.json` and `quarantine_files.txt`). Without it the only record of why those files were set aside would be the Cromwell execution directory, which is exactly what gets cleaned up. The separate names preserve the invariant above: this path does re-run on a successful retry, so it must not be able to overwrite a preserved failure diagnostic.

## Key Design Points

**Idempotency** is based on what BigQuery actually contains, not a secondary tracking table. `parse_and_group_files.py` queries `INFORMATION_SCHEMA.PARTITIONS` for non-empty integer partitions in the `vet_%` and `ref_ranges_%` tables, and queries `sample_chromosome_ploidy` directly for distinct `sample_id` values. Files belonging to already-loaded pairs are silently skipped on every run.

**Preemption recovery** relies on deterministic BigQuery job IDs (a SHA-1 hash of the project, dataset, table name, and sorted batch file list). When a retry VM submits the same batch, BigQuery returns a `Conflict` exception; `load_parquet_to_bq.py` catches this and fetches the pre-existing job to wait on, so no data is written twice and no local state file is required.

## Further Reading

- [BigQuery: Loading Parquet data from Cloud Storage](https://cloud.google.com/bigquery/docs/loading-data-cloud-storage-parquet)
- [BigQuery quotas and limits — load jobs](https://cloud.google.com/bigquery/quotas#load_jobs)
