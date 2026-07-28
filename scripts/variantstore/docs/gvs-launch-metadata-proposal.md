# Proposal: recording GVS launch metadata (VS-1961)

GVS does not record the inputs it was launched with. This proposal identifies which inputs matter,
proposes a schema for storing them alongside the data they describe, and sketches how the recorded
values would be used to fail fast when an input that invalidates joint calling changes.

## Scope

This proposal targets `GvsJointVariantCalling.wdl`, the top-level entry point for GVS Beta, as the
first workflow to record metadata.

That choice covers Beta but *not* AoU: AoU callsets never run `GvsJointVariantCalling.wdl`, they run
`GvsBulkIngestGenomes.wdl`, `GvsPopulateAltAllele.wdl`, `GvsCreateFilterSet.wdl`,
`GvsExtractAvroFilesForHail.wdl` and `GvsCreateVDS.wdl` individually (see `AOU_DELIVERABLES.md`). A
`GvsJointVariantCalling`-only implementation would therefore record nothing at all for the largest
and most operationally sensitive callsets we produce. The schema below is consequently
*workflow-agnostic*: one row per workflow launch, one row per input, with the workflow name as a
column. Wiring in the AoU entry points afterwards requires no schema change and no new tooling — just
one more task call per workflow.

## Why this is worth doing

Three concrete observations from the current code:

1. **Invariants are reverse-engineered from data, and only when they happen to leave a fingerprint.**
   `IsUsingCompressedReferences` (`GvsUtils.wdl:1084`) queries `INFORMATION_SCHEMA.COLUMNS` for
   `ref_ranges_001` to work out whether references were compressed at ingest. `IsVETS`
   (`GvsUtils.wdl:1014`) counts non-NULL `calibration_sensitivity` versus `vqslod` rows in
   `filter_set_info` to work out whether a filter set is VETS or VQSR. `GetExtractVetTableVersion`
   (`GvsUtils.wdl:1147`) does the same kind of archaeology for vet table layout. These work because
   those particular choices are visible in the physical schema. Choices that leave no fingerprint —
   `drop_state`, the interval list actually used, `is_wgs` — are simply unrecoverable after the fact.

2. **Defaults already disagree between entry points, silently.** `drop_state` defaults to `"FORTY"` in
   `GvsJointVariantCalling.wdl:28` but to `"NONE"` in `GvsBulkIngestGenomes.wdl:44` and
   `GvsImportGenomes.wdl`. `maximum_alternate_alleles` defaults to `1000` in
   `GvsJointVariantCalling.wdl:83` but `100` in `GvsExtractCallset.wdl:67`. Nothing detects a dataset
   whose samples were ingested under two different `drop_state` values.

3. **Provenance is currently a manual, human process.** `docs/aou/AOU_PROVENANCE_TRACKING.md`
   describes recording release/branch/workflow provenance by hand for each AoU callset. Everything
   recorded there by hand is available programmatically at launch time.

## Guiding principle: record everything, enforce a subset

Recording an input costs one row (~50 rows per workflow launch, i.e. nothing). The expensive and
interesting judgment call is not *what to store* but *what to enforce*. So rather than dropping
inputs we currently think are unimportant — and losing them permanently when we turn out to be wrong
— every input is recorded and each carries two attributes:

- `scope`: the blast radius over which the value must stay constant.
- `is_enforced`: whether a mismatch should stop the pipeline.

Scopes:

| Scope | Meaning | Consequence of a change |
| --- | --- | --- |
| `DATASET` | Must hold for the life of the BigQuery dataset | The store becomes internally inconsistent: samples are not comparable to each other |
| `CALLSET` | Must agree between filter set creation and extract for one callset | The callset is filtered or extracted inconsistently with how it was modeled |
| `RUN` | Provenance only | None; recorded for forensics and cost/attribution |

## Which inputs

### `DATASET` scope, enforced

These are the "invalidates joint calling" inputs the ticket asks about. Every one of them makes
samples ingested before the change non-comparable to samples ingested after it.

| Input | Consumed by | What breaks if it changes mid-dataset |
| --- | --- | --- |
| `project_id`, `dataset_name` | everything | Dataset identity; recorded as the run's subject |
| `reference_name` | ingest, filter set, extract | Coordinates and reference bases differ |
| `is_wgs` | ingest, extract | Selects WGS vs exome interval list, and exome-specific extract behavior |
| effective interval list (`interval_list_to_use`, i.e. `interval_list` or the `GetReference` default) | ingest, filter set, extract | Samples cover different regions; absence of data becomes ambiguous |
| `drop_state` | ingest, extract | Which GQ reference blocks were discarded; mixing yields heterogeneous reference data and silently wrong reference genotypes |
| `use_compressed_references` | ingest | Physical `ref_ranges_%` schema (`location` vs `packed_ref_data`); not mixable within a dataset |
| `load_vet_and_ref_ranges` | ingest | Which core tables were populated for which samples |
| `load_vcf_headers` | ingest | Whether per-sample header provenance exists |
| `use_parquet_ingest`, `parquet_output_gcs_dir` | ingest | Ingest path and where the intermediate data lives |

When `GvsBulkIngestGenomes.wdl` is wired in, `samples_are_controls`, `data_table_name` and
`bulk_ingest_fofn` join this list as sample-provenance inputs (which samples entered the store, and
whether they were flagged as controls).

### `CALLSET` scope, enforced

| Input | Consumed by | What breaks if it changes |
| --- | --- | --- |
| effective `filter_set_name` | filter set, extract | Extract applies a different model than the one created |
| `use_VQSR` (→ `use_VETS`) | filter set | VETS vs VQSR scores are not interchangeable |
| `add_additional_annotations_to_sites_only_vcf` | filter set | Model training features |
| `INDEL_VQSR_max_gaussians_override`, `SNP_VQSR_max_gaussians_override` | filter set | Model shape |
| `training_python_script`, `scoring_python_script` | filter set | The VETS model *code* itself — highest-leverage override we have |
| effective `extract_table_prefix` | prepare, extract | Extract reads different prepare tables than were written |
| `sample_names_to_extract` | prepare | Cohort definition |
| `extract_do_not_filter_override` | extract | Filters silently not applied |
| `target_interval_list` | extract | Exome/BGE target regions |
| `maximum_alternate_alleles` | extract | Truncation of multi-allelic sites (see the 1000 vs 100 default mismatch above) |
| `drop_state`, `is_wgs`, effective interval list | extract | Re-asserted at extract; must equal the dataset-scoped values |

### `RUN` scope, recorded only

`call_set_identifier`, `extract_output_gcs_dir`, `extract_output_file_base_name`,
`merge_output_vcfs`, `bgzip_output_vcfs`, `collect_variant_calling_metrics`, `sample_set_name`,
`sample_id_column_name`, `vcf_files_column_name`, `vcf_index_files_column_name`,
`billing_project_id`, `query_labels`, `tighter_gcp_quotas`, `gatk_override`, and the performance
knobs (`extract_scatter_count`, `extract_preemptible_override`, `extract_maxretries_override`,
`load_data_scatter_width`, `load_data_preemptible_override`, `load_data_maxretries_override`,
`split_intervals_disk_size_override`, `split_intervals_mem_override`, `INDEL_VQSR_mem_gb_override`,
`SNP_VQSR_mem_gb_override`).

`gatk_override` deserves a callout: it is recorded like anything else, but because a defined
`gatk_override` means the run did not use the production GATK binary, it is worth surfacing in a view
(`gatk_override IS NOT NULL` → not a production run).

### Non-input run context

Beyond the WDL inputs, each launch records:

| Field | Source |
| --- | --- |
| workflow name | `GetToolVersions`, see implementation notes — already in the path it scrapes |
| workflow ID (root), submission ID | `GetToolVersions.workflow_id` / `.submission_id` |
| workspace ID, bucket, Google project | `GetToolVersions.workspace_id` / `.workspace_bucket` / `.google_project` |
| workspace name, namespace | Rawls, via existing `scripts/get_workspace_name_for_import.py` |
| BigQuery dataset location, workspace bucket location | `bq show` / `gcloud storage buckets describe` |
| launch timestamp | `GetToolVersions.date_as_secs_since_unix_epoch` |
| `gvs_version`, `git_branch_or_tag`, `git_hash` | `GetToolVersions` |
| `gatk_docker`, `variants_docker`, `cloud_sdk_docker`, `basic_docker` | `GetToolVersions` (effective values, so caller overrides are captured) |

## Proposed schema

Two tables, created in the GVS dataset itself so that metadata travels with the data it describes
(and so a dataset copy carries its own provenance). Both are tiny and unpartitioned, and can be
created with the existing `GvsCreateTables.CreateTables` task with `partitioned = "false"` and
`superpartitioned = "false"`.

### `gvs_workflow_run` — one row per workflow launch

| Column | Type | Mode | Notes |
| --- | --- | --- | --- |
| `run_id` | STRING | REQUIRED | Cromwell root workflow ID (UUID); primary key by convention |
| `workflow_name` | STRING | REQUIRED | e.g. `GvsJointVariantCalling`, `GvsBulkIngestGenomes` |
| `call_set_identifier` | STRING | NULLABLE | Callset key where the workflow has one |
| `started` | TIMESTAMP | REQUIRED | Launch time |
| `submission_id` | STRING | REQUIRED | Terra submission ID |
| `terra_workspace_name` | STRING | NULLABLE | From Rawls |
| `terra_workspace_namespace` | STRING | NULLABLE | From Rawls |
| `terra_workspace_id` | STRING | NULLABLE | |
| `terra_workspace_bucket` | STRING | NULLABLE | `gs://fc-...` |
| `terra_google_project` | STRING | NULLABLE | |
| `terra_workspace_bucket_location` | STRING | NULLABLE | Region/multi-region of the workspace bucket |
| `bq_dataset_location` | STRING | NULLABLE | Region of the GVS BigQuery dataset |
| `gvs_version` | STRING | REQUIRED | `"unspecified"` for non-release runs |
| `git_branch_or_tag` | STRING | NULLABLE | |
| `git_hash` | STRING | REQUIRED | |
| `gatk_docker` | STRING | REQUIRED | Effective value |
| `variants_docker` | STRING | REQUIRED | Effective value |
| `cloud_sdk_docker` | STRING | REQUIRED | Effective value |
| `basic_docker` | STRING | REQUIRED | Effective value |
| `gatk_override` | STRING | NULLABLE | Non-NULL means a non-production run |

GVS-style schema JSON:

```
[{"name":"run_id","type":"STRING","mode":"REQUIRED"},{"name":"workflow_name","type":"STRING","mode":"REQUIRED"},{"name":"call_set_identifier","type":"STRING","mode":"NULLABLE"},{"name":"started","type":"TIMESTAMP","mode":"REQUIRED"},{"name":"submission_id","type":"STRING","mode":"REQUIRED"},{"name":"terra_workspace_name","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_namespace","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_id","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_bucket","type":"STRING","mode":"NULLABLE"},{"name":"terra_google_project","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_bucket_location","type":"STRING","mode":"NULLABLE"},{"name":"bq_dataset_location","type":"STRING","mode":"NULLABLE"},{"name":"gvs_version","type":"STRING","mode":"REQUIRED"},{"name":"git_branch_or_tag","type":"STRING","mode":"NULLABLE"},{"name":"git_hash","type":"STRING","mode":"REQUIRED"},{"name":"gatk_docker","type":"STRING","mode":"REQUIRED"},{"name":"variants_docker","type":"STRING","mode":"REQUIRED"},{"name":"cloud_sdk_docker","type":"STRING","mode":"REQUIRED"},{"name":"basic_docker","type":"STRING","mode":"REQUIRED"},{"name":"gatk_override","type":"STRING","mode":"NULLABLE"}]
```

### `gvs_workflow_run_input` — one row per input per launch

| Column | Type | Mode | Notes |
| --- | --- | --- | --- |
| `run_id` | STRING | REQUIRED | Joins to `gvs_workflow_run` |
| `input_name` | STRING | REQUIRED | WDL input name, e.g. `drop_state`; effective values use the input's own name, not `effective_*` |
| `effective_value` | STRING | NULLABLE | Canonical string rendering; NULL means an optional input with no value |
| `value_type` | STRING | REQUIRED | `STRING` \| `BOOLEAN` \| `INT` \| `FLOAT` \| `FILE` \| `ARRAY` |
| `scope` | STRING | REQUIRED | `DATASET` \| `CALLSET` \| `RUN` |
| `is_enforced` | BOOLEAN | REQUIRED | Whether a mismatch should fail a subsequent run |
| `was_specified` | BOOLEAN | NULLABLE | `false` when the WDL default applied — distinguishes "chose the default" from "default changed under us" |
| `file_generation` | STRING | NULLABLE | GCS object generation, for `FILE` values |
| `file_md5` | STRING | NULLABLE | GCS `md5Hash`, for `FILE` values |
| `file_size_bytes` | INTEGER | NULLABLE | For `FILE` values |

GVS-style schema JSON:

```
[{"name":"run_id","type":"STRING","mode":"REQUIRED"},{"name":"input_name","type":"STRING","mode":"REQUIRED"},{"name":"effective_value","type":"STRING","mode":"NULLABLE"},{"name":"value_type","type":"STRING","mode":"REQUIRED"},{"name":"scope","type":"STRING","mode":"REQUIRED"},{"name":"is_enforced","type":"BOOLEAN","mode":"REQUIRED"},{"name":"was_specified","type":"BOOLEAN","mode":"NULLABLE"},{"name":"file_generation","type":"STRING","mode":"NULLABLE"},{"name":"file_md5","type":"STRING","mode":"NULLABLE"},{"name":"file_size_bytes","type":"INTEGER","mode":"NULLABLE"}]
```

The GCS fingerprint columns are the reason `File` inputs need more than a path. An interval list, a
`sample_names_to_extract` file or a `training_python_script` can be overwritten in place at the same
URI; the object generation is the only reliable answer to "is this the same interval list we ingested
against?".

### Optional: `gvs_workflow_run_event`

A row is written to `gvs_workflow_run` at launch, before any real work starts, so that failed and
aborted runs are visible too. Terminal status, if wanted, is best appended rather than updated —
mirroring the existing `sample_load_status` pattern (`sample_id`, `status`, `event_timestamp`) which
avoids DML against rows that may still be in the streaming buffer:

```
[{"name":"run_id","type":"STRING","mode":"REQUIRED"},{"name":"event","type":"STRING","mode":"REQUIRED"},{"name":"event_timestamp","type":"TIMESTAMP","mode":"REQUIRED"}]
```

### Why one row per input rather than one column per input

A wide table (one column per WDL input) is more pleasant to query but a poor fit here:

- `GvsJointVariantCalling` alone has ~45 inputs, and they change from release to release. A wide
  table needs a schema migration each time an input is added or renamed.
- Made workflow-agnostic, a wide table becomes the union of every entry point's inputs — hundreds of
  mostly-NULL columns.
- The consistency check (below) can be written once, generically, without enumerating column names.

The cost is that values are stringly typed. Canonical rendering rules keep that manageable: booleans
as `true`/`false`, arrays as JSON array text, files as `gs://` URIs, NULL reserved for "no value".
Views can present the typed, pivoted shape where it's convenient.

## How the values get recorded

A new `WriteWorkflowRunMetadata` task in `GvsUtils.wdl`, called from `GvsJointVariantCalling.wdl`
immediately after `GetToolVersions`/`GetReference` and before `GvsBulkIngestGenomes`, with its `done`
output gating ingest (the established GVS ordering idiom).

Implementation details that matter:

- **Declare the task's inputs as `String`, never `File`.** `interval_list_to_use`,
  `target_interval_list`, `sample_names_to_extract` and the VETS scripts are `File` at the workflow
  level; passing them to a `File` task input would localize them and record a container-local path.
  Coercing to `String` at the call site records the `gs://` URI and costs nothing.
- **Compute effective values in the workflow, not the task.** `GvsJointVariantCalling.wdl` already
  computes `effective_*` for dockers, git hash, workspace and the extract names; the same pattern
  applies to everything else, and `defined(x)` supplies `was_specified`.
- **Load with `bq load`, not DML.** Emit newline-delimited JSON and
  `bq --apilog=false load --source_format=NEWLINE_DELIMITED_JSON`, following the precedent at
  `GvsUtils.wdl:859`. No `INSERT` means no quoting hazards and no streaming-buffer interactions.
- **File fingerprints** come from `gcloud storage objects describe <uri> --format="value(generation,md5_hash,size)"`,
  tolerating failure (a missing or inaccessible object should not fail the run; record the URI with
  NULL fingerprints).
- **`GetToolVersions` needs only one addition, `workflow_name`.** The delocalization path it already
  scrapes is `gs://fc-<workspace id>/submissions/<submission id>/<workflow name>/<workflow id>/...`,
  and the existing regex at `GvsUtils.wdl:96` already captures the workflow name as group 4 — the new
  output is the same `sed` with `\4` instead of `\5`.
- **Workspace name/namespace and region belong in the new task, not `GetToolVersions`.**
  `GetToolVersions` runs in the alpine `cloud-sdk` image, which has neither `requests` nor
  `terra_notebook_utils`. The new task can run in `variants_docker`, where
  `/app/get_workspace_name_for_import.py` already returns workspace name and namespace given a
  workspace ID. Region lookups (`bq show --format=prettyjson`, `gcloud storage buckets describe`) can
  live in the same task.

## How the recorded metadata gets used

The payoff is a `CheckDatasetInvariants` task that runs before ingest: read the enforced
`DATASET`-scoped values from the most recent prior run against this dataset, compare to this run's,
and fail with a readable diff. An `allow_dataset_invariant_change` escape hatch input covers
deliberate changes — and being an input, it is itself recorded, so a deliberate change leaves a
trail.

Drift detection across a dataset's whole history is a single generic query:

```sql
SELECT
  i.input_name,
  ARRAY_AGG(DISTINCT IFNULL(i.effective_value, '<none>')) AS distinct_values,
  COUNT(DISTINCT i.run_id) AS runs
FROM `PROJECT.DATASET.gvs_workflow_run_input` i
WHERE i.scope = 'DATASET' AND i.is_enforced
GROUP BY i.input_name
HAVING COUNT(DISTINCT IFNULL(i.effective_value, '<none>')) > 1
```

Longer term this replaces inference with lookup: `IsUsingCompressedReferences` and `IsVETS` become
metadata reads, with the existing inference retained as a fallback for datasets that predate the
table. Two notes on that transition:

- **Backfill.** Existing datasets have no rows, so "no prior run recorded" must be treated as a pass,
  not a failure. Where it's worth it, a single synthetic row per legacy dataset can be seeded from
  what the existing inference tasks can still determine.
- **Warn before failing.** See open questions.

## Open questions for review

1. **Which entry points in the first implementation?** Recommendation: build the workflow-agnostic
   table now and wire up `GvsJointVariantCalling` and `GvsBulkIngestGenomes` first, since the latter
   is what AoU actually runs and is where the dataset-scoped ingest invariants are set.
2. **Warn or fail in v1?** Recommendation: record and warn in the first release, and flip enforced
   mismatches to hard failures once datasets have accumulated metadata — enforcement against an
   almost-empty table mostly produces false negatives, and any bug in the check produces a hard
   outage.
3. **Where do the tables live?** Per-GVS-dataset (recommended: provenance travels with the data,
   no cross-project permissions) versus one central table across all datasets (easier fleet-wide
   reporting). These are not mutually exclusive — a central view can union per-dataset tables.
4. **Is `call_set_identifier` a good enough callset key** for `CALLSET`-scoped enforcement, or do we
   need an explicit callset ID? It is user-supplied and gets mangled (spaces and underscores to
   hyphens) for naming purposes.
5. **Do we want submitter identity?** Available from the Rawls submission API and useful for
   forensics, but worth a deliberate decision for AoU workspaces before recording it.
