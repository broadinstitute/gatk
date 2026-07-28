# Proposal: recording GVS launch metadata (VS-1961)

GVS does not record the inputs it was launched with. This proposal identifies which inputs matter,
proposes a schema for storing them, and separates the two distinct needs that get conflated when this
is described as a single piece of work.

## Two motivations, two deliverables

The ticket describes the goal as tracking inputs "that invalidate the joint calling process if they
are changed". That is a *drift* concern: it presumes a dataset that is used more than once, so that a
second run can contradict the first. Every GVS use case has that property except the one this epic is
about — on TSPS the BigQuery dataset is created per run and deleted when the run finishes, so there is
never a second run to disagree with the first.

That does not make the work TSPS-irrelevant; it means two different needs are in play, and only one of
them is a TSPS dependency:

|                       | Deliverable A: record and export                                       | Deliverable B: enforce                                                     |
|-----------------------|------------------------------------------------------------------------|----------------------------------------------------------------------------|
| **Need**              | Explain, troubleshoot and cost-attribute a run after the fact          | Refuse to corrupt a dataset by changing an invalidating input between runs |
| **Driven by**         | TSPS, where the dataset is deleted and the record is all that survives | GVS Beta and AoU, where datasets persist and multi-batch ingest is normal  |
| **Belongs in**        | This epic                                                              | A Beta/AoU follow-on ticket                                                |
| **Testable on TSPS?** | Yes                                                                    | No — no TSPS run can ever exercise it                                      |

The recommendation is to build A now and specify B without building it. The schema is the same either
way: `scope` and `is_enforced` remain columns even though nothing consumes `is_enforced` in the first
release, so the follow-on adds behavior rather than migrating data.

The practical reason to split rather than build both: a drift check implemented under this epic would
ship without any TSPS run able to test it, against datasets that have no recorded history yet. It
would be exercised for the first time by a Beta user. Pairing it with a Beta/AoU ticket puts it where
it can be validated.

## Scope

`GvsJointVariantCalling.wdl` is the entry point of interest, as the top-level workflow TSPS will use.
GVS Beta uses the same entry point and so is covered incidentally.

AoU is kept in mind but is not the target: AoU callsets never run `GvsJointVariantCalling.wdl`, they
run `GvsBulkIngestGenomes.wdl`, `GvsCreateFilterSet.wdl`, `GvsExtractAvroFilesForHail.wdl` and
`GvsCreateVDS.wdl` individually. The schema below is consequently *workflow-agnostic* — one row per
workflow launch, one row per input, workflow name as a column — so those entry points can be wired in
later with no schema change, one task call each.

## Deliverable A: why recording inputs is worth doing for TSPS

TSPS executes its workflows in Terra, so the Terra run context that `GetToolVersions` already scrapes
(`workspace_id`, `workspace_bucket`, `submission_id`, `workflow_id`) is available; the workspace is
simply owned by the TSPS service rather than by the requester. Four properties of the TSPS model shape
what recording is for. (See `scripts/variantstore/docs/tsps/gvs-tsps-bigquery-requirements.md` on the
`vs_1941_bigquery_needs` branch for the fuller picture.)

1. **A delivered callset becomes unexplainable once the dataset is gone.** GVS today reverse-engineers
   its own configuration from data: `IsUsingCompressedReferences` (`GvsUtils.wdl:1084`) reads
   `INFORMATION_SCHEMA.COLUMNS` for `ref_ranges_001`, `IsVETS` (`GvsUtils.wdl:1014`) counts non-NULL
   `calibration_sensitivity` versus `vqslod` rows in `filter_set_info`, `GetExtractVetTableVersion`
   (`GvsUtils.wdl:1147`) does the same for vet table layout. On Beta and AoU that archaeology still
   works months later. On TSPS the dataset is deleted, the user never had BigQuery access, and only
   the VCFs remain — so a question like "which interval list, filter model and GATK version produced
   this?" has no answer at all unless it was recorded at launch. This is the one respect in which the
   ephemeral dataset makes the case *stronger* rather than weaker.

2. **Cost rows are hard to interpret without the inputs beside them.** `cost_observability` records
   what a run cost; sample count, `is_wgs`, interval list and scatter widths are what make a
   3x-more-expensive run explicable. On a service where every run is a billing event inside one shared
   project, that pairing is the difference between cost data and cost understanding.

3. **A service accumulates runs, and troubleshooting is comparative.** "This request behaved
   differently from the same request last month" is answerable only if both runs' pinned versions —
   `gvs_version`, `git_hash`, the four Docker images — were written down.

4. **Inputs are service-pinned, which makes the record small and stable.** TSPS pins nearly everything
   from its pipeline definition; the values that vary per request are the sample FOFN, WGS versus
   exome, the interval lists and the callset name. Note that WDL cannot distinguish "supplied by the
   user's request" from "pinned by the pipeline definition" — both arrive as inputs. If that
   distinction is wanted, TSPS has to supply it; `was_specified` below only separates "a value
   arrived" from "the WDL default applied".

## Record effective values, not declared defaults

Defaults differ between entry points by design — `drop_state` is `"FORTY"` in
`GvsJointVariantCalling.wdl:28` but `"NONE"` in `GvsBulkIngestGenomes.wdl:44`, and
`maximum_alternate_alleles` is `1000` in `GvsJointVariantCalling.wdl:83` but `100` in
`GvsExtractCallset.wdl:67` — because each suits how its own workflow is normally invoked, and the
top-level workflow threads its own values down to the subworkflows that need them.

The consequence for this proposal is that only *effective* values, resolved at the top level, are worth
recording; a declared default would be actively misleading about what a run did. The schema below
records effective values under the input's own name for that reason.

`drop_state` is the clearest illustration of why the record matters more than any default: the value AoU
requires is `"ZERO"`, which `AOU_DELIVERABLES.md:108` instructs an operator to type in and which no
workflow default supplies, so the correct value is determined by protocol documentation rather than by
code and nothing confirms it was followed. See `drop-state-handling.md` for the full picture, including
the follow-on tickets it suggests.

## Which inputs

Recording an input costs one row (~50 rows per launch). The expensive judgment call is not what to
store but what would be *enforced*, so everything is recorded and each input carries two attributes
that Deliverable B would later consume:

- `scope`: the blast radius over which the value must stay constant.
- `is_enforced`: whether a mismatch should stop the pipeline.

| Scope     | Meaning                                            | Consequence of a change                                                     |
|-----------|----------------------------------------------------|-----------------------------------------------------------------------------|
| `DATASET` | Must hold for the life of the BigQuery dataset     | Samples in the store are not comparable to each other                       |
| `CALLSET` | Must agree between filter set creation and extract | The callset is filtered or extracted inconsistently with how it was modeled |
| `RUN`     | Provenance only                                    | None; recorded for forensics, support and cost attribution                  |

On TSPS all three scopes collapse into the single run, and the classification is carried for the
record's own sake and for Deliverable B's benefit.

### `DATASET` scope

Each of these makes samples ingested before a change non-comparable to samples ingested after it, and
each is something a user would need in order to interpret a delivered callset.

| Input                                                                    | Consumed by                 | Why it matters                                                                                                                                                                                                |
|--------------------------------------------------------------------------|-----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `project_id`, `dataset_name`                                             | everything                  | Dataset identity; on TSPS the dataset name is also the per-run isolation boundary, so recording it gives an audit trail                                                                                       |
| `reference_name`                                                         | ingest, filter set, extract | Coordinates and reference bases                                                                                                                                                                               |
| `is_wgs`                                                                 | ingest, extract             | Selects WGS vs exome interval list and exome-specific extract behavior                                                                                                                                        |
| effective interval list (`interval_list`, or the `GetReference` default) | ingest, filter set, extract | Which regions samples cover; whether absence of data is meaningful                                                                                                                                            |
| `drop_state`                                                             | ingest, extract             | Which GQ reference blocks were discarded; determines reference genotype fidelity                                                                                                                              |
| `use_compressed_references`                                              | ingest                      | Physical `ref_ranges_%` schema (`location` vs `packed_ref_data`)                                                                                                                                              |
| `load_vet_and_ref_ranges`                                                | ingest                      | Which core tables were populated                                                                                                                                                                              |
| `load_vcf_headers`                                                       | ingest                      | Whether per-sample header provenance exists                                                                                                                                                                   |
| `use_parquet_ingest`, `parquet_output_gcs_dir`                           | ingest                      | Ingest path. TSPS is Parquet-only, so the flag should always be `true` — recording lets that be asserted rather than assumed, and the GCS directory outlives the dataset, so it matters for cleanup and audit |

### `CALLSET` scope

| Input                                                                  | Consumed by         | Why it matters                                                                                      |
|------------------------------------------------------------------------|---------------------|-----------------------------------------------------------------------------------------------------|
| effective `filter_set_name`                                            | filter set, extract | Which model extract applied                                                                         |
| `use_VQSR` (→ `use_VETS`)                                              | filter set          | VETS and VQSR scores are not interchangeable                                                        |
| `add_additional_annotations_to_sites_only_vcf`                         | filter set          | Model training features                                                                             |
| `INDEL_VQSR_max_gaussians_override`, `SNP_VQSR_max_gaussians_override` | filter set          | Model shape                                                                                         |
| `training_python_script`, `scoring_python_script`                      | filter set          | The VETS model *code* itself — the highest-leverage override we expose, and invisible in the output |
| effective `extract_table_prefix`                                       | prepare, extract    | Which prepare tables extract read                                                                   |
| `sample_names_to_extract`                                              | prepare             | Cohort definition                                                                                   |
| `extract_do_not_filter_override`                                       | extract             | Whether filters were applied at all                                                                 |
| `target_interval_list`                                                 | extract             | Exome/BGE target regions                                                                            |
| `maximum_alternate_alleles`                                            | extract             | Truncation of multi-allelic sites                                                                   |

### `RUN` scope

`call_set_identifier`, `extract_output_gcs_dir`, `extract_output_file_base_name`,
`merge_output_vcfs`, `bgzip_output_vcfs`, `collect_variant_calling_metrics`, `sample_set_name`,
`sample_id_column_name`, `vcf_files_column_name`, `vcf_index_files_column_name`,
`billing_project_id`, `query_labels`, `tighter_gcp_quotas`, `gatk_override`, and the performance knobs
(`extract_scatter_count`, `extract_preemptible_override`, `extract_maxretries_override`,
`load_data_scatter_width`, `load_data_preemptible_override`, `load_data_maxretries_override`,
`split_intervals_disk_size_override`, `split_intervals_mem_override`, `INDEL_VQSR_mem_gb_override`,
`SNP_VQSR_mem_gb_override`).

The performance knobs are `RUN` scope but are exactly what Deliverable A's cost-correlation use case
needs, which is a good illustration of why "not enforced" should not mean "not recorded".

`gatk_override` deserves a callout: a defined `gatk_override` means the run did not use the production
GATK binary. On TSPS that should never happen, so recording it makes `gatk_override IS NOT NULL` a
cheap assertion rather than an assumption.

### Proposed new input: `external_run_id`

An optional `String? external_run_id` on `GvsJointVariantCalling.wdl`, recorded in the run row, gives
TSPS somewhere to stamp its job/run ID so a delivered callset can be tied back to the request that
produced it. Because the Terra workspace is service-owned, workspace identity says nothing about the
requester, and `call_set_identifier` is user-supplied and mangled for naming purposes, so it is a poor
join key. Cheap to add now and awkward to retrofit.

### Non-input run context

| Field                                                                | Source                                                                                                                                  |
|----------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|
| workflow name                                                        | `GetToolVersions` — already present in the path it scrapes, see implementation notes                                                    |
| workflow ID (root), submission ID                                    | `GetToolVersions.workflow_id` / `.submission_id`                                                                                        |
| workspace ID, bucket, Google project                                 | `GetToolVersions.workspace_id` / `.workspace_bucket` / `.google_project`                                                                |
| workspace name, namespace                                            | Rawls, via the existing `scripts/get_workspace_name_for_import.py`                                                                      |
| BigQuery dataset location, workspace bucket location                 | `bq show` / `gcloud storage buckets describe`. The dataset location doubles as a check that TSPS runs land in the intended multi-region |
| launch timestamp                                                     | `GetToolVersions.date_as_secs_since_unix_epoch`                                                                                         |
| `gvs_version`, `git_branch_or_tag`, `git_hash`                       | `GetToolVersions`                                                                                                                       |
| `gatk_docker`, `variants_docker`, `cloud_sdk_docker`, `basic_docker` | `GetToolVersions`, effective values so caller overrides are captured                                                                    |

## Proposed schema

Two tables, created with the existing `GvsCreateTables.CreateTables` task with `partitioned = "false"`
and `superpartitioned = "false"`. Both are tiny.

### `gvs_workflow_run` — one row per workflow launch

| Column                            | Type      | Mode     | Notes                                                       |
|-----------------------------------|-----------|----------|-------------------------------------------------------------|
| `run_id`                          | STRING    | REQUIRED | Cromwell root workflow ID (UUID); primary key by convention |
| `workflow_name`                   | STRING    | REQUIRED | e.g. `GvsJointVariantCalling`                               |
| `external_run_id`                 | STRING    | NULLABLE | TSPS job/run ID when supplied                               |
| `call_set_identifier`             | STRING    | NULLABLE |                                                             |
| `started`                         | TIMESTAMP | REQUIRED | Launch time                                                 |
| `submission_id`                   | STRING    | REQUIRED | Terra submission ID                                         |
| `terra_workspace_name`            | STRING    | NULLABLE | Service-owned under TSPS                                    |
| `terra_workspace_namespace`       | STRING    | NULLABLE |                                                             |
| `terra_workspace_id`              | STRING    | NULLABLE |                                                             |
| `terra_workspace_bucket`          | STRING    | NULLABLE |                                                             |
| `terra_google_project`            | STRING    | NULLABLE |                                                             |
| `terra_workspace_bucket_location` | STRING    | NULLABLE |                                                             |
| `bq_dataset_location`             | STRING    | NULLABLE | Region/multi-region of the GVS dataset                      |
| `gvs_version`                     | STRING    | REQUIRED | `"unspecified"` for non-release runs                        |
| `git_branch_or_tag`               | STRING    | NULLABLE |                                                             |
| `git_hash`                        | STRING    | REQUIRED |                                                             |
| `gatk_docker`                     | STRING    | REQUIRED | Effective value                                             |
| `variants_docker`                 | STRING    | REQUIRED | Effective value                                             |
| `cloud_sdk_docker`                | STRING    | REQUIRED | Effective value                                             |
| `basic_docker`                    | STRING    | REQUIRED | Effective value                                             |
| `gatk_override`                   | STRING    | NULLABLE | Non-NULL means a non-production run                         |

```
[{"name":"run_id","type":"STRING","mode":"REQUIRED"},{"name":"workflow_name","type":"STRING","mode":"REQUIRED"},{"name":"external_run_id","type":"STRING","mode":"NULLABLE"},{"name":"call_set_identifier","type":"STRING","mode":"NULLABLE"},{"name":"started","type":"TIMESTAMP","mode":"REQUIRED"},{"name":"submission_id","type":"STRING","mode":"REQUIRED"},{"name":"terra_workspace_name","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_namespace","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_id","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_bucket","type":"STRING","mode":"NULLABLE"},{"name":"terra_google_project","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_bucket_location","type":"STRING","mode":"NULLABLE"},{"name":"bq_dataset_location","type":"STRING","mode":"NULLABLE"},{"name":"gvs_version","type":"STRING","mode":"REQUIRED"},{"name":"git_branch_or_tag","type":"STRING","mode":"NULLABLE"},{"name":"git_hash","type":"STRING","mode":"REQUIRED"},{"name":"gatk_docker","type":"STRING","mode":"REQUIRED"},{"name":"variants_docker","type":"STRING","mode":"REQUIRED"},{"name":"cloud_sdk_docker","type":"STRING","mode":"REQUIRED"},{"name":"basic_docker","type":"STRING","mode":"REQUIRED"},{"name":"gatk_override","type":"STRING","mode":"NULLABLE"}]
```

### `gvs_workflow_run_input` — one row per input per launch

| Column            | Type    | Mode     | Notes                                                                                                        |
|-------------------|---------|----------|--------------------------------------------------------------------------------------------------------------|
| `run_id`          | STRING  | REQUIRED | Joins to `gvs_workflow_run`                                                                                  |
| `input_name`      | STRING  | REQUIRED | WDL input name, e.g. `drop_state`; effective values recorded under the input's own name, not `effective_*`   |
| `effective_value` | STRING  | NULLABLE | Canonical string rendering; NULL means an optional input with no value                                       |
| `value_type`      | STRING  | REQUIRED | `STRING` \| `BOOLEAN` \| `INT` \| `FLOAT` \| `FILE` \| `ARRAY`                                               |
| `scope`           | STRING  | REQUIRED | `DATASET` \| `CALLSET` \| `RUN`                                                                              |
| `is_enforced`     | BOOLEAN | REQUIRED | Recorded now, consumed by Deliverable B                                                                      |
| `was_specified`   | BOOLEAN | NULLABLE | `false` when the WDL default applied — distinguishes "chose the default" from "the default changed under us" |
| `file_generation` | STRING  | NULLABLE | GCS object generation, for `FILE` values                                                                     |
| `file_md5`        | STRING  | NULLABLE | GCS `md5Hash`, for `FILE` values                                                                             |
| `file_size_bytes` | INTEGER | NULLABLE | For `FILE` values                                                                                            |

```
[{"name":"run_id","type":"STRING","mode":"REQUIRED"},{"name":"input_name","type":"STRING","mode":"REQUIRED"},{"name":"effective_value","type":"STRING","mode":"NULLABLE"},{"name":"value_type","type":"STRING","mode":"REQUIRED"},{"name":"scope","type":"STRING","mode":"REQUIRED"},{"name":"is_enforced","type":"BOOLEAN","mode":"REQUIRED"},{"name":"was_specified","type":"BOOLEAN","mode":"NULLABLE"},{"name":"file_generation","type":"STRING","mode":"NULLABLE"},{"name":"file_md5","type":"STRING","mode":"NULLABLE"},{"name":"file_size_bytes","type":"INTEGER","mode":"NULLABLE"}]
```

The GCS fingerprint columns are why `File` inputs need more than a path. An interval list, a
`sample_names_to_extract` file or a `training_python_script` can be overwritten in place at the same
URI, so the object generation is the only reliable answer to "is this the same interval list that run
used?" — and on TSPS the interval list is the input most likely to vary per request.

### Optional: `gvs_workflow_run_event`

The run row is written at launch, before any real work starts, so failed and aborted runs are visible
too — which matters for Deliverable A, since a failed run is exactly the one someone will ask about.
Terminal status, if wanted, is best appended rather than updated, mirroring the existing
`sample_load_status` pattern (`sample_id`, `status`, `event_timestamp`) and avoiding DML against rows
that may still be in the streaming buffer:

```
[{"name":"run_id","type":"STRING","mode":"REQUIRED"},{"name":"event","type":"STRING","mode":"REQUIRED"},{"name":"event_timestamp","type":"TIMESTAMP","mode":"REQUIRED"}]
```

### Why one row per input rather than one column per input

A wide table (one column per WDL input) is more pleasant to query but a poor fit:

- `GvsJointVariantCalling` alone has ~45 inputs and they change from release to release. A wide table
  needs a schema migration each time an input is added or renamed.
- Made workflow-agnostic, a wide table becomes the union of every entry point's inputs — hundreds of
  mostly-NULL columns.
- Deliverable B's checks can then be written once, generically, without enumerating column names.

The cost is stringly-typed values. Canonical rendering rules keep that manageable: booleans as
`true`/`false`, arrays as JSON array text, files as `gs://` URIs, NULL reserved for "no value". Views
can present a typed, pivoted shape where convenient.

## Where the record lives, given that TSPS deletes the dataset

Three parts, in priority order:

1. **Emit the rows as a workflow-level `File` output** (newline-delimited JSON). This is what survives
   teardown on TSPS and asks nothing of TSPS beyond retaining a small file — no BigQuery access, no
   schema knowledge. It is also the natural artifact to hand to a user asking how their callset was
   produced. For Deliverable A this is the essential piece.
2. **Write the two tables in the run's own dataset.** One code path for TSPS, Beta and later AoU, and
   for Beta the durable home. Also where Deliverable B would read from.
3. **Optionally copy to a long-lived GVS-owned dataset** in the shared TSPS project, enabling cross-run
   questions ("which interval lists have ever been run through this service?", "how does cost track
   with sample count across all runs?"). This is new infrastructure, not required by the AC, and is
   the piece that would also make Deliverable B possible on TSPS.

Whatever is chosen should match how `cost_observability` is extracted before teardown — the same
constraint, and ideally the same mechanism, so there is one answer to "get the run's records out
before the dataset disappears".

## How the values get recorded

A new `WriteWorkflowRunMetadata` task in `GvsUtils.wdl`, called from `GvsJointVariantCalling.wdl`
immediately after `GetToolVersions`/`GetReference` and before `GvsBulkIngestGenomes`, with its `done`
output gating ingest (the established GVS ordering idiom).

Implementation details that matter:

- **Declare the task's inputs as `String`, never `File`.** `interval_list`, `target_interval_list`,
  `sample_names_to_extract` and the VETS scripts are `File` at the workflow level; passing them to a
  `File` task input would localize them and record a container-local path. Coercing to `String` at the
  call site preserves the `gs://` URI and costs nothing.
- **Compute effective values in the workflow, not the task.** `GvsJointVariantCalling.wdl` already
  computes `effective_*` values for dockers, git hash, workspace and the extract names; the same
  pattern covers the rest, and `defined(x)` supplies `was_specified`.
- **Load with `bq load`, not DML.** Emit newline-delimited JSON and
  `bq --apilog=false load --source_format=NEWLINE_DELIMITED_JSON`, following the precedent at
  `GvsUtils.wdl:859`. No `INSERT` means no quoting hazards and no streaming-buffer interactions — and
  the same NDJSON file is the workflow output from part 1 above, so there is one artifact, not two.
- **File fingerprints** come from
  `gcloud storage objects describe <uri> --format="value(generation,md5_hash,size)"`, tolerating
  failure: a missing or inaccessible object should not fail the run, just record the URI with NULL
  fingerprints.
- **`GetToolVersions` needs only one addition, `workflow_name`.** The delocalization path it already
  scrapes is `gs://fc-<workspace id>/submissions/<submission id>/<workflow name>/<workflow id>/...`,
  and the existing regex at `GvsUtils.wdl:96` already captures the workflow name as group 4 — the new
  output is the same `sed` with `\4` instead of `\5`.
- **Workspace name/namespace and region lookups belong in the new task, not `GetToolVersions`.**
  `GetToolVersions` runs in the alpine `cloud-sdk` image, which has neither `requests` nor
  `terra_notebook_utils`. The new task can run in `variants_docker`, where
  `/app/get_workspace_name_for_import.py` already returns workspace name and namespace from a
  workspace ID; `bq show --format=prettyjson` and `gcloud storage buckets describe` supply the
  locations.

## Deliverable B, specified but not built: stop asking operators for facts we already know

Recorded here so the follow-on ticket has a starting point, and so it is clear what the columns above
are for.

The obvious framing for enforcement is drift detection — compare this run to the last one and complain.
The more valuable framing is narrower and more actionable: several GVS inputs are not choices at all,
they are *assertions about how the data was already ingested*, and an operator supplies them by hand
long after the fact. Those inputs should be looked up rather than asked for.

### Worked example: `drop_state` at extract

`drop_state` reaches extract as `--inferred-reference-state`, and
`ExtractCohortEngine.java:687` stamps that state's GQ onto every reference genotype extract synthesizes
for a sample with no row at a site. It is therefore not an independent knob — it is a claim about what
ingest dropped, and the delivered genotype qualities depend on the claim being true:

| Ingested with | Extract told     | Result                                                                                      |
|---------------|------------------|---------------------------------------------------------------------------------------------|
| `ZERO`        | `FORTY`          | GQ40 confident reference calls invented where no gVCF supported any confidence              |
| `FORTY`       | `ZERO`           | GQ0 where GQ40 was intended; conservative, still wrong                                      |
| anything      | argument omitted | GATK's own default `SIXTY` (`ExtractCohort.java:208`), the most confident value in the enum |

Nothing checks this today. Under `GvsJointVariantCalling.wdl` the two ends cannot diverge, because one
input feeds both. They can diverge for AoU, for sub-cohort extracts (`SUB_COHORT_WORKFLOW.md`), and for
any direct `GvsExtractCallset` run — exactly the cases where the value is typed in by hand, possibly
months after ingest, with nothing recording what ingest used.

With launch metadata recorded, extract should read the ingest value (`scope = 'DATASET'`,
`input_name = 'drop_state'`) instead of accepting it as an input, retaining the input as an explicit
override and failing when an override disagrees with the record. For datasets that predate the table
there is a serviceable fallback: which GQ band is wholly absent from `ref_ranges` reveals what was
dropped, the same archaeology as `IsUsingCompressedReferences`, and it need only run once per dataset.
Where `use_compressed_references` is set the state is packed into `packed_ref_data`, so that fallback
needs the unpack logic (`SchemaUtils#encodeCompressedRefBlock` / `UnpackRefRangeInfo`).

One detail worth preserving in any such change: `"NONE"` asserts the same thing at both ends — no band
was dropped — and extract remaps it to `"ZERO"` (`GvsExtractCallset.wdl:440`) because
`GQStateEnum.NONE` carries a null `referenceGQ` and extract needs a concrete GQ to stamp. That remap is
consistent with GVS's convention that absent reference data reads as GQ0, which `GQStateEnum.java:16-18`
documents directly.

This is not hypothetical. The AoU small callsets are extracted without specifying `drop_state` at all,
and they come out correct only because a hand-typed ingest value agrees with an extract default that was
set for an unrelated reason and then remapped in a third change. Nothing in the system knows those facts
agree, and because no artifact records what ingest used, the agreement can only be re-derived from git
history and the data. `SUB_COHORT_WORKFLOW.md:27` already states the rule — the invariant is known, it
just has nowhere to live. `drop-state-handling.md` documents the case in full.

### The same pattern applies to every value extract currently has to be told

`IsUsingCompressedReferences` (`GvsUtils.wdl:1084`), `IsVETS` (`GvsUtils.wdl:1014`) and
`GetExtractVetTableVersion` (`GvsUtils.wdl:1147`) all exist because a downstream workflow needs a fact
about ingest or filtering that nobody wrote down. Each becomes a metadata read, with the existing
inference retained as the fallback for datasets that predate the table. `reference_name` and
`use_compressed_references` are the same shape of problem.

### Drift detection, where a dataset really is reused

A `CheckDatasetInvariants` task, run before ingest, reads the most recent prior run against this
dataset, diffs its enforced `DATASET`-scope values against this run's, and fails with a readable diff.
An `allow_dataset_invariant_change` escape hatch covers deliberate changes and, being an input, is
itself recorded — so a deliberate change leaves a trail. Datasets with no prior recorded run must pass
rather than fail, which is also what makes the check a harmless no-op on TSPS.

Drift across a dataset's whole history is a single generic query:

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

### Input allowlists

The VS-1941 requirements doc raises the product questions of runs using unvalidated or mixed interval
lists. Recording the declared interval list is the prerequisite for the "unvalidated" half: a run whose
interval list is not on an allowlist of known-good lists can be rejected or flagged. The *mixed* half is
a property of the input gVCF headers rather than of a WDL input, so this metadata is the declared value
to validate those headers against, not a substitute for reading them.

## Open questions for review

1. **Does the epic owner agree with the A/B split?** Specifically, that the ticket's stated AC — inputs
   that invalidate joint calling if changed — describes a Beta/AoU need, and what this epic actually
   needs is the record-and-export half.
2. **Was the drift framing motivated by resumability?** VS-1962 covers runtime state. If TSPS ever
   resumes or retries a run against a dataset that outlived the first attempt, drift between the
   original launch and the resumed launch *is* the invalidation concern, and Deliverable B becomes a
   genuine TSPS dependency. This is the one story in which the AC as written is TSPS-motivated, and it
   is worth confirming before deferring B.
3. **Which retention mechanism does TSPS want?** Recommendation: the NDJSON workflow output as the
   surviving record, with the in-dataset tables as the canonical write. Is per-run retention on the
   TSPS side sufficient, or do we want the long-lived cross-run registry (which is also what would
   make Deliverable B possible on TSPS)?
4. **Should `external_run_id` be plumbed through now?** Without it, tying a callset back to a TSPS
   request depends on the Terra submission record in a service-owned workspace.
5. **When do the AoU entry points get wired in?** The schema supports them now; the work is one task
   call per workflow and can follow the TSPS work rather than gate it. Note that AoU is where the
   `drop_state` recording argument above bites hardest, since the required `"ZERO"` is supplied by an
   operator following a document rather than by any default.
