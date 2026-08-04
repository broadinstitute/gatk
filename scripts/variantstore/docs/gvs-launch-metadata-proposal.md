# Proposal: recording GVS launch metadata (VS-1961)

GVS currently does not record the input parameters it was launched with. This proposal identifies which inputs
matter and proposes a schema for storing them. It also separates two needs that a single piece of work cannot
serve at once.

## Two motivations, two deliverables

Tracking inputs "that invalidate the joint calling process if they are changed" is a *drift* concern: it
presumes a dataset used more than once, so that a second run can contradict the first. That describes Beta
and AoU rather than TSPS, where the BigQuery dataset is created per run and deleted when the run finishes,
leaving nothing for a later run to disagree with apart from TSPS-controlled retries. So two needs are in
play, and only one of them is a TSPS dependency:

|                       | Deliverable A: record and export                                       | Deliverable B: enforce                                                     |
|-----------------------|------------------------------------------------------------------------|----------------------------------------------------------------------------|
| **Need**              | Explain, troubleshoot and cost-attribute a run after the fact          | Refuse to corrupt a dataset by changing an invalidating input between runs |
| **Driven by**         | TSPS, where the dataset is deleted and the record is all that survives | GVS Beta and AoU, where datasets persist and multi-batch ingest is normal  |
| **Belongs in**        | This epic                                                              | A Beta/AoU follow-on ticket                                                |
| **Testable on TSPS?** | Yes                                                                    | No — no TSPS run can ever exercise it                                      |

This document will describe how to build deliverable A now in a way that should be compatible with building
deliverable B later.

## Scope

`GvsJointVariantCalling.wdl` is expected to be the TSPS entry point; GVS Beta uses this same WDL as its
entry point and so is covered incidentally.

AoU is kept in mind but is not the target: AoU callsets never run `GvsJointVariantCalling.wdl`, they
run `GvsBulkIngestGenomes.wdl`, `GvsPopulateAltAllele.wdl`, etc. individually. The schema below is
consequently *workflow-agnostic* — one row per workflow launch, one row per input, workflow name as
a column — so those entry points can be wired in later with no schema change, one task call each.

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
   this?" has no answer at all unless it was recorded at launch.

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

## Recording inputs

| Input                                                                    | Consumed by                  | Why it matters                                                                                                          |
|--------------------------------------------------------------------------|------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| `external_run_id` (new)                                                  | ingest                       | External job/run ID for auditing and traceability                                                                       |
| `project_id`, `dataset_name`                                             | everything                   | Dataset identity; on TSPS the dataset name is also the per-run isolation boundary, so recording it gives an audit trail |
| `bulk_ingest_fofn` (new)                                                 | ingest                       | Paths to the input data: VCFs and indexes                                                                               |
| `reference_name`                                                         | ingest, filter set, extract  | Coordinates and reference bases                                                                                         |
| `is_wgs`                                                                 | ingest, extract              | Selects WGS vs exome interval list and exome-specific extract behavior                                                  |
| effective interval list (`interval_list`, or the `GetReference` default) | ingest, filter set, extract  | Which regions samples cover; whether absence of data is meaningful                                                      |
| `drop_state`                                                             | ingest, extract              | Which GQ reference blocks were discarded; determines reference genotype fidelity                                        |
| `use_compressed_references`                                              | ingest, prepare, Avro export | Physical `ref_ranges_%` schema (`location` vs `packed_ref_data`); every reader of `ref_ranges` has to know which it is  |
| `extract_output_gcs_dir`                                                 | extract                      | Where the delivered VCFs were written; the first thing to look up when asked what produced a given callset              |

### Proposed new inputs

#### `bulk_ingest_fofn`

This optional String parameter would serve the same function as the parameter of the same name in
`GvsBulkIngestGenomes.wdl`: A FOFN of the VCF and VCF index data to load. We want this because:

- We really don't want to create and populate Terra data tables for every Joint Calling run in
  the shared TSPS workspace
- This is a low-cost option that reduces GVS's coupling to Terra and will serve this and many other
  use cases well

We would likely also want to capture the paths of the files enumerated within, as well as their hashes.
That per-sample provenance does not belong in `gvs_workflow_run_input`, whose grain is one row per
parameter — a TSPS run enumerates up to 10,000 samples. It belongs in the sample-level table the FOFN
contents are loaded into, keyed by `workflow_id`, with `gvs_workflow_run_input` recording only the FOFN
itself as a single row carrying its own generation and CRC32c, which pins the manifest exactly. AoU does
this loading by hand today; a canonical table should arrive with the ticket adding `bulk_ingest_fofn` to
`GvsJointVariantCalling.wdl`.

Collecting the hashes is cheap: CRC32c is present on every GCS object, and object listing returns it
alongside generation and size a thousand objects per page, so a whole delivery costs a handful of
requests rather than one per file.

#### `external_run_id`

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

| Column                            | Type      | Mode     | Notes                                                                                              |
|-----------------------------------|-----------|----------|----------------------------------------------------------------------------------------------------|
| `workflow_id`                     | STRING    | REQUIRED | Cromwell root workflow ID (UUID); primary key by convention                                        |
| `workflow_name`                   | STRING    | REQUIRED | e.g. `GvsJointVariantCalling`                                                                      |
| `external_run_id`                 | STRING    | NULLABLE | TSPS job/run ID when supplied                                                                      |
| `call_set_identifier`             | STRING    | NULLABLE |                                                                                                    |
| `started`                         | TIMESTAMP | REQUIRED | Launch time                                                                                        |
| `submission_id`                   | STRING    | REQUIRED | Terra submission ID                                                                                |
| `terra_workspace_name`            | STRING    | NULLABLE | Service-owned under TSPS                                                                           |
| `terra_workspace_namespace`       | STRING    | NULLABLE |                                                                                                    |
| `terra_workspace_id`              | STRING    | NULLABLE |                                                                                                    |
| `terra_workspace_bucket`          | STRING    | NULLABLE |                                                                                                    |
| `terra_google_project`            | STRING    | NULLABLE |                                                                                                    |
| `terra_workspace_bucket_location` | STRING    | NULLABLE |                                                                                                    |
| `bq_dataset_location`             | STRING    | NULLABLE | Region/multi-region of the GVS dataset                                                             |
| `gvs_version`                     | STRING    | REQUIRED | `"unspecified"` for non-release runs; TSPS should launch from a released tag so this is never that |
| `git_branch_or_tag`               | STRING    | NULLABLE |                                                                                                    |
| `git_hash`                        | STRING    | REQUIRED |                                                                                                    |
| `gatk_docker`                     | STRING    | REQUIRED | Effective value                                                                                    |
| `variants_docker`                 | STRING    | REQUIRED | Effective value                                                                                    |
| `cloud_sdk_docker`                | STRING    | REQUIRED | Effective value                                                                                    |
| `basic_docker`                    | STRING    | REQUIRED | Effective value                                                                                    |
| `gatk_override`                   | STRING    | NULLABLE | Non-NULL means a non-production run                                                                |

```
[{"name":"workflow_id","type":"STRING","mode":"REQUIRED"},{"name":"workflow_name","type":"STRING","mode":"REQUIRED"},{"name":"external_run_id","type":"STRING","mode":"NULLABLE"},{"name":"call_set_identifier","type":"STRING","mode":"NULLABLE"},{"name":"started","type":"TIMESTAMP","mode":"REQUIRED"},{"name":"submission_id","type":"STRING","mode":"REQUIRED"},{"name":"terra_workspace_name","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_namespace","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_id","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_bucket","type":"STRING","mode":"NULLABLE"},{"name":"terra_google_project","type":"STRING","mode":"NULLABLE"},{"name":"terra_workspace_bucket_location","type":"STRING","mode":"NULLABLE"},{"name":"bq_dataset_location","type":"STRING","mode":"NULLABLE"},{"name":"gvs_version","type":"STRING","mode":"REQUIRED"},{"name":"git_branch_or_tag","type":"STRING","mode":"NULLABLE"},{"name":"git_hash","type":"STRING","mode":"REQUIRED"},{"name":"gatk_docker","type":"STRING","mode":"REQUIRED"},{"name":"variants_docker","type":"STRING","mode":"REQUIRED"},{"name":"cloud_sdk_docker","type":"STRING","mode":"REQUIRED"},{"name":"basic_docker","type":"STRING","mode":"REQUIRED"},{"name":"gatk_override","type":"STRING","mode":"NULLABLE"}]
```

### `gvs_workflow_run_input` — one row per input per launch

| Column            | Type    | Mode     | Notes                                                                                                        |
|-------------------|---------|----------|--------------------------------------------------------------------------------------------------------------|
| `workflow_id`     | STRING  | REQUIRED | Joins to `gvs_workflow_run`                                                                                  |
| `input_name`      | STRING  | REQUIRED | WDL input name, e.g. `drop_state`; effective values recorded under the input's own name, not `effective_*`   |
| `effective_value` | STRING  | NULLABLE | Canonical string rendering; NULL means an optional input with no value                                       |
| `value_type`      | STRING  | REQUIRED | One of STRING, BOOLEAN, INT, FLOAT, FILE, ARRAY                                                              |
| `was_specified`   | BOOLEAN | NULLABLE | `false` when the WDL default applied — distinguishes "chose the default" from "the default changed under us" |
| `file_generation` | STRING  | NULLABLE | GCS object generation, for `FILE` values                                                                     |
| `file_crc32c`     | STRING  | NULLABLE | GCS `crc32c`, for `FILE` values. Preferred over `md5Hash`, which composite objects do not have               |
| `file_size_bytes` | INTEGER | NULLABLE | For `FILE` values                                                                                            |

```
[{"name":"workflow_id","type":"STRING","mode":"REQUIRED"},{"name":"input_name","type":"STRING","mode":"REQUIRED"},{"name":"effective_value","type":"STRING","mode":"NULLABLE"},{"name":"value_type","type":"STRING","mode":"REQUIRED"},{"name":"was_specified","type":"BOOLEAN","mode":"NULLABLE"},{"name":"file_generation","type":"STRING","mode":"NULLABLE"},{"name":"file_crc32c","type":"STRING","mode":"NULLABLE"},{"name":"file_size_bytes","type":"INTEGER","mode":"NULLABLE"}]
```

The GCS hash columns are why `File` inputs need more than a path. An interval list or a FOFN can
be overwritten in place at the same URI, so the object generation is the only reliable answer to "is this
the same interval list that run used?" — and on TSPS the interval list is the input most likely to vary
per request.

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

The sample-level table loaded from `bulk_ingest_fofn` should be exported the same way and at the same
time. It is the per-sample half of the provenance — which gVCFs went in, and their CRC32c hashes — and it
dies with the dataset like everything else. Exporting it also answers what sample count a run had, which
matters for the cost correlation above and which no launch input records: a TSPS run's 100 samples are a
row count, not a parameter.

## Restarts, and datasets that are reused

Outside TSPS the dataset persists, so a workflow can be re-submitted against it — most often after failing
partway through. A restart is a new Terra submission and therefore a new `workflow_id`, so it adds a row
rather than amending one. That is the intended behavior: two attempts did run, and the first attempt's
record survives. Attempts are correlated by the dataset they targeted, since the record lives in that
dataset: ordering `gvs_workflow_run` by `started` is the attempt history, and each row's `submission_id`
and `workflow_id` lead back to the corresponding Cromwell run. `call_set_identifier` is the human-facing
label for a callset, but it is free text, so a restart that mistypes it still correlates by dataset and
timestamp.

Nothing about resuming depends on this record, and nothing should. GVS re-derives what is already loaded
from the data: `CreateSampleDataViews` (`GvsImportGenomes.wdl:895`) builds `samples_with_reference_data`
and `samples_with_variant_data` from `INFORMATION_SCHEMA` partition ids `UNION DISTINCT`
`sample_load_status`, deliberately, because the Parquet flow writes no status rows while
`INFORMATION_SCHEMA` lags Storage Write API ingest by hours (`:911-921`); `CurateInputLists` (`:857`) then
trims the input lists to whatever `samples_with_all_data` (`:1041`) reports as missing. Because those views
read data rather than the `is_loaded` flag, a crash before that flag is set does not cause double-loading.

The distinction worth preserving is to infer facts about **data** and record facts about **configuration**.
Which samples are loaded is ground truth in BigQuery and self-correcting after a crash; `drop_state` and
the interval list leave no trace and cannot be re-derived. It is also why restarts are the likeliest source
of drift on a reusable dataset — an operator re-filling a form after a failure, or re-submitting after a
workflow version changed a default underneath them, can load the second batch under different invariants.
Two rows are what make that visible, which is the strongest practical argument for the enforcement
follow-on.

## Out of scope: runtime state (VS-1962)

This record describes a launch. It is written once, before work starts, and never updated. Tracking how far
a run got — which phases completed, which scattered shards succeeded, what a restart may skip — is a
different problem with a different write pattern, and belongs to VS-1962. Two decisions here are meant to
serve it. `gvs_workflow_run.workflow_id` is the intended parent key, so state rows can hang off a launch
without having to re-derive run identity. And whatever mechanism exports this record before a TSPS dataset
is torn down will have to carry state rows too — the constraint is identical, so it should not be solved
twice.

Worth carrying into that work: much of what looks like runtime state is already inferable from the data, as
the previous section describes. The useful output there is the list of things that genuinely are not.

## How the values get recorded

A new `WriteWorkflowRunMetadata` task in `GvsUtils.wdl`, called from `GvsJointVariantCalling.wdl`
immediately after `GetToolVersions`/`GetReference` and before `GvsBulkIngestGenomes`, with its `done`
output gating ingest (the established GVS ordering idiom).

Implementation details that matter:

- **Keep `File` types and set `localization_optional: true`.** `interval_list` and `bulk_ingest_fofn` are
  `File` at the workflow level, and the metadata task should declare them the same way under a
  `parameter_meta { <arg>: { localization_optional: true } }` entry. That skips localization while leaving
  the interpolated value as the `gs://` URI — the same mechanism that lets `CopyFile` run
  `gsutil cp ~{input_file}` (`GvsUtils.wdl:1635-1638`), and that `MergeVCFs` and `SplitIntervals` rely on.
  This is preferable to coercing to `String`, which would also carry the URI but would drop the input out
  of Cromwell's content-based call-cache hashing, so a file overwritten in place would no longer
  invalidate a cached result. Note that the `parameter_meta` entry is not merely a performance
  optimization: without it the task interpolates a localized container path, and the recorded value would
  silently become that path instead of the URI.
- **Compute effective values in the workflow, not the task.** `GvsJointVariantCalling.wdl` already
  computes `effective_*` values for dockers, git hash, workspace and the extract names; the same
  pattern covers the rest, and `defined(x)` supplies `was_specified`.
- **Load with `bq load`, not DML.** Emit newline-delimited JSON and
  `bq --apilog=false load --source_format=NEWLINE_DELIMITED_JSON`, following the precedent at
  `GvsUtils.wdl:859`. No `INSERT` means no quoting hazards and no streaming-buffer interactions — and
  the same NDJSON file is the workflow output from part 1 above, so there is one artifact, not two.
- **File hashes** come from
  `gcloud storage objects describe <uri> --format="value(generation,crc32c_hash,size)"`, and a failure
  should fail the run rather than record NULL hashes. If the interval list or the FOFN cannot be
  read, ingest is going to fail moments later with a worse message, so failing here is the better
  outcome: earlier, cheaper, and naming the object that could not be read. Recording NULLs instead would
  produce a record that looks complete while being unverifiable, which is the failure mode this proposal
  exists to prevent.
- Requester-pays inputs look unlikely to be supported on TSPS at all: TSPS collects no billing project
  from the user, so it would have to charge the egress to its own project, which sits badly with the TSPS
  credit model. If that is the decision, the cheapest place to enforce it is while reading file hashes —
  this task already reads each input object's metadata, so one `gcloud storage buckets describe` per
  distinct bucket can reject a requester-pays input at launch with a precise message.
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
