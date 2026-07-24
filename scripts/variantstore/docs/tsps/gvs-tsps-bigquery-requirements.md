# BigQuery Requirements for GVS on TSPS

## Purpose and scope

This document describes the requirements for a minimum viable GVS on TSPS with
respect to Google BigQuery. This is purely a requirements-gathering document: it
identifies *what* needs to happen and *why*, and calls out gaps between
where GVS is today and where it needs to be for GVS on TSPS.

## Decisions needed from product/scientific ownership

There are potentially problematic input data conditions that GVS *could* detect,
assuming the required information is present in the input gVCFs. GVS already has
a manual BigQuery-based VCF header scanning process used for AoU which could be
adapted for these purposes, which is why these checks are discussed in this
document.

But whether GVS on TSPS *should* try to detect these conditions, and what it
should do if it does detect them, are product/scientific questions that need to
be answered:

1. Mixed DRAGEN versions within a callset
1. Unvalidated DRAGEN versions (i.e. versions that have not been tested and are
   not known to work with GVS)
1. Unsupported DRAGEN versions (i.e. versions known not to work with GVS, if
   any)
1. Mixed interval lists (e.g. WGS mixed with exome)
1. Unvalidated interval lists (e.g. T2T)

Some possible ways of handling these conditions might be:

1. Issuing an error message and terminating
2. Requiring the user to acknowledge the condition before proceeding
3. Something else

## Execution model contrast

|                                           | GVS Beta (today)                                                                                     | GVS on TSPS (target)                                                                                                                             |
|-------------------------------------------|------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| Who runs the workflow                     | The requesting user, in their own Terra workspace                                                    | The TSPS service account, on the user's behalf                                                                                                   |
| Who creates/owns the GCP project          | The user                                                                                             | N/A. Today, one GCP project shared by all TSPS requests (a round-robin pool of projects has been discussed but is not a current TSPS capability) |
| Who creates the BigQuery dataset          | The user, manually, before launching the workflow (`gvs-quickstart.md` Step 4)                       | TSPS, automatically, per pipeline run                                                                                                            |
| Who grants BigQuery IAM roles             | The user, manually, to their own Terra proxy group (`gvs-quickstart.md` Step 6)                      | TSPS, automatically, to its own service account                                                                                                  |
| Dataset lifetime                          | Persists indefinitely until the user deletes it                                                      | Should be created and destroyed automatically per run                                                                                            |
| Who sees/interacts with BigQuery directly | The user (dataset is visible in their GCP project)                                                   | Nobody. The dataset is an internal implementation detail                                                                                         |
| Failure recovery                          | User (or Variants team member) manually inspects BigQuery/Cromwell and decides whether/how to resume | TSPS must decide automatically whether to resume, and how                                                                                        |

Everything below expands on the requirements this contrast implies.

## Dataset lifecycle automation

Today, dataset creation is entirely a human, one-time, manual action
(`gvs-quickstart.md` Step 4: "Create a BigQuery dataset"). The only
programmatic dataset creation in the codebase is test/dev scaffolding, not
production automation.

**Requirements for TSPS:**

- Programmatically create a new BigQuery dataset per pipeline run (naming,
  location, labeling: see below).
- Programmatically and reliably delete the dataset when a run finishes
  successfully. This is new engineering, there is no reusable
  automated-deletion code path for BigQuery datasets anywhere in the GVS
  codebase today.
- **The shared-project model.** TSPS runs every request inside a single
  pre-provisioned GCP project shared across all requests; there is no
  per-request project. Each run's dataset is short-lived, created and destroyed
  within that shared project. GVS on TSPS can therefore never assume it has the
  project to itself; its job is simply to work inside whatever project it is
  handed. (TSPS has discussed replacing the single project with a round-robin
  pool of pre-prepared projects. That is not a current capability, and it would
  not change anything here: a run would still share whichever pool project it
  landed in.)
- Because isolation between concurrent runs now lives entirely at the
  **dataset** level rather than the project level, dataset naming must be
  collision-proof (e.g. keyed by TSPS run/job ID) and the code path that
  threads a run's dataset name through every BigQuery call needs to be
  correct; a bug that points one run's query at another run's dataset name
  is a cross-tenant data exposure, not just a wasted query, since the
  TSPS service account's project-level IAM grant gives it access to every
  dataset in the project.
- Dataset location: GVS Beta docs recommend "Multi-region" location to avoid
  restrictive quotas (`gvs-quickstart.md` Step 4). TSPS should carry this
  recommendation forward.
- No dataset labeling scheme (of the kind existing test automation uses:
  `service`, `team`, `environment`, `managedby` in `GvsUtils.wdl`) is needed
  purely for BigQuery cost attribution: because each pipeline run already gets
  its own dataset, the dataset name itself is the join key for per-run BigQuery
  cost (e.g. via `INFORMATION_SCHEMA.JOBS`/billing export filtered by
  destination dataset).

## IAM / permissions for the TSPS service account

The Beta docs are explicit about what a human must manually grant to their
Terra proxy group (`gvs-quickstart.md` Step 6, `TERRA_QUICKSTART.md`):

- **BigQuery Data Editor**
- **BigQuery Job User**
- **BigQuery Read Session User**

These are granted at the **project** level. Under the shared-project model,
this grant belongs to the TSPS service account and is part of the shared
project's one-time setup — not something done per run, nor something the GVS on
TSPS pipeline sets up or verifies at runtime. That makes it simpler to reproduce
than Beta's per-user manual grant.

It also means the grant is a standing, project-wide capability rather than
something scoped per run: the service account can read and write every dataset
in the project for its whole lifetime. Per-run isolation therefore has to come
from correct dataset-name handling in code, not from IAM boundaries (see "The
shared-project model" above).

There is no existing GVS code to reuse here, since Beta assumes a human performs
the grant once per project via the Google Cloud Console or `gcloud` CLI
(`gvs-quickstart.md` Step 6). The shared project needs this binding in place
before it is ever used to run a pipeline.

## BigQuery quotas and scale limits

**GVS on TSPS will use Parquet ingest exclusively**. This significantly
narrows the quota surface TSPS needs to worry about:

- `GvsImportGenomes.wdl`'s `is_rate_limited_beta_customer` flag exists
  specifically to stay on the good side of a particular **Storage Write API**
  quota. Since Storage Write API ingest will not be used in GVS on TSPS, this
  specific quota concern goes away.
- The quota regime that does apply to TSPS's Parquet-based load path is
  ordinary BigQuery load job quotas (e.g. jobs per project per day, concurrent
  load jobs), which are far more generous than the Storage Write API's
  throughput limits. `load_parquet_to_bq.py` already batches files per load
  job (default 10,000 files/job) and retries with exponential backoff
  on quota/rate-limit errors, so day-to-day handling of this quota is already
  built. Because the project is shared, concurrent runs draw on the *same*
  project-level quotas (load jobs, concurrent queries, etc.), so quota
  contention is an ongoing operational concern, not a hypothetical one.
- GVS has long utilized on-demand slots exclusively for its BigQuery compute
  capacity, as opposed to flat-rate/flex-slot reservations. It is currently
  expected that this will remain the case for GVS on TSPS, particularly given
  the credits billing model used by TSPS.

## Recoverability and resumability

GVS today has only *manual* recovery (re-run `GvsJointVariantCalling.wdl` with
call caching), which is not automated and not always sufficient. TSPS will need
automatic decide-and-resume behavior (see the "Failure recovery" row above).
This is being explored under VS-1962 ("Spike: GVS Metadata, storing runtime
state"), but at this stage it appears to be a workflow-state problem rather than
a BigQuery-specific one, so it is out of scope for this document beyond this
pointer.

There does not yet appear to be a follow-on ticket for *automatically* resuming
a recoverable failed run.

## Cost observability

- The BigQuery-side cost data (`cost_observability` table) should be
  reusable as-is, but it only exists **inside** the per-run dataset. Since
  that dataset will be deleted once a run completes, TSPS should extract this
  table's contents out to wherever per-run records are stored *before* deleting
  the dataset.
- The same applies to any BigQuery-native job/billing data TSPS might want
  beyond what `cost_observability` captures (e.g. bytes processed or slot
  time from `INFORMATION_SCHEMA.JOBS` for that dataset/project); that too
  only exists as long as the underlying jobs are queryable, so TSPS should
  decide what it wants to retain and pull it before teardown.
- The compute-cost portion of `GvsCallsetCost.wdl` likely needs a TSPS-specific
  equivalent rather than being reused directly, since it's currently coupled
  to Terra workspace/submission semantics and a single tenant project model.
- Since the user never sees the underlying BQ project/dataset, TSPS is the
  only place cost data can be surfaced. It need not be shown to the *user*, but
  TSPS may want to record it internally.
