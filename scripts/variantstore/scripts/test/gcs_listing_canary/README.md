# GCS-listing canary (VS-1990)

`DiscoverParquetFiles` in `scripts/variantstore/wdl/GvsImportGenomes.wdl` lists the Parquet output
directory with `gcloud storage ls` and treats **exactly one** failure mode — stderr containing
`One or more URLs matched no objects` — as a legitimately empty directory (proceed with an empty
file list). Every other non-zero exit fails the task. This is what makes re-runs idempotent when a
sample produced no output, while still failing loudly on a real listing error.

That behavior is coupled to gcloud's error wording. `gcloud` is version-pinned in the Variants
Docker image, so the wording can only change at a deliberate cloud-sdk bump — a reviewed,
integration-tested event, not something that drifts under a running pipeline. This canary is the
automated "watch it": run it against a real `gcloud` and real GCS and it asserts, on live output,
that the empty-vs-error classification still holds.

## What it checks

A **preflight** runs first: it lists a known non-empty prefix. If that fails we cannot reach GCS
(no credentials, no network, or no access), so the empty-vs-error contract cannot be evaluated and
the canary **skips** (exit 77) rather than reporting a false drift. If the preflight succeeds, two
checks run, each using the same command shape as the WDL task (`gcloud storage ls --recursive
[--billing-project P] "<url>/"`):

1. **Empty prefix under an existing bucket** → non-zero exit **and** stderr matches the sentinel.
   This is the drift canary. If a cloud-sdk bump reworded the message, the WDL would stop
   recognizing the empty case. The WDL failure is fail-safe (a legitimately empty directory would
   fail the task loudly rather than silently load nothing), but this catches the drift earlier, in
   a test, and tells you to update both places.
2. **Nonexistent bucket (a genuine error)** → non-zero exit **and** stderr does **not** match the
   sentinel. This proves a real listing error stays distinguishable, so the WDL still fails fast on
   it — the direction that would otherwise be silent data loss.

## When it runs

`build_docker.sh` runs this canary **inside every freshly built Variants image**, after the unit
tests and before the push, so a cloud-sdk bump that reworded the message is caught at rebuild time
— which is the only time the pinned `gcloud` can change. It runs against the image's `gcloud` (not
the host's), using the building developer's mounted credentials; a drift (exit 1) aborts the build,
while an unreachable/unauthenticated environment (exit 77) only warns so an unrelated rebuild is not
blocked. You can also run it standalone as below.

## The coupling you must maintain

The `SENTINEL` in `run_gcs_listing_canary.sh` must stay byte-for-byte identical to the `grep`
pattern in `DiscoverParquetFiles`. If you change one, change the other — and re-run this canary.

## Prerequisites

- A real `gcloud` (the one whose wording you want to validate — normally the version pinned in the
  Variants image).
- Credentials that can list a GCS bucket. The defaults list a public Broad bucket
  (`gs://gcp-public-data--broad-references`), so an authenticated Broad account works out of the box.

## Usage

```shell
./run_gcs_listing_canary.sh
```

Optional env vars (defaults shown):

| Var               | Default                                           | Meaning                                                                      |
|-------------------|---------------------------------------------------|------------------------------------------------------------------------------|
| `BILLING_PROJECT` | `$(gcloud config get-value project)`              | Passed as `--billing-project`; if empty the flag is omitted (as in the WDL). |
| `EMPTY_BUCKET`    | `gs://gcp-public-data--broad-references`          | Any bucket you can list; a UUID subpath under it is the empty prefix.        |
| `NONEMPTY_URL`    | `gs://gcp-public-data--broad-references/hg38/v0/` | A known non-empty prefix for the happy-path check.                           |

Exit status: **0** when all checks pass, **1** when the contract has drifted, **77** when the canary
skipped because it could not reach GCS (no credentials, no network, or no access).
`run_gcs_listing_canary.txt` in this directory is a captured passing run for reference.

The canary creates nothing — it lists throwaway UUID paths that are guaranteed empty or nonexistent.
