# Local End-to-End Header-Loading Test (VS-1803)

`run_header_loading_e2e.sh` is a **manual developer harness** that validates VCF
header loading across both GVS ingest paths - **Parquet** and **BQ**. It is **not run in CI** — it creates
real BigQuery datasets and requires authenticated `bq`/`gcloud` access, so it must
be run by hand against a throwaway project.

## What it validates

For each iteration it runs the full flow on both paths and compares them:

| Path | Flow |
|---|---|
| **Parquet** | `CreateVariantIngestFiles --output-type PARQUET` → `bq load` the local parquet into `vcf_header_lines_scratch` → promote (`process_sample_vcf_headers.py`) → verify |
| **BQ (Write API)** | `CreateVariantIngestFiles --output-type BQ` (writes scratch directly) → promote → verify |

It then:

- **A/B compares** the two datasets' final tables (`vcf_header_lines` and
  `sample_vcf_header`) — they must be byte-for-byte identical.
- **Repeats** the whole flow (2× by default) to prove **idempotency**: the final
  table counts must not change on re-runs. The parquet path stays put because the
  promotion is an anti-join INSERT; the BQ path stays put because
  `CreateVariantIngestFiles` short-circuits via `HEADERS_LOADED` in
  `sample_load_status`.

## Prerequisites

- Authenticated `gcloud` / `bq` (`gcloud auth login`, `gcloud auth application-default login`).
- A **throwaway** GCP project you can create and drop datasets in.
- A built `gatk` (`./gradlew clean shadowJar`) — or point `GATK` at one.
- `python3` on `PATH` for the promotion script.

## Running

Three env vars are required; the rest have defaults.

```bash
PROJECT=my-throwaway-project \
INPUT_VCF=gs://.../sample.reblocked.g.vcf.gz \
INTERVALS=gs://.../wgs_calling_regions.hg38.interval_list \
  ./run_header_loading_e2e.sh          # 2 iterations, datasets KEPT for inspection
```

Optional overrides (defaults in parentheses):

| Var | Default | Purpose |
|---|---|---|
| `DATASET_PREFIX` | `header_loading_e2e` | parquet ds = `<prefix>`, BQ ds = `<prefix>_bq` |
| `SAMPLE_ID` | `1` | gVCS sample id used for ingest |
| `REF_VERSION` | `38` | reference version passed to the tool |
| `WORK_DIR` | `./header_loading_e2e_work` | where local parquet is generated |
| `GATK` | `<repo>/gatk` | path to the gatk launcher |
| `PROMOTE` | `<repo>/scripts/variantstore/scripts/process_sample_vcf_headers.py` | promotion script |
| `TEARDOWN` | `0` | set `1` to drop both datasets + work dir at the end |

Pass an iteration count as the first argument (e.g. `./run_header_loading_e2e.sh 3`).

## Expected output

- Each `verify` reports **`3 / 3 / 0`** (`header_lines` / `sample_assocs` /
  `scratch_remaining`) and now **hard-fails** if the counts are off. These are
  chunk/row counts, not VCF-header-line counts — each `vcf_header_lines` row is a
  comma-joined blob of many header lines. Override the expected count with
  `EXPECT_HEADER_CHUNKS` for a different input VCF. On the **BQ path** only,
  `scratch_remaining > 0` is tolerated (streaming-buffer wrinkle, below); the
  **parquet path** must drain to `0`.
- Each `ab_compare` prints **`0 / 0`** for both `vcf_header_lines` and
  `sample_vcf_header`.
- The final step reconstructs the sample's full header from
  `sample_vcf_header` ⋈ `vcf_header_lines` and asserts it is **non-empty** (proves
  the split is losslessly reversible; exact equality to the BQ path is already
  covered by `ab_compare`).
- The counts are identical across every iteration.

### Known benign wrinkle: BQ scratch_remaining

On the **BQ path**, the scratch cleanup `DELETE` can fail against rows still in
BigQuery's streaming buffer ("would affect rows in the streaming buffer"). The
script tolerates this (prints a `WARN` and continues) because the final tables are
already committed — only `scratch_remaining` may read `> 0` until the buffer
flushes (a few minutes, worst case ~90). The parquet path uses batch `bq load`, so
it is unaffected and drains to `0` immediately.

## Cleanup

If run with `TEARDOWN=0` (the default), the script prints the exact `bq rm` /
`rm -rf` commands to remove the two datasets and the work dir when you're done.
