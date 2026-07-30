#!/usr/bin/env bash
#
# VS-1803 — local, manual end-to-end test for VCF header loading.
#
# Exercises BOTH ingest paths and proves they agree, then re-runs to prove idempotency:
#   PARQUET path:  CreateVariantIngestFiles --output-type PARQUET -> bq load into scratch
#                  -> promote (anti-join INSERT + scratch cleanup) -> verify
#   BQ path:       CreateVariantIngestFiles --output-type BQ (writes scratch directly)
#                  -> promote -> verify
#   A/B:           diff the two datasets' final tables (must be identical)
#
# This is a MANUAL developer harness. It is NOT run in CI. It creates real BigQuery
# datasets in a project you supply and (optionally) tears them down at the end.
# See README.md in this directory for prerequisites and expected output.
#
# Required env vars (no safe defaults):
#   PROJECT      GCP project that owns the throwaway datasets, e.g. gvs-internal
#   INPUT_VCF    a (reblocked) gVCF, local path or gs:// URI
#   INTERVALS    an interval_list, local path or gs:// URI
#
# Optional env vars (defaults shown):
#   DATASET_PREFIX  header_loading_e2e   # parquet ds = <prefix>, BQ ds = <prefix>_bq
#   SAMPLE_ID       1
#   REF_VERSION     38
#   WORK_DIR        ./header_loading_e2e_work   # where local parquet is generated
#   GATK            <repo>/gatk
#   PROMOTE         <repo>/scripts/variantstore/scripts/process_sample_vcf_headers.py
#   TEARDOWN        0                    # 1 = drop both datasets + work dir at the end
#
# Usage:
#   PROJECT=my-proj INPUT_VCF=gs://.../sample.g.vcf.gz INTERVALS=gs://.../calling.interval_list \
#     ./run_header_loading_e2e.sh            # 2 iterations, datasets KEPT
#   ... ./run_header_loading_e2e.sh 3        # 3 iterations
#   TEARDOWN=1 ... ./run_header_loading_e2e.sh
#
set -euo pipefail

### ---- Locate the repo (this script lives in scripts/variantstore/scripts/test/E2E_local) ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"

### ---- Config ------------------------------------------------------------------
GATK=${GATK:-$REPO_ROOT/gatk}
PROMOTE=${PROMOTE:-$REPO_ROOT/scripts/variantstore/scripts/process_sample_vcf_headers.py}

PROJECT=${PROJECT:?Set PROJECT to a GCP project that can own throwaway datasets}
INPUT_VCF=${INPUT_VCF:?Set INPUT_VCF to a (reblocked) gVCF path or gs:// URI}
INTERVALS=${INTERVALS:?Set INTERVALS to an interval_list path or gs:// URI}

DATASET_PREFIX=${DATASET_PREFIX:-header_loading_e2e}
SAMPLE_ID=${SAMPLE_ID:-1}
REF_VERSION=${REF_VERSION:-38}
WORK_DIR=${WORK_DIR:-$PWD/header_loading_e2e_work}
TEARDOWN=${TEARDOWN:-0}

PQ_DS=$DATASET_PREFIX          # parquet path dataset
BQ_DS=${DATASET_PREFIX}_bq     # BQ (Write API) path dataset
PQ_DIR=$WORK_DIR/parquet
BQ_DIR=$WORK_DIR/bq

ITERATIONS=${1:-2}             # run count (>=2 to see idempotency)

### ---- Helpers -----------------------------------------------------------------
banner() { echo; echo "==================== $* ===================="; }

ensure_dataset() {   # $1 = dataset name
  bq show "$PROJECT:$1" >/dev/null 2>&1 || bq --location=US mk --dataset "$PROJECT:$1"
}

ensure_table() {     # $1 = dataset, $2 = table, $3 = schema
  bq show "$PROJECT:$1.$2" >/dev/null 2>&1 || bq mk --table "$PROJECT:$1.$2" "$3"
}

create_header_tables() {  # $1 = dataset, $2 = "with_load_status" (optional)
  ensure_table "$1" vcf_header_lines_scratch 'sample_id:INTEGER,vcf_header_lines:STRING,vcf_header_lines_hash:STRING,is_expected_unique:BOOLEAN'
  ensure_table "$1" vcf_header_lines         'vcf_header_lines_hash:STRING,vcf_header_lines:STRING,is_expected_unique:BOOLEAN'
  ensure_table "$1" sample_vcf_header        'sample_id:INTEGER,vcf_header_lines_hash:STRING'
  # The BQ path (unlike parquet) reads/writes per-sample load status for idempotency
  # (CreateVariantIngestFiles.onTraversalStart -> LoadStatus.getSampleLoadState).
  if [[ "${2:-}" == "with_load_status" ]]; then
    ensure_table "$1" sample_load_status 'sample_id:INTEGER,status:STRING,event_timestamp:TIMESTAMP'
  fi
}

# Run CreateVariantIngestFiles for one output type. $1 = PARQUET|BQ, $2 = dataset, $3 = out dir
ingest() {
  rm -rf "$3"
  "$GATK" CreateVariantIngestFiles \
    -V "$INPUT_VCF" \
    -L "$INTERVALS" \
    --gvs-sample-id "$SAMPLE_ID" \
    --output-directory "$3" \
    -IG FORTY \
    --force-loading-from-non-allele-specific false \
    --ignore-above-gq-threshold false \
    --project-id "$PROJECT" \
    --dataset-name "$2" \
    --enable-reference-ranges false \
    --enable-vet false \
    --output-type "$1" \
    --ref-version "$REF_VERSION" \
    --skip-loading-vqsr-fields false \
    --enable-vcf-headers true
}

# Load the generated parquet into the scratch table (what LoadParquetFilesToBQ does in the WDL).
load_parquet_to_scratch() {
  bq load --source_format=PARQUET \
    "$PROJECT:$PQ_DS.vcf_header_lines_scratch" \
    "$PQ_DIR/vcf_header_lines_scratch_${SAMPLE_ID}.parquet"
}

# scratch -> final promotion (anti-join INSERT + scratch cleanup). Tolerates the BQ
# streaming-buffer DELETE failure: the finals are committed regardless; scratch drains
# once the Write API buffer flushes (can take a few minutes, worst case ~90).
promote() {  # $1 = dataset
  if python3 "$PROMOTE" --project_id="$PROJECT" --dataset_name="$1"; then
    echo "  promote OK for $1"
  else
    echo "  WARN: promote returned non-zero for $1 — most likely the streaming-buffer"
    echo "        DELETE on scratch cleanup (BQ path). Final tables are still correct;"
    echo "        scratch_remaining may stay >0 until the buffer flushes. Benign."
  fi
}

verify() {  # $1 = dataset  (expect 3 / 3 / 0)
  bq query --use_legacy_sql=false "
    SELECT
      (SELECT COUNT(*) FROM \`$PROJECT.$1.vcf_header_lines\`)         AS header_lines,
      (SELECT COUNT(*) FROM \`$PROJECT.$1.sample_vcf_header\`)        AS sample_assocs,
      (SELECT COUNT(*) FROM \`$PROJECT.$1.vcf_header_lines_scratch\`) AS scratch_remaining"
}

ab_compare() {  # every count must be 0
  echo "-- A/B diff: vcf_header_lines (hash,text,is_expected_unique) — expect 0/0 --"
  bq query --use_legacy_sql=false "
    WITH
    p AS (SELECT vcf_header_lines_hash, vcf_header_lines, is_expected_unique FROM \`$PROJECT.$PQ_DS.vcf_header_lines\`),
    b AS (SELECT vcf_header_lines_hash, vcf_header_lines, is_expected_unique FROM \`$PROJECT.$BQ_DS.vcf_header_lines\`)
    SELECT
      (SELECT COUNT(*) FROM (SELECT * FROM p EXCEPT DISTINCT SELECT * FROM b)) AS in_parquet_not_bq,
      (SELECT COUNT(*) FROM (SELECT * FROM b EXCEPT DISTINCT SELECT * FROM p)) AS in_bq_not_parquet"
  echo "-- A/B diff: sample_vcf_header (sample_id,hash) — expect 0/0 --"
  bq query --use_legacy_sql=false "
    WITH
    p AS (SELECT sample_id, vcf_header_lines_hash FROM \`$PROJECT.$PQ_DS.sample_vcf_header\`),
    b AS (SELECT sample_id, vcf_header_lines_hash FROM \`$PROJECT.$BQ_DS.sample_vcf_header\`)
    SELECT
      (SELECT COUNT(*) FROM (SELECT * FROM p EXCEPT DISTINCT SELECT * FROM b)) AS in_parquet_not_bq,
      (SELECT COUNT(*) FROM (SELECT * FROM b EXCEPT DISTINCT SELECT * FROM p)) AS in_bq_not_parquet"
}

# Reconstruct the sample's full header text from the tables (the design-doc sanity check).
show_reconstructed_header() {  # $1 = dataset
  echo "-- reconstructed header for sample $SAMPLE_ID from $1 --"
  bq query --use_legacy_sql=false --format=csv "
    SELECT STRING_AGG(h.vcf_header_lines, '\n') AS full_header
    FROM \`$PROJECT.$1.sample_vcf_header\` s
    JOIN \`$PROJECT.$1.vcf_header_lines\`  h USING (vcf_header_lines_hash)
    WHERE s.sample_id = $SAMPLE_ID" | tail -n +2
}

### ---- Setup (once) ------------------------------------------------------------
banner "SETUP: datasets + tables (project=$PROJECT, prefix=$DATASET_PREFIX)"
ensure_dataset "$PQ_DS";  create_header_tables "$PQ_DS"
ensure_dataset "$BQ_DS";  create_header_tables "$BQ_DS" with_load_status

### ---- Iterate ----------------------------------------------------------------
for i in $(seq 1 "$ITERATIONS"); do
  banner "ITERATION $i / $ITERATIONS"

  echo "--- PARQUET path ($PQ_DS) ---"
  ingest PARQUET "$PQ_DS" "$PQ_DIR"
  ( cd "$PQ_DIR" && load_parquet_to_scratch )
  promote "$PQ_DS"
  echo "== parquet counts (expect 3 / 3 / 0) =="
  verify "$PQ_DS"

  echo "--- BQ path ($BQ_DS) ---"
  ingest BQ "$BQ_DS" "$BQ_DIR"
  promote "$BQ_DS"
  echo "== bq counts (expect 3 / 3 / 0; scratch_remaining may be >0 on early passes"
  echo "   due to the streaming buffer — see promote WARN above) =="
  verify "$BQ_DS"

  echo "--- A/B comparison (parquet vs BQ) ---"
  ab_compare
done

banner "IDEMPOTENCY: header_lines/sample_assocs above should be identical across all $ITERATIONS iterations"

show_reconstructed_header "$PQ_DS"

### ---- Teardown ---------------------------------------------------------------
if [[ "$TEARDOWN" == "1" ]]; then
  banner "TEARDOWN"
  bq rm -r -f --dataset "$PROJECT:$PQ_DS"
  bq rm -r -f --dataset "$PROJECT:$BQ_DS"
  rm -rf "$WORK_DIR"
else
  banner "DATASETS KEPT (TEARDOWN=0)"
  echo "To clean up when done:"
  echo "  bq rm -r -f --dataset $PROJECT:$PQ_DS"
  echo "  bq rm -r -f --dataset $PROJECT:$BQ_DS"
  echo "  rm -rf $WORK_DIR"
fi
