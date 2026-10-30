#!/usr/bin/env bash
#
# VS-1990 — canary for the "empty GCS listing is OK, everything else fails fast" logic in
# GvsImportGenomes.wdl's DiscoverParquetFiles task.
#
# That task treats exactly one gcloud failure mode -- stderr containing
#   "One or more URLs matched no objects"
# -- as a legitimately empty directory (proceed with an empty file list); every other non-zero
# exit fails the task. This couples GVS idempotency to gcloud's wording. gcloud is version-pinned
# in the Variants image, so the wording can only change at a deliberate cloud-sdk bump -- and a
# cloud-sdk bump forces a Variants image rebuild, so build_docker.sh runs this canary against the
# freshly built image's gcloud (see build_docker.sh). It also runs standalone. Given a real gcloud
# + real GCS it asserts, on live output, that
#
#   1. an EMPTY prefix under an existing bucket        -> non-zero exit AND stderr matches SENTINEL
#      (the drift canary: if a cloud-sdk bump reworded this, the WDL would stop recognizing the
#       empty case -- the failure is fail-safe/loud, but this catches it earlier, in a test)
#   2. a NONEXISTENT bucket (a genuine listing error)  -> non-zero exit AND stderr does NOT match
#      SENTINEL (proves a real error stays distinguishable, so the WDL still fails fast on it --
#      this is the direction that would otherwise be silent data loss)
#
# A PREFLIGHT lists a known non-empty prefix first. If that fails we cannot reach GCS (no
# credentials, no network, no access), so the empty-vs-error contract cannot be evaluated and the
# canary SKIPS rather than reporting a false drift.
#
# SENTINEL below MUST stay identical to the grep pattern in DiscoverParquetFiles
# (scripts/variantstore/wdl/GvsImportGenomes.wdl). If you change one, change the other.
#
# Exit status: 0 = contract holds; 1 = contract drifted; 77 = SKIPPED (GCS unreachable/unauthed).
#
# This is a MANUAL / build-time harness, NOT a CI unit test: it needs a real gcloud and credentials
# that can list a GCS bucket. It creates nothing (it lists throwaway UUID paths that are guaranteed
# empty / nonexistent) and cleans up nothing.
#
# Optional env vars (defaults shown):
#   BILLING_PROJECT   $(gcloud config get-value project)   # passed as --billing-project; if empty
#                                                           # the flag is omitted (matches the WDL,
#                                                           # where billing_project_id is optional)
#   EMPTY_BUCKET      gs://gcp-public-data--broad-references  # any bucket you can list; a UUID
#                                                           # subpath under it is the empty prefix
#   NONEMPTY_URL      gs://gcp-public-data--broad-references/hg38/v0/  # a known non-empty prefix
#
# Usage:
#   ./run_gcs_listing_canary.sh
#   BILLING_PROJECT=my-proj EMPTY_BUCKET=gs://my-bucket ./run_gcs_listing_canary.sh
#
set -o errexit -o nounset -o pipefail

# MUST match GvsImportGenomes.wdl DiscoverParquetFiles.
SENTINEL='One or more URLs matched no objects'

# Exit code for "environment cannot run the check" (automake's conventional SKIP code).
readonly SKIP_EXIT=77

BILLING_PROJECT=${BILLING_PROJECT:-$(gcloud config get-value project 2>/dev/null || true)}
EMPTY_BUCKET=${EMPTY_BUCKET:-gs://gcp-public-data--broad-references}
NONEMPTY_URL=${NONEMPTY_URL:-gs://gcp-public-data--broad-references/hg38/v0/}

BILLING_FLAG=()
if [[ -n "${BILLING_PROJECT}" ]]; then
  BILLING_FLAG=(--billing-project "${BILLING_PROJECT}")
fi

UUID=$(python3 -c 'import uuid; print(uuid.uuid4())')
EMPTY_URL="${EMPTY_BUCKET%/}/gvs-1990-canary-empty-${UUID}/"
MISSING_URL="gs://gvs-1990-canary-missing-${UUID}/some/dir/"

banner() { echo; echo "==================== $* ===================="; }

# Run the exact command shape DiscoverParquetFiles uses. Sets globals RC and STDERR_FILE.
run_ls() {  # $1 = url
  local url=$1
  STDERR_FILE=$(mktemp)
  set +o errexit
  gcloud storage ls --recursive "${BILLING_FLAG[@]}" "${url}" >/dev/null 2>"${STDERR_FILE}"
  RC=$?
  set -o errexit
}

FAILURES=0
pass() { echo "  PASS: $*"; }
fail() { echo "  FAIL: $*"; FAILURES=$((FAILURES + 1)); }

banner "gcloud version (for the GCS listing canary)"
gcloud --version | head -1
echo "billing project: ${BILLING_PROJECT:-<none>}"
echo "empty bucket:    ${EMPTY_BUCKET}"
echo "non-empty url:   ${NONEMPTY_URL}"

banner "PREFLIGHT: reach GCS by listing a known non-empty prefix"
run_ls "${NONEMPTY_URL}"
echo "  url=${NONEMPTY_URL}  exit=${RC}"
if [[ ${RC} -ne 0 ]]; then
  echo "  stderr: $(cat "${STDERR_FILE}")"
  echo "  SKIP: could not list a known non-empty prefix. This is an environment problem (no gcloud"
  echo "  credentials, no network, or no access) -- NOT a wording drift. The empty-vs-error contract"
  echo "  was not evaluated."
  exit "${SKIP_EXIT}"
fi
pass "reached GCS and listed a non-empty prefix"

banner "CASE 1: empty prefix under an existing bucket -> non-zero + SENTINEL"
run_ls "${EMPTY_URL}"
echo "  url=${EMPTY_URL}"
echo "  exit=${RC}  stderr: $(cat "${STDERR_FILE}")"
[[ ${RC} -ne 0 ]] && pass "non-zero exit as expected" || fail "expected non-zero exit, got 0"
if grep -qF "${SENTINEL}" "${STDERR_FILE}"; then
  pass "stderr still matches SENTINEL -- the WDL will treat this as an empty directory"
else
  fail "stderr NO LONGER matches SENTINEL ('${SENTINEL}') -- gcloud reworded the empty-listing"
  fail "message. DiscoverParquetFiles would now FAIL FAST on a legitimately empty directory."
  fail "Update the grep pattern in GvsImportGenomes.wdl AND the SENTINEL in this script."
fi

banner "CASE 2: nonexistent bucket (genuine error) -> non-zero + NO SENTINEL"
run_ls "${MISSING_URL}"
echo "  url=${MISSING_URL}"
echo "  exit=${RC}  stderr: $(cat "${STDERR_FILE}")"
[[ ${RC} -ne 0 ]] && pass "non-zero exit as expected" || fail "expected non-zero exit, got 0"
if grep -qF "${SENTINEL}" "${STDERR_FILE}"; then
  fail "a genuine listing error NOW matches SENTINEL -- DiscoverParquetFiles would mis-read this"
  fail "real failure as 'empty' and silently load nothing. This is the dangerous direction."
else
  pass "stderr does not match SENTINEL -- the WDL will fail fast on this real error"
fi

banner "RESULT"
if [[ ${FAILURES} -eq 0 ]]; then
  echo "ALL CHECKS PASSED -- the empty/error classification still holds on this gcloud."
  exit 0
else
  echo "${FAILURES} CHECK(S) FAILED -- see above. The gcloud listing contract has drifted."
  exit 1
fi
