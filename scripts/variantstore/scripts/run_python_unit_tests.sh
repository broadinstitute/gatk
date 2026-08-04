#!/usr/bin/env bash
#
# Runs the GVS Python unit tests (test/test_*.py) together with a throwaway BigQuery emulator so the
# header-load integration test (test_process_sample_vcf_headers.py::TestHeaderLoadIntegration) runs
# hermetically -- no GCP project or secrets required.
#
# Two execution modes so the emulator wiring and test loop live in one place, shared by both callers:
#
#   ./run_python_unit_tests.sh <runtime_docker_image>
#       DOCKER MODE. Runs each test file inside the given image (which supplies python3 + the
#       dependencies); the working-copy source is bind-mounted at /in. Used by build_docker.sh with
#       the freshly built Variants image. Test containers join the emulator's Docker network and
#       reach it by container name.
#
#   ./run_python_unit_tests.sh
#       NATIVE MODE. Runs each test file with the host's python3 -- dependencies must already be
#       installed (e.g. by a CI "install dependencies" step). The emulator's port is published to
#       localhost. Used by the GitHub Actions workflow (VS-1983).
#
# Optional environment overrides:
#   EMULATOR_IMAGE   BigQuery emulator image (default ghcr.io/goccy/bigquery-emulator:latest).
#   EMULATOR_PORT    Host port for the emulator in native mode (default 9050).
#   PYTHON_BIN       Python interpreter for native mode (default python3).
#
# Exit status is non-zero if any test file fails or the emulator does not start.

set -o errexit -o nounset -o pipefail

RUNTIME_IMAGE="${1:-}"
if [[ -n "${RUNTIME_IMAGE}" ]]; then
    MODE="docker"
else
    MODE="native"
fi

# Resolve to this script's directory (the scripts dir) so the test glob, bind mount, and relative
# data-file paths are stable regardless of the caller's working directory.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# TODO(VS-1983): pin to a specific emulator release once this stabilizes.
EMULATOR_IMAGE="${EMULATOR_IMAGE:-ghcr.io/goccy/bigquery-emulator:latest}"
EMULATOR_PORT="${EMULATOR_PORT:-9050}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
# $$-suffixed so concurrent runs (e.g. a local build while CI runs) do not collide.
EMULATOR_NAME="gvs-bq-emulator-$$"
EMULATOR_NETWORK="gvs-python-test-net-$$"
EMULATOR_PROJECT="gvs-header-emulator"
EMULATOR_DATASET="header_it"

cleanup_emulator() {
    docker rm -f "${EMULATOR_NAME}" > /dev/null 2>&1 || true
    if [[ "${MODE}" == "docker" ]]; then
        docker network rm "${EMULATOR_NETWORK}" > /dev/null 2>&1 || true
    fi
}
trap cleanup_emulator EXIT

# A plain `docker run` (not an Actions `services:` block) is used because we need to pass
# --project/--dataset args, which service containers can't override. Docker mode puts the emulator on
# a private network (test containers join it and reach it by name); native mode publishes the port so
# the host python client can reach it at localhost.
if [[ "${MODE}" == "docker" ]]; then
    docker network create "${EMULATOR_NETWORK}" > /dev/null
    docker run -d --platform linux/amd64 --name "${EMULATOR_NAME}" --network "${EMULATOR_NETWORK}" \
        "${EMULATOR_IMAGE}" --project="${EMULATOR_PROJECT}" --dataset="${EMULATOR_DATASET}" > /dev/null
    EMULATOR_ENDPOINT="http://${EMULATOR_NAME}:9050"
else
    docker run -d --platform linux/amd64 --name "${EMULATOR_NAME}" -p "${EMULATOR_PORT}:9050" \
        "${EMULATOR_IMAGE}" --project="${EMULATOR_PROJECT}" --dataset="${EMULATOR_DATASET}" > /dev/null
    EMULATOR_ENDPOINT="http://localhost:${EMULATOR_PORT}"
fi

echo "Waiting for the BigQuery emulator to start..."
emulator_ready=0
for _ in $(seq 1 60)
do
    if docker logs "${EMULATOR_NAME}" 2>&1 | grep -q "REST server listening"; then
        emulator_ready=1
        break
    fi
    sleep 1
done
if [[ "${emulator_ready}" -ne 1 ]]; then
    echo "BigQuery emulator failed to start:" >&2
    docker logs "${EMULATOR_NAME}" >&2 || true
    exit 1
fi
echo "BigQuery emulator is up (mode=${MODE}, endpoint=${EMULATOR_ENDPOINT})."

# Run each test file. The GVS_HEADER_IT_* env vars are consumed only by the header-load integration
# test; other test files ignore them.
fail=0
for test in test/test_*.py
do
    set +o errexit
    if [[ "${MODE}" == "docker" ]]; then
        docker run --platform linux/amd64 --rm -v "${SCRIPT_DIR}":/in \
            --network "${EMULATOR_NETWORK}" \
            -e GVS_HEADER_IT_PROJECT="${EMULATOR_PROJECT}" \
            -e GVS_HEADER_IT_DATASET="${EMULATOR_DATASET}" \
            -e GVS_HEADER_IT_ENDPOINT="${EMULATOR_ENDPOINT}" \
            -t "${RUNTIME_IMAGE}" bash -c "cd /in; python3 -m unittest ${test}"
        rc=$?
    else
        # Invoke by bare module name with the test dir on PYTHONPATH so the stdlib `test` package
        # does not shadow our local test/ directory.
        module="$(basename "${test}" .py)"
        GVS_HEADER_IT_PROJECT="${EMULATOR_PROJECT}" \
        GVS_HEADER_IT_DATASET="${EMULATOR_DATASET}" \
        GVS_HEADER_IT_ENDPOINT="${EMULATOR_ENDPOINT}" \
        PYTHONPATH="${SCRIPT_DIR}:${SCRIPT_DIR}/test" "${PYTHON_BIN}" -m unittest "${module}"
        rc=$?
    fi
    set -o errexit
    if [[ ${rc} -ne 0 ]]; then
        fail=1
        echo "${test} has failed"
    fi
done

if [[ ${fail} -ne 0 ]]; then
    echo "One or more Python unit tests have failed."
    exit 1
fi
echo "All Python unit tests passed."
