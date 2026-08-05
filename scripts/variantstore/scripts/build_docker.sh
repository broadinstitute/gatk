set -o errexit -o nounset -o pipefail

usage() {
    echo "

USAGE: ./build_docker.sh

Build a Variants Docker image with an appropriate tag and push to GAR.
The repo name will be 'us-central1-docker.pkg.dev/broad-dsde-methods/gvs' and the image name will be 'variants'.
Tags will be of the form <ISO 8601 Date>-alpine-<Docker image ID>.

e.g. 2024-04-22-alpine-f000ba44
"
    exit 1
}

if [[ $# -ne 0 ]]
then
    usage
fi

set +o errexit
docker buildx ls | grep singlebuilder
RET_VAL=$?
set -o errexit
if [ $RET_VAL -eq 1 ]; then
  docker buildx create --platform linux/amd64 --name singlebuilder
  echo "Created the single platform builder"
fi

docker buildx use singlebuilder
echo "Using the single platform builder"

# Build and write the full image id to an id file.
docker buildx build --platform linux/amd64 --load . --iidfile idfile.txt

# Scrape out the first 12 hex digits of the SHA256 hash to use for the tag.
IMAGE_ID=$(cut -d : -f 2 idfile.txt | cut -c 1-12)

rm idfile.txt

# The Variants Docker image is alpine-based.
IMAGE_TYPE="alpine"

# Build the image tag using the image type and Docker image ID:
TAG=$(python3 ./build_docker_tag.py --image-id "${IMAGE_ID}" --image-type "${IMAGE_TYPE}")

BASE_REPO="broad-dsde-methods/gvs"
REPO_WITH_TAG="${BASE_REPO}/variants:${TAG}"
docker tag "${IMAGE_ID}" "${REPO_WITH_TAG}"

# Run the Python unit tests (including the BigQuery-emulator-backed header-load integration test)
# inside the freshly built image before pushing. The runner also manages the emulator lifecycle and
# is shared with CI (see run_python_unit_tests.sh); errexit aborts the build here if any test fails.
./run_python_unit_tests.sh "${REPO_WITH_TAG}"

GAR_TAG="us-central1-docker.pkg.dev/${REPO_WITH_TAG}"
docker tag "${REPO_WITH_TAG}" "${GAR_TAG}"

# Docker must be configured for GAR before pushes will work:
# gcloud auth configure-docker us-central1-docker.pkg.dev
docker push "${GAR_TAG}"

echo "Docker image pushed to \"${GAR_TAG}\""
