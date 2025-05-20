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

docker buildx build --platform linux/amd64 --load .
IMAGE_ID=$(docker image ls -q | head -n 1)

# The Variants Docker image is alpine-based.
IMAGE_TYPE="alpine"

# Build the image tag using the image type and Docker image ID:
TAG=$(python3 ./build_docker_tag.py --image-id "${IMAGE_ID}" --image-type "${IMAGE_TYPE}")

BASE_REPO="broad-dsde-methods/gvs"
REPO_WITH_TAG="${BASE_REPO}/variants:${TAG}"
docker tag "${IMAGE_ID}" "${REPO_WITH_TAG}"

# Run unit tests before pushing.
set +o errexit
fail=0
for test in test/test_*.py
do
    docker run --rm -v "$PWD":/in -t "${REPO_WITH_TAG}" bash -c "cd /in; python3 -m unittest $test"
    if [ $? -ne 0 ]; then
        fail=1
        echo "$test has failed"
    fi
done

if [ $fail -ne 0 ]; then
    echo "One or more unit tests have failed, exiting."
    exit $fail
fi

set -o errexit

GAR_TAG="us-central1-docker.pkg.dev/${REPO_WITH_TAG}"
docker tag "${REPO_WITH_TAG}" "${GAR_TAG}"

# Docker must be configured for GAR before pushes will work:
# gcloud auth configure-docker us-central1-docker.pkg.dev
docker push "${GAR_TAG}"

echo "Docker image pushed to \"${GAR_TAG}\""
