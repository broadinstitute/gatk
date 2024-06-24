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

# Write the full Docker image ID to a file. This will look something like:
# sha256:5286e46648c595295dcb58a4cc2ec0b8893c9f26d7d49393908e5ae6d4dea188
docker build . --iidfile idfile.txt
FULL_IMAGE_ID=$(cat idfile.txt)

# Take the slice of this full Docker image ID that corresponds with the output of `docker images`:
IMAGE_ID=${FULL_IMAGE_ID:7:12}

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
for test in test_*.py
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
#docker push "${GAR_TAG}"

echo "Docker image NOT pushed to \"${GAR_TAG}\""

rm idfile.txt
