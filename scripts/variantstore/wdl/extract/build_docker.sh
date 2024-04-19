set -o errexit -o nounset -o pipefail

usage() {
    echo "

USAGE: ./build_docker.sh

Generate a tag suitable for publication of a Variants Docker image.
Tags will be of the form <ISO 8601 Date>-alpine-<Docker image ID>.

e.g. 2024-04-19-alpine-f000ba44
"
    exit 1
}

if [[ $# -ne 0 ]]
then
    usage
fi

docker build .

TAG=$(python3 ./build_docker_tag.py)
# Take everything after the last dash to recover the Docker image ID from the tag.
IMAGE_ID=${TAG##*-}
BASE_REPO="broad-dsde-methods/variantstore"
REPO_WITH_TAG="${BASE_REPO}/variantstore:${TAG}"
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
docker push "${GAR_TAG}"

echo "Docker image pushed to \"${GAR_TAG}\""
