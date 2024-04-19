usage() {
    echo "

USAGE: ./build_docker.sh

Generate a tag suitable for publication of a Variants Docker image.
Tags will be of the form <ISO 8601 Date>-alpine-<Docker image ID>.

e.g. 2024-04-18-alpine-f000ba44
"
    exit 1
}

if [[ $? -ne 0 ]]
then
    usage
fi

docker build .

TAG=$(python3 ./build_docker_tag.py $*)
# Take everything after the last slash to recover the Docker image ID from the tag.
IMAGE_ID=${TAG##*-}
BASE_REPO="broad-dsde-methods/variantstore"
REPO_WITH_TAG="${BASE_REPO}:${TAG}"
docker tag "${IMAGE_ID}" "${REPO_WITH_TAG}"
GCR_TAG="us.gcr.io/${REPO_WITH_TAG}"

# Run unit tests before pushing to GCR.
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

docker tag "${REPO_WITH_TAG}" "${GCR_TAG}"
docker push "${GCR_TAG}"

echo "Docker image pushed to \"${GCR_TAG}\""
