if [ $# -lt 1 ]; then
    echo "USAGE: ./build_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g.: ./build_docker.sh mybranch_20210403"
    exit 1
fi

set -o xtrace

BASE_REPO="broad-dsde-methods/variantstore"
REPO_WITH_TAG="${BASE_REPO}:${1}"
GCR_TAG="us.gcr.io/${REPO_WITH_TAG}"

docker build . -t "${REPO_WITH_TAG}"

# Run unit tests before pushing to GCR.
fail=0
for test in test_*.py
do
    docker run -v "$PWD":/in -t "${REPO_WITH_TAG}" bash -c "cd /in; python3 -m unittest $test"
    if [ $? -ne 0 ]; then
        fail=1
        echo "$test has failed"
    fi
done

if [ $fail -ne 0 ]; then
    echo "One or more unit test has failed, exiting."
    exit $fail
fi

docker tag "${REPO_WITH_TAG}" "${GCR_TAG}"
docker push "${GCR_TAG}"

echo "docker image pushed to \"${GCR_TAG}\""
