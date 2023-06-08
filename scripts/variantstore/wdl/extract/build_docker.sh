usage() {
     echo "

USAGE: ./build_docker.sh [--release|--branch]

One of the options --release or --branch must be given, which determines the format of the Docker image tag.
Release tags look like <ISO 8601 Date>-alpine, branch tags look like <ISO 8601 Date>-alpine-<short git hash>.

e.g. 2023-06-06-alpine or 2023-06-06-alpine-ed338e48e
"
     exit 1
}


VALID_ARGS=$(getopt --options r,b --longoptions release,branch -- "$@")
if [[ $? -ne 0 ]]
then
    usage
fi


TAG=""
eval set -- "$VALID_ARGS"
while true
do
    case "$1" in
        -r|--release)
            TAG="$(date -Idate)-alpine"
            shift
            ;;
        -b|--branch)
            TAG="$(date -Idate)-alpine-$(git rev-parse --short HEAD)"
            shift
            ;;
        --) shift;
            break
            ;;
    esac
done

if [[ -z "$TAG" ]]
then
    usage
fi

BASE_REPO="broad-dsde-methods/variantstore"
REPO_WITH_TAG="${BASE_REPO}:${TAG}"
GCR_TAG="us.gcr.io/${REPO_WITH_TAG}"

docker build . -t "${REPO_WITH_TAG}"

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
    echo "One or more unit test has failed, exiting."
    exit $fail
fi

set -o errexit

docker tag "${REPO_WITH_TAG}" "${GCR_TAG}"
#docker push "${GCR_TAG}"

echo "docker image pushed to \"${GCR_TAG}\""
