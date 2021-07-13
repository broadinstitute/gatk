if [ $# -lt 1 ]; then
    echo "USAGE: ./build_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g.: ./build_docker.sh mybranch_20210403"
    exit 1
fi

INFO=$1
BASE_REPO="broad-dsde-methods/variantstore"
REPO_WITH_TAG="${BASE_REPO}:${INFO}"
GCR_TAG="us.gcr.io/${REPO_WITH_TAG}"

docker build . -t ${REPO_WITH_TAG}
docker tag ${REPO_WITH_TAG} ${GCR_TAG}
docker push ${GCR_TAG}

echo "docker image pushed to \"${GCR_TAG}\""
