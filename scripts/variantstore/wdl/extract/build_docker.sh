if [ $# -lt 1 ]; then
    echo "USAGE: ./build_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g. to build and tag only with your specific tag: ./build_docker.sh mybranch_2021_04_03"
    echo " e.g. to build and tag with your specific tag and latest tag: ./build_docker.sh mybranch_2021_04_03 latest"
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


if [ $# -eq 2 ] && [ $2 == "latest" ]; then
  echo "additionally pushing this version to 'latest'"
  GCR_BASE="us.gcr.io/${BASE_REPO}"
  docker tag ${REPO_WITH_TAG} ${GCR_BASE}
  docker push ${GCR_BASE}

  echo "docker image pushed to \"${GCR_BASE}\""
fi


