if [ $# -lt 1 ]; then
    echo "USAGE: ./build_docker.sh [DOCKER_TAG_STRING]"
    echo " e.g.: ./build_docker.sh mybranch_2021_04_03"
    exit 1
fi

INFO=$1
GCR_TAG="us.gcr.io/broad-dsde-methods/variantstore:${INFO}"

docker build . -t broad-dsde-methods/variantstore:${INFO}
docker tag broad-dsde-methods/variantstore:${INFO} ${GCR_TAG}
docker push ${GCR_TAG}

echo "docker image pushed to \"${GCR_TAG}\""
