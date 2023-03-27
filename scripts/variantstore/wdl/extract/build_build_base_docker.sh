#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "USAGE: ./build_build_base_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g.: ./build_build_base_docker.sh 2022-11-04-alpine-build-base"
    exit 1
fi

if [[ ! "$1" == *-build-base ]]; then
    echo "Specified tag '$1' does end with '-build-base'."
    echo "build_build_base_docker.sh is intended for building build base images only."
    exit 1
fi

set -o errexit -o nounset -o pipefail -o xtrace

BASE_REPO="broad-dsde-methods/variantstore"
REPO_WITH_TAG="${BASE_REPO}:${1}"
GCR_TAG="us.gcr.io/${REPO_WITH_TAG}"

docker build --file build_base.Dockerfile . -t "${REPO_WITH_TAG}"

docker tag "${REPO_WITH_TAG}" "${GCR_TAG}"
docker push "${GCR_TAG}"

echo "docker image pushed to \"${GCR_TAG}\""
