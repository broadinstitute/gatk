#!/usr/bin/env bash

if [ $# -lt 1 ]; then
    echo "USAGE: ./build_build_base_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g.: ./build_build_base_docker.sh \$(date -Idate)-alpine-build-base"
    exit 1
fi

if [[ ! "$1" == *-build-base ]]; then
    echo "Specified tag '$1' does end with '-build-base'."
    echo "build_build_base_docker.sh is intended for building build base images only."
    exit 1
fi

set -o errexit -o nounset -o pipefail -o xtrace

BASE_REPO="broad-dsde-methods/variantstore"
TAG="${1}"
REPO_WITH_TAG="${BASE_REPO}:${TAG}"
GCR_TAG="us.gcr.io/${REPO_WITH_TAG}"

set +o errexit
grep "${TAG}" <(docker images | grep "^${BASE_REPO}" | awk '{print $2}') > /dev/null
RC="$?"
set -o errexit

if [[ $RC -eq 0 ]]
then
  echo "Error: Tag '${TAG}' already exists for image '${BASE_REPO}', refusing to tag again."
  exit 1
fi

docker build --file build_base.Dockerfile . -t "${REPO_WITH_TAG}"

# Smoke test: make sure the gsutil command does not crash with a no-args invocation.
set +o errexit
docker run --rm -it "${REPO_WITH_TAG}" gsutil > /dev/null
RC="$?"
set -o errexit

if [[ $RC -ne 0 ]]
then
  echo "Error: image '${REPO_WITH_TAG}' did not pass smoke test, will need to be deleted before attempting to tag again."
  exit 1
fi

docker tag "${REPO_WITH_TAG}" "${GCR_TAG}"
docker push "${GCR_TAG}"

echo "docker image pushed to \"${GCR_TAG}\""
