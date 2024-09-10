#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

if [ $# -lt 1 ]; then
    echo "USAGE: ./build_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g.: ./build_docker.sh myrepo/nirvana:v3.18.1"
    exit 1
fi

DOCKER_TAG="$1"

docker build . -t "${DOCKER_TAG}"
