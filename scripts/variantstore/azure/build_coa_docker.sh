if [ $# -lt 1 ]; then
    echo "USAGE: ./build_coa_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g.: ./build_coa_docker.sh $(date -I)"
    exit 1
fi

set -o errexit -o nounset -o pipefail -o xtrace

BASE_REPO="broad-dsde-methods/variantstore"
REPO_WITH_TAG="${BASE_REPO}:coa-${1}"
GCR_TAG="us.gcr.io/${REPO_WITH_TAG}"

docker build . -t "${REPO_WITH_TAG}" -f cromwell_on_azure.Dockerfile

docker tag "${REPO_WITH_TAG}" "${GCR_TAG}"
docker push "${GCR_TAG}"

echo "Docker image pushed to \"${GCR_TAG}\""
