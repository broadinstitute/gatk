set -o errexit -o nounset -o pipefail -o xtrace

usage() {
    echo "

USAGE: ./build_plink_docker.sh

Build a PLINK 2 (https://www.cog-genomics.org/plink/2.0/) Docker image with an appropriate tag and push to GAR.
Tags will be of the form <ISO 8601 Date>-slim-<Docker image ID>.

e.g. 2024-04-22-slim-f000ba44
"
    exit 1
}

if [[ $# -ne 0 ]]
then
    usage
fi

docker build . --iidfile idfile.txt

FULL_IMAGE_ID=$(cat idfile.txt)
IMAGE_ID=${FULL_IMAGE_ID:7:12}
IMAGE_TYPE="slim"
TAG=$(python3 ../build_docker_tag.py --image-id "${IMAGE_ID}" --image-type "${IMAGE_TYPE}")

BASE_REPO="broad-dsde-methods/gvs"
REPO_WITH_TAG="${BASE_REPO}/plink2:${TAG}"
docker tag "${IMAGE_ID}" "${REPO_WITH_TAG}"


GAR_TAG="us-central1-docker.pkg.dev/${REPO_WITH_TAG}"
docker tag "${REPO_WITH_TAG}" "${GAR_TAG}"

# Docker must be configured for GAR before pushes will work:
# gcloud auth configure-docker us-central1-docker.pkg.dev
docker push "${GAR_TAG}"

echo "Docker image pushed to \"${GAR_TAG}\""

rm idfile.txt
