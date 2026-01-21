# The Dockerfile in this directory is used to build a Docker image for Ensembl VEP 115 with GRCh38 LOFTEE support.
# Copy the Dockerfile over the `docker/Dockerfile` in the Ensembl VEP repo and copy this script to the root of the repo.
#
# This script is best run on an x86 VM. I was able to run it on my Apple Silicon Mac, but there were random crashes at
# several different points and I had to restart each time.
set -o errexit -o nounset -o pipefail -o xtrace

# Build and write the full image id to an id file.
docker buildx build --platform linux/amd64 --file docker/Dockerfile --load . --iidfile idfile.txt

# Scrape out the first 12 hex digits of the SHA256 hash to use for the tag.
IMAGE_ID=$(cut -d : -f 2 idfile.txt | cut -c 1-12)

rm idfile.txt

TAG="$(date -Idate)-${IMAGE_ID}"
BASE_REPO="broad-dsde-methods/gvs"
REPO_WITH_TAG="${BASE_REPO}/loftee:${TAG}"
docker tag "${IMAGE_ID}" "${REPO_WITH_TAG}"

# Tag and push
GAR_TAG="us-central1-docker.pkg.dev/${REPO_WITH_TAG}"
docker tag "${REPO_WITH_TAG}" "${GAR_TAG}"

docker push "${GAR_TAG}"

echo "Docker image pushed to \"${GAR_TAG}\""