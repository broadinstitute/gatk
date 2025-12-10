The Dockerfile in this directory is used to build a Docker image for Ensembl VEP 115 with GRCh38 LOFTEE support.
This file is intended to replace the Dockerfile at `docker/Dockerfile` in the Ensembl VEP repo. On an x86 VM
the image is can be built from the root of the Ensembl VEP repo with the command:

```
docker build -f docker/Dockerfile .
```

Get the image id from `docker images` and assign:

```
IMAGE_ID=<image id>
```

Then:

```
TAG="$(date -Idate)-${IMAGE_ID}"
BASE_REPO="broad-dsde-methods/gvs"
REPO_WITH_TAG="${BASE_REPO}/loftee:${TAG}"
docker tag "${IMAGE_ID}" "${REPO_WITH_TAG}"

# Tag and push
GAR_TAG="us-central1-docker.pkg.dev/${REPO_WITH_TAG}"
docker tag "${REPO_WITH_TAG}" "${GAR_TAG}"

docker push "${GAR_TAG}"

echo "Docker image pushed to \"${GAR_TAG}\""
```
