INFO=$1
GCR_TAG="us.gcr.io/broad-dsde-methods/variantstore:${INFO}"

docker build . -t broad-dsde-methods/variantstore:${INFO}
docker tag broad-dsde-methods/variantstore:${INFO} ${GCR_TAG}
docker push ${GCR_TAG}

echo "docker pushed to ${GCR_TAG}"
