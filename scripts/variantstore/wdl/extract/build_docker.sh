if [ $# -lt 1 ]; then
    echo "USAGE: ./build_docker.sh [DOCKER_TAG_STRING] [OPTIONAL:LATEST]"
    echo " e.g.: ./build_docker.sh mybranch_20210403"
    exit 1
fi

# Test that the VAT python code has not been broken
python -m unittest create_variant_annotation_table_test.py  > /dev/null
VAT_TEST_RESULTS=$?
if [ $VAT_TEST_RESULTS -ne 0 ]; then
    echo "TestMakeAnnotatedJsonRow python test has failed"
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
