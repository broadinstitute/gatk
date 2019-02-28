#!/usr/bin/env bash

# Have script stop if there is an error
set -e

# Gatk Base features
BASE_REPO=jamesemery
BASE_PROJECT=gatknightly
BASE_VERSION=2.0.3

# Gatk Python and R image features
INTERMEDIATE_REPO=jamesemery
INTERMEDIATE_PROJECT=gatknightly
INTERMEDIATE_VERSION=0.0.1

BASE_FULL_PATH=${BASE_REPO}/${BASE_PROJECT}:gatkbase-${BASE_VERSION}
INTERMEDIATE_FULL_PATH=${INTERMEDIATE_REPO}/${INTERMEDIATE_PROJECT}:gatklibs-${INTERMEDIATE_VERSION}
#################################################
# Parsing arguments
#################################################
while getopts "ph" option; do
	case "$option" in
	    h) IS_HELP=true ;;
		p) IS_PUSH=true ;;
	esac
done

if [ -n "$IS_HELP" ]; then
	printf "Build the GATK4 base image.\n \
Usage: %s: [-p] [-h] \n \
Optional arguments:  \n \
-p \t push image to docker hub once complete. \n \n" $0
	exit 1
fi

# Output the parameters
echo -e "\n"
echo -e "docker hub repo, project, and tag: ${FULL_PATH}\n\n"
echo "Other options (Blank is false)"
echo "---------------"
echo "Push to dockerhub: ${IS_PUSH}"

# Login to dockerhub
if [ -n "${IS_PUSH}" ]; then
	echo "Please login to dockerhub"
	docker login
fi

# Build
echo "Building image to tag ${BASE_FULL_PATH}..."
docker build -f scripts/docker/gatkbase/Dockerfile --squash -t ${BASE_FULL_PATH} scripts/docker/gatkbase/

echo "Building image to tag ${INTERMEDIATE_FULL_PATH}..."
docker build -f scripts/docker/gatkPythonR/Dockerfile -t ${INTERMEDIATE_FULL_PATH} --build-arg ZIPPATH=./unzippedJar --build-arg BASEIMAGE=${BASE_FULL_PATH} .

## Push
if [ -n "${IS_PUSH}" ]; then
	docker push ${BASE_FULL_PATH}
	docker push ${INTERMEDIATE_FULL_PATH}
else
	echo "Not pushing to dockerhub"
fi

