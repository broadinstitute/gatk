#!/usr/bin/env bash

# Have script stop if there is an error
set -e


REPO=broadinstitute
PROJECT=gatk
VERSION=1.2
FULL_PATH=${REPO}/${PROJECT}:gatkbase-${VERSION}

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

# Make sure sudo or root was used.
if [ "$(whoami)" != "root" ]; then
	echo "You must have superuser privileges (through sudo or root user) to run this script"
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
echo "Building image to tag ${FULL_PATH}..."
docker build -t ${FULL_PATH} .

## Push
if [ -n "${IS_PUSH}" ]; then
	docker push ${FULL_PATH}
else
	echo "Not pushing to dockerhub"
fi

