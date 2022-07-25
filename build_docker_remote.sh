#!/usr/bin/env bash
#
# Build (and optionally push) a GATK docker image to GCR using Google Cloud Build. Images are built in the cloud rather than locally. Pushing to dockerhub is not supported by this script.
#
# If you are pushing an image to our release repositories, be sure that you've followed
# the setup instructions here:
# https://github.com/broadinstitute/gatk/wiki/How-to-release-GATK4#setup_docker
# and here:
# https://github.com/broadinstitute/gatk/wiki/How-to-release-GATK4#setup_gcloud
#

# Have script stop if there is an error
set -e

REPO=broadinstitute
PROJECT=gatk
GCR_REPO="us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds"
STAGING_CLONE_DIR=${PROJECT}_staging_temp

#################################################
# Parsing arguments
#################################################
while getopts "e:sd:t:" option; do
	case "$option" in
		e) GITHUB_TAG="$OPTARG" ;;
		s) IS_HASH=true ;;
		d) STAGING_DIR="$OPTARG" ;;
		t) DOCKER_IMAGE_TAG="$OPTARG" ;;
	esac
done

if [ -z "$GITHUB_TAG" ]; then
	printf "Option -e requires an argument.\n \
Usage: %s: -e <GITHUB_TAG> [-sdt] \n \
where <GITHUB_TAG> is the github tag (or hash when -s is used) to use in building the docker image\n \
(e.g. bash build_docker_remote.sh -e 4.2.6.1 )\n \
Optional arguments:  \n \
-s \t The GITHUB_TAG (-e parameter) is actually a github hash, not tag.  \n \
-d <STAGING_DIR> \t staging directory to grab code from repo and build the docker image.  If unspecified, then use whatever is in current dir (do not go to the repo).  NEVER SPECIFY YOUR WORKING DIR \n \
-t <IMAGE_TAG>\t  The tag to assign image once it is finished constructing.  NOTE: currently this MUST be on either GCR or the Google Artifact Registry. \n" $0
	exit 1
fi

# Output the parameters
echo -e "\n"
echo -e "github tag/hash: ${GITHUB_TAG}"
echo -e "github project: ${REPO}/${PROJECT}:${GITHUB_TAG}\n\n"
echo "Other options (Blank is false)"
echo "---------------"
echo "This is a git hash: ${IS_HASH}"
echo "Staging directory: ${STAGING_DIR}"

ORIGINAL_WORKING_DIRECTORY=$(pwd)

if [ -n "$STAGING_DIR" ]; then
    GITHUB_DIR="tags/"

    if [ -n "${IS_HASH}" ]; then
        GITHUB_DIR=" "
    fi

    mkdir -p ${STAGING_DIR}
    cd ${STAGING_DIR}
    set +e
    rm -Rf ${STAGING_CLONE_DIR}
    set -e
    GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/${REPO}/${PROJECT}.git ${STAGING_CLONE_DIR}
    cd ${STAGING_CLONE_DIR}
    STAGING_ABSOLUTE_PATH=$(pwd)

    echo "Now in $(pwd)"
    GIT_CHECKOUT_COMMAND="git checkout ${GITHUB_DIR}${GITHUB_TAG}"
    echo "${GIT_CHECKOUT_COMMAND}"
    ${GIT_CHECKOUT_COMMAND}
fi

GIT_HASH_FOR_TAG=$(git describe --tags --long)

## generate the tag if it wasn't explicitly specified
if [ -z "$DOCKER_IMAGE_TAG" ]; then
  DOCKER_IMAGE_TAG=${GCR_REPO}:$(whoami)-${GITHUB_TAG}-${GIT_HASH_FOR_TAG}
fi

echo "Building image with the tag ${DOCKER_IMAGE_TAG}..."

SUBMIT_COMMAND="gcloud builds submit --tag ${DOCKER_IMAGE_TAG} --timeout=24h --machine-type n1_highcpu_8"
echo "running the following gcloud command: ${SUBMIT_COMMAND}"
## We need to override the default .gcloudignore to preserve the .git directory which we need in order to download LFS files in the remote build.
echo -n "" >> .gcloudignore
${SUBMIT_COMMAND}

cd ${ORIGINAL_WORKING_DIRECTORY}
if [ -n "$STAGING_DIR" ] ; then
    rm -Rf ${STAGING_DIR}/${STAGING_CLONE_DIR}
fi

exit 0