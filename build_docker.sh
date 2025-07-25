#!/usr/bin/env bash
#
# Build (and optionally push) a GATK docker image.
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
REPO_PRJ=${REPO}/${PROJECT}
GCR_REPO="us.gcr.io/broad-gatk/gatk"
STAGING_CLONE_DIR=${PROJECT}_staging_temp

#################################################
# Parsing arguments
#################################################
while getopts "e:pslrud:t:m" option; do
	case "$option" in
		e) GITHUB_TAG="$OPTARG" ;;
		p) IS_PUSH=true ;;
		s) IS_HASH=true ;;
		l) IS_NOT_LATEST=true ;;
		r) IS_NOT_REMOVE_UNIT_TEST_CONTAINER=true ;;
		u) IS_NOT_RUN_UNIT_TESTS=true ;;
		d) STAGING_DIR="$OPTARG" ;;
		t) PULL_REQUEST_NUMBER="$OPTARG" ;;
        m) LITE=true ;;
	esac
done

if [ -z "$GITHUB_TAG" ]; then
	printf "Option -e requires an argument.\n \
Usage: %s: -e <GITHUB_TAG> [-psl] \n \
where <GITHUB_TAG> is the github tag (or hash when -s is used) to use in building the docker image\n \
(e.g. bash build_docker.sh -e 1.0.0.0-alpha1.2.1)\n \
Optional arguments:  \n \
-s \t The GITHUB_TAG (-e parameter) is actually a github hash, not tag.  git hashes cannot be pushed as latest, so -l is implied.  \n \
-l \t Do not also push the image to the 'latest' tag. \n \
-u \t Do not run the unit tests. \n \
-m \t Build the lite image (which does not contain the conda environment). \n \
-d <STAGING_DIR> \t staging directory to grab code from repo and build the docker image.  If unspecified, then use whatever is in current dir (do not go to the repo).  NEVER SPECIFY YOUR WORKING DIR \n \
-p \t (GATK4 developers only) push image to docker hub once complete.  This will use the GITHUB_TAG in dockerhub as well. \n \
\t\t Unless -l is specified, this will also push this image to the 'latest' tag. \n \
-r \t (GATK4 developers only) Do not remove the unit test docker container.  This is useful for debugging failing unit tests. \n \
-t <PULL_REQUEST_NUMBER>\t (Travis CI only) The pull request number.  This is only used during pull request builds on Travis CI. \n" $0
	exit 1
fi

# -z is like "not -n"
if [ -z ${IS_NOT_LATEST} ] && [ -n "${IS_HASH}" ] && [ -n "${IS_PUSH}" ]; then
	echo -e "\n##################"
	echo " WARNING:  Refusing to push a hash as latest to dockerhub. "
	echo "##################"
	IS_NOT_LATEST=true
fi


# Output the parameters
echo -e "\n"
echo -e "github tag/hash: ${GITHUB_TAG}"
echo -e "docker hub repo, project, and tag: ${REPO_PRJ}:${GITHUB_TAG}\n\n"
echo "Other options (Blank is false)"
echo "---------------"
echo "This is a git hash: ${IS_HASH}"
echo "Push to dockerhub: ${IS_PUSH}"
echo "Do NOT remove the unit test container: ${IS_NOT_REMOVE_UNIT_TEST_CONTAINER}"
echo "Staging directory: ${STAGING_DIR}"
echo -e "Do NOT make this the default docker image: ${IS_NOT_LATEST}"
echo -e "Fetch from this remote path: ${PULL_REQUEST_NUMBER}\n\n\n"

# Login to dockerhub
if [ -n "${IS_PUSH}" ]; then
	echo "Please login to dockerhub"
	docker login
fi

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
    if [ ${PULL_REQUEST_NUMBER} ]; then
        GIT_FETCH_COMMAND="git fetch origin +refs/pull/${PULL_REQUEST_NUMBER}/merge"
        echo "${GIT_FETCH_COMMAND}"
        ${GIT_FETCH_COMMAND}
    fi
    GIT_CHECKOUT_COMMAND="git checkout ${GITHUB_DIR}${GITHUB_TAG}"
    echo "${GIT_CHECKOUT_COMMAND}"
    ${GIT_CHECKOUT_COMMAND}
    # Since the large runtime resources are compiled into the jar, they have to be available to
    # the build when it's done as part of the Docker build. Although the build itself will pull
    # them from the lfs server if they're not present, the Docker doesn't have git-lfs installed
    # so that would fail. So pull them into the staging areas so they'll be copied directly.
    GIT_PULL_LARGE_COMMAND="git lfs pull --include src/main/resources/large/"
    echo ${GIT_PULL_LARGE_COMMAND}
    ${GIT_PULL_LARGE_COMMAND}
fi

# Build
if [ -n "${IS_PUSH}" ]; then
    RELEASE=true
else
    RELEASE=false
fi

echo "Building image to tag ${REPO_PRJ}:${GITHUB_TAG}..."
if [ -n "${IS_PUSH}" ]; then
    if [ -n "${LITE}" ]; then
        docker build -t ${REPO_PRJ}:${GITHUB_TAG} --build-arg RELEASE=${RELEASE} --target gatk-lite --squash .
    else
        docker build -t ${REPO_PRJ}:${GITHUB_TAG} --build-arg RELEASE=${RELEASE}  --squash .
    fi
else
    if [ -n "${LITE}" ]; then
        docker build -t ${REPO_PRJ}:${GITHUB_TAG} --build-arg RELEASE=${RELEASE} --target gatk-lite .
    else
        docker build -t ${REPO_PRJ}:${GITHUB_TAG} --build-arg RELEASE=${RELEASE} .
    fi
fi
# Since we build the docker image with stages, the first build stage for GATK will be leftover in
# the local docker context after executing the above commands. This step reclaims that space automatically
# by removing the intermediate image based on a tag we assign it.
docker image prune -f --filter label=stage=gatkIntermediateBuildImage

if [ -z "${IS_NOT_RUN_UNIT_TESTS}" ] ; then

	# Run unit tests
	echo "Running unit tests..."
	REMOVE_CONTAINER_STRING=" --rm "
	if [ -n "${IS_NOT_REMOVE_UNIT_TEST_CONTAINER}" ] ; then
		REMOVE_CONTAINER_STRING=" "
	fi

	git lfs pull
    chmod -R a+w ${STAGING_ABSOLUTE_PATH}/src/test/resources

    cp build.gradle build.gradle.backup
    cp scripts/docker/dockertest.gradle .

	echo docker run ${REMOVE_CONTAINER_STRING} -v  ${STAGING_ABSOLUTE_PATH}:/gatkCloneMountPoint -v  ${STAGING_ABSOLUTE_PATH}/testJars:/jars -t ${REPO_PRJ}:${GITHUB_TAG} bash /root/run_unit_tests.sh
    docker run ${REMOVE_CONTAINER_STRING} -v  ${STAGING_ABSOLUTE_PATH}:/gatkCloneMountPoint -v  ${STAGING_ABSOLUTE_PATH}/testJars:/jars -t ${REPO_PRJ}:${GITHUB_TAG} bash /root/run_unit_tests.sh
	echo " Unit tests passed..."
	mv build.gradle.backup build.gradle
fi

## Push
if [ -n "${IS_PUSH}" ]; then

	echo "Pushing to ${REPO_PRJ}"
	docker push ${REPO_PRJ}:${GITHUB_TAG}

	echo "Pushing to ${GCR_REPO} (if this fails, run \"gcloud auth configure-docker\" and try again)"
	docker tag ${REPO_PRJ}:${GITHUB_TAG} ${GCR_REPO}:${GITHUB_TAG}
	docker push ${GCR_REPO}:${GITHUB_TAG}

	if [ -z "${IS_NOT_LATEST}" ] && [ -z "${IS_HASH}" ] ; then
		echo "Updating latest tag in ${REPO_PRJ}"
		docker tag ${REPO_PRJ}:${GITHUB_TAG} ${REPO_PRJ}:latest
		docker push ${REPO_PRJ}:latest
		
		echo "Updating latest tag in ${GCR_REPO} (if this fails, run \"gcloud auth configure-docker\" and try again)"
		docker tag ${GCR_REPO}:${GITHUB_TAG} ${GCR_REPO}:latest
		docker push ${GCR_REPO}:latest
	fi

else
	echo "Not pushing to dockerhub"
fi

cd ${ORIGINAL_WORKING_DIRECTORY}
if [ -n "$STAGING_DIR" ] ; then
    rm -Rf ${STAGING_DIR}/${STAGING_CLONE_DIR}
fi
