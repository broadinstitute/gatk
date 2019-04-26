#!/usr/bin/env bash

# Script to enable running Python modules within Docker containers
# Note: If 'nvidia-docker' is specified as $DOCKER_COMMAND, this script must be run on a dl-image machine
# (rather than a 'ukbb-image' machine) on a GPU-enabled machine.

################### VARIABLES ############################################

# The default images are based on ufoym/deepo:all-py36-jupyter
DOCKER_IMAGE_GPU="gcr.io/broad-ml4cvd/deeplearning:latest-gpu"
DOCKER_IMAGE_NO_GPU="gcr.io/broad-ml4cvd/deeplearning:latest-cpu"
DOCKER_IMAGE=${DOCKER_IMAGE_GPU}
DOCKER_COMMAND="nvidia-docker"
INTERACTIVE=""
SCRIPT_NAME=$( echo $0 | sed 's#.*/##g' )

################### HELP TEXT ############################################

usage()
{
    cat <<USAGE_MESSAGE

    This script can be used to run a Python module within a Docker container.

    Usage: ${SCRIPT_NAME} [-nth] [-i <image>] module [arg ...]

    Example: ./${SCRIPT_NAME} -n -t -i gcr.io/broad-ml4cvd/deeplearning:latest-cpu recipes.py --mode tensorize ...

        -n                  Assume non-GPU-enabled machine and use the regular 'docker' launcher.
                            By default, 'nvidia-docker' wrapper is used to launch Docker assuming the machine is GPU-enabled.

        -t                  Run Docker container interactively.

        -h                  Print this help text.

        -i      <image>     Run Docker with the specified custom <image>. The default image is '${DOCKER_IMAGE}'.
USAGE_MESSAGE
}

################### OPTION PARSING #######################################

while getopts ":i:nth" opt ; do
    case ${opt} in
        h)
            usage
            exit 1
            ;;
        i)
            DOCKER_IMAGE=$OPTARG
            ;;
        n)
            DOCKER_IMAGE=${DOCKER_IMAGE_NO_GPU}
            DOCKER_COMMAND=docker
            ;;
        t)
            INTERACTIVE_RUN="-it"
            ;;
        :)
            echo "ERROR: Option -${OPTARG} requires an argument." 1>&2
            usage
            exit 1
            ;;
        *)
            echo "ERROR: Invalid option: -${OPTARG}" 1>&2
            usage
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

if [[ $# -eq 0 ]]; then
    echo "ERROR: No Python module was specified." 1>&2
    usage
    exit 1
fi

################### SCRIPT BODY ##########################################

if ! docker pull ${DOCKER_IMAGE}; then
    echo "ERROR: Could not pull the image ${DOCKER_IMAGE}. Aborting..."
    exit 1;
fi

# Get your external IP directly from a DNS provider
WANIP=$(dig +short myip.opendns.com @resolver1.opendns.com)

# Let anyone run this script
USER=$(whoami)

mkdir -p /home/${USER}/jupyter/
chmod o+w /home/${USER}/jupyter/
mkdir -p /home/${USER}/jupyter/root/
mkdir -p /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/

cat <<LAUNCH_MESSAGE
Attempting to run Docker with
    ${DOCKER_COMMAND} run ${INTERACTIVE}
        --rm
        --ipc=host
        -v /home/${USER}/jupyter/root/:/root/
        -v /home/${USER}/:/home/${USER}/
        -v /mnt/:/mnt/
        ${DOCKER_IMAGE} python $@
LAUNCH_MESSAGE

${DOCKER_COMMAND} run ${INTERACTIVE} \
--rm \
--ipc=host \
-v /home/${USER}/jupyter/root/:/root/ \
-v /home/${USER}/:/home/${USER}/ \
-v /mnt/:/mnt/ \
${DOCKER_IMAGE} python "$@"
