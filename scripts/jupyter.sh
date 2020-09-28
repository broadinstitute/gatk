#!/usr/bin/env bash

# Script to enable running Python modules within Docker containers
# Note: If 'nvidia-docker' is specified as $DOCKER_COMMAND, this script must be run on a dl-image machine
# (rather than a 'ukbb-image' machine) on a GPU-enabled machine.

################### VARIABLES ############################################

# The default images are based on ufoym/deepo:all-py36-jupyter
DOCKER_IMAGE_GPU="gcr.io/broad-ml4cvd/deeplearning:tf2-latest-gpu"
DOCKER_IMAGE_NO_GPU="gcr.io/broad-ml4cvd/deeplearning:tf2-latest-cpu"
DOCKER_IMAGE=${DOCKER_IMAGE_GPU}
DOCKER_COMMAND="docker"
PORT="8888"
SCRIPT_NAME=$( echo $0 | sed 's#.*/##g' )
GPU_DEVICE="--gpus all"

################### HELP TEXT ############################################

usage()
{
    cat <<USAGE_MESSAGE

    This script can be used to run a Jupyter server ona VM with a tunnel to specific port.

    Usage: ${SCRIPT_NAME} [-nth] [-i <image>] module [arg ...]

    Example: ./${SCRIPT_NAME} -n -p 8889  -i gcr.io/broad-ml4cvd/deeplearning:latest-cpu

        -c                  Use CPU docker image and use the regular 'docker' launcher.
                            By default, 'nvidia-docker' wrapper is used to launch Docker assuming the machine is GPU-enabled.

        -h                  Print this help text.

        -p                  Port to use, by default '${PORT}'

        -i      <image>     Run Docker with the specified custom <image>. The default image is '${DOCKER_IMAGE}'.
USAGE_MESSAGE
}

################### OPTION PARSING #######################################

while getopts ":ip:ch" opt ; do
    case ${opt} in
        h)
            usage
            exit 1
            ;;
        i)
            DOCKER_IMAGE=$OPTARG
            ;;
        p)
            PORT=$OPTARG
            ;;
        c)
            DOCKER_IMAGE=${DOCKER_IMAGE_NO_GPU}
            GPU_DEVICE=""
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


################### SCRIPT BODY ##########################################

if ! docker pull ${DOCKER_IMAGE}; then
   echo "ERROR: Could not pull the image ${DOCKER_IMAGE}. Aborting..."
   exit 1;
fi

# Get your external IP directly from a DNS provider
WANIP=$(dig +short myip.opendns.com @resolver1.opendns.com)

# Let anyone run this script
USER=$(whoami)

echo 'Launching docker.'
echo 'Locally, you should now run:'
echo $(tput setaf 1)$(tput setab 7)"gcloud compute ssh ${USER}@$(hostname) -- -NnT -L ${PORT}:localhost:${PORT}"$(tput sgr 0)
echo 'or:'
echo $(tput setaf 1)$(tput setab 7)"ssh -nNT -L ${PORT}:localhost:${PORT} ${WANIP}"$(tput sgr 0)

mkdir -p /home/${USER}/jupyter/
chmod o+w /home/${USER}/jupyter/

mkdir -p /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/

${DOCKER_COMMAND} run -it \
${GPU_DEVICE} \
--rm \
--ipc=host \
--hostname=$(hostname) \
-v /home/${USER}/:/home/${USER}/ \
-v /mnt/:/mnt/ \
-p 0.0.0.0:${PORT}:${PORT} \
${DOCKER_IMAGE} /bin/bash -c "pip install -e /home/${USER}/ml; jupyter notebook --no-browser --ip=0.0.0.0 --port=${PORT} --NotebookApp.token= --allow-root --notebook-dir=/home/${USER}"


# Automatically back up any local notebooks and artifacts non-recursively (no subfolders)
echo 'Backing up local files to' /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/
cp /home/${USER}/jupyter/* /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/
