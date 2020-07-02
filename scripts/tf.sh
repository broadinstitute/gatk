#!/usr/bin/env bash

# Script to enable running Python modules within Docker containers
# Note: If 'nvidia-docker' is specified as $DOCKER_COMMAND, this script must be run on a dl-image machine
# (rather than a 'ukbb-image' machine) on a GPU-enabled machine.

################### VARIABLES ############################################

# The default images are based on ufoym/deepo:all-py36-jupyter
DOCKER_IMAGE_GPU="gcr.io/broad-ml4cvd/deeplearning:tf2-latest-gpu"
DOCKER_IMAGE_NO_GPU="gcr.io/broad-ml4cvd/deeplearning:tf2-latest-cpu"
DOCKER_IMAGE=${DOCKER_IMAGE_GPU}
GPU_DEVICE="--gpus all"
INTERACTIVE=""
MOUNTS=""
PYTHON_COMMAND="python"
TEST_COMMAND="python -m pytest"
SCRIPT_NAME=$( echo $0 | sed 's#.*/##g' )

################### USERNAME & GROUPS ####################################

# Get group names
GROUP_NAMES=$(groups ${USER} | sed -e 's/.*:\ //')

# Get group names as array to iterate through
GROUP_NAMES_ARR=( $GROUP_NAMES )

# Iterate through array, get group ID for each group name, append to string
GROUP_IDS=""
for GROUP_TO_ADD in "${GROUP_NAMES_ARR[@]}"; do
    GROUP_IDS="$GROUP_IDS $(getent group ${GROUP_TO_ADD} | grep -o -P '(?<=x:).*(?=:)')"
done

# Export environment variables so they can be passed into Docker and accessed in bash
export GROUP_NAMES GROUP_IDS

# Create string to be called in Docker's bash shell via eval;
# this creates a user, adds groups, adds user to groups, then calls the Python script
CALL_DOCKER_AS_USER="
    apt-get -y install sudo;
    useradd -u $(id -u) ${USER};
    GROUP_NAMES_ARR=( \${GROUP_NAMES} );
    GROUP_IDS_ARR=( \${GROUP_IDS} );
    for (( i=0; i<\${#GROUP_NAMES_ARR[@]}; ++i )); do
        echo \"Creating group\" \${GROUP_NAMES_ARR[i]} \"with gid\" \${GROUP_IDS_ARR[i]};
        groupadd -f -g \${GROUP_IDS_ARR[i]} \${GROUP_NAMES_ARR[i]};
        echo \"Adding user ${USER} to group\" \${GROUP_NAMES_ARR[i]}
        usermod -aG \${GROUP_NAMES_ARR[i]} ${USER}
    done;
    sudo -u ${USER}"

################### HELP TEXT ############################################

usage()
{
    cat <<USAGE_MESSAGE

    This script can be used to run a Python module within a Docker container.

    Usage: ${SCRIPT_NAME} [-nth] [-i <image>] module [arg ...]

    Example: ./${SCRIPT_NAME} -n -t -i gcr.io/broad-ml4cvd/deeplearning:latest-cpu recipes.py --mode tensorize ...

        -c                  if set use CPU docker image and machine and use the regular 'docker' launcher.
                            By default, we assume the machine is GPU-enabled.

        -d                  Select a particular GPU device on multi GPU machines

        -m                  Directories to mount at the same path in the docker image

        -t                  Run Docker container interactively.

        -j                  Set up Jupyter directory

        -r                  Call Python script as root. If this flag is not specified,
                            the owner and group of the output directory will be those
                            of the user who called the script.

        -h                  Print this help text.

        -i      <image>     Run Docker with the specified custom <image>. The default image is '${DOCKER_IMAGE}'.
        -T                  Run tests
USAGE_MESSAGE
}

################### OPTION PARSING #######################################

while getopts ":i:d:m:ctjrhT" opt ; do
    case ${opt} in
        h)
            usage
            exit 1
            ;;
        i)
            DOCKER_IMAGE=$OPTARG
            ;;
        d)
            GPU_DEVICE="--gpus device=${OPTARG}"
            ;;
        m)
            MOUNTS="-v ${OPTARG}:${OPTARG}"
            ;;
        c)
            DOCKER_IMAGE=${DOCKER_IMAGE_NO_GPU}
            GPU_DEVICE=
            ;;
        t)
            INTERACTIVE="-it"
            ;;
        j)  # Set up Jupyter
            mkdir -p /home/${USER}/jupyter/
            chmod o+w /home/${USER}/jupyter/
            mkdir -p /home/${USER}/jupyter/root/
            mkdir -p /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/
            ;;
        r) # Output owned by root
            CALL_DOCKER_AS_USER=""
            ;;
        T)
            PYTHON_COMMAND=${TEST_COMMAND}
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
    echo "Could not pull the image ${DOCKER_IMAGE}. Will try anyway..."
fi

if [[ -d "/data" ]] ; then
    echo "Found /data folder will try to mount it."
    MOUNTS="${MOUNTS} -v /data/:/data/"
fi

if [[ -d "/mnt" ]] ; then
    echo "Found /mnt folder will try to mount it."
    MOUNTS="${MOUNTS} -v /mnt/:/mnt/"
fi

# Get your external IP directly from a DNS provider
WANIP=$(dig +short myip.opendns.com @resolver1.opendns.com)

# Let anyone run this script
WORKDIR=$(pwd)

PYTHON_ARGS="$@"
cat <<LAUNCH_MESSAGE
Attempting to run Docker with
    docker run ${INTERACTIVE} \
    ${GPU_DEVICE} \
    --env GROUP_NAMES \
    --env GROUP_IDS \
    --rm \
    --ipc=host \
    -v ${WORKDIR}/:${WORKDIR}/ \
    -v ${HOME}/:${HOME}/ \
    ${MOUNTS} \
    ${DOCKER_IMAGE} /bin/bash -c "pip install ${WORKDIR};
        eval ${CALL_DOCKER_USER} ${PYTHON_COMMAND} ${PYTHON_ARGS}"
LAUNCH_MESSAGE

docker run ${INTERACTIVE} \
${GPU_DEVICE} \
--env GROUP_NAMES \
--env GROUP_IDS \
--rm \
--ipc=host \
-v ${WORKDIR}/:${WORKDIR}/ \
-v ${HOME}/:${HOME}/ \
${MOUNTS} \
${DOCKER_IMAGE} /bin/bash -c "pip install ${WORKDIR}; eval ${CALL_DOCKER_AS_USER} ${PYTHON_COMMAND} ${PYTHON_ARGS}"
