#!/usr/bin/env bash

# Note: This must be run on a dl-image machine, rather than a ukbb-image
# machine. It expects to be run on a GPU-enabled machine.
PORT=${1:-8888}
shift 1

DOCKER_IMAGE=gcr.io/broad-ml4cvd/deeplearning:latest
docker pull ${DOCKER_IMAGE}
#Based on ufoym/deepo:all-py36-jupyter

# Get your external IP directly from a DNS provider
WANIP=$(dig +short myip.opendns.com @resolver1.opendns.com)

# Let anyone run this script
USER=$(whoami)

echo 'Launching docker.'
echo 'Locally, you should now run:'
echo $(tput setaf 1)$(tput setab 7)"gcloud compute ssh ${USER}@$(hostname) -- -NnT -L ${PORT}:localhost:8888"$(tput sgr 0)
echo 'or:'
echo $(tput setaf 1)$(tput setab 7)"ssh -i ~/.ssh/google_compute_engine -nNT -L ${PORT}:localhost:8888 ${WANIP}"$(tput sgr 0)

mkdir -p /home/${USER}/jupyter/
chmod o+w /home/${USER}/jupyter/

mkdir -p /home/${USER}/jupyter/root/

mkdir -p /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/

nvidia-docker run -it \
--rm \
--ipc=host \
-v /home/${USER}/:/home/${USER}/ \
-v /mnt/:/mnt/ \
-p ${PORT}:${PORT} \
${DOCKER_IMAGE} jupyter notebook --no-browser --ip=0.0.0.0 --NotebookApp.token= --allow-root --notebook-dir=/home/${USER}

# Automatically back up any local notebooks and artifacts non-recursively (no subfolders)
echo 'Backing up local files to' /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/
cp /home/${USER}/jupyter/* /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/
