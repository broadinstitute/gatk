#!/usr/bin/env bash

docker pull jupyter/datascience-notebook
  
# Get your external IP directly from a DNS provider
WANIP=$(dig +short myip.opendns.com @resolver1.opendns.com)

# Let anyone run this script
USER=$(whoami)

echo 'Launching docker.'
echo 'Locally, you should now run:'
echo $(tput setaf 1)$(tput setab 7)"gcloud compute ssh ${USER}@$(hostname) -- -NnT -L 8888:localhost:8888"$(tput sgr 0)
echo 'or:'
echo $(tput setaf 1)$(tput setab 7)"ssh -nNT -L 8888:localhost:8888 ${WANIP}"$(tput sgr 0)

mkdir -p /home/${USER}/jupyter/
chmod o+w /home/${USER}/jupyter/

mkdir -p /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/

docker run \
--rm \
-v /mnt/:/mnt/ \
-v /home/${USER}/:/home/${USER}/ \
-p 127.0.0.1:8888:8888 \
jupyter/datascience-notebook

# Automatically back up any local notebooks and artifacts non-recursively (no subfolders)
echo 'Backing up local files to' /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/
cp /home/${USER}/jupyter/* /mnt/ml4cvd/projects/${USER}/projects/jupyter/auto/
