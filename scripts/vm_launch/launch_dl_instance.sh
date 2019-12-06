#!/usr/bin/env bash

NAME=${1:-sam-p4}
shift 1
INSTANCE_TYPE=${1:-n1-standard-4}
shift 1
DISK_SIZE=${1:-100GB}
shift 1
ACCEL=${1:-nvidia-tesla-k80}
shift 1

echo "Creating GPU instance ${NAME} from family dl-image of type ${INSTANCE_TYPE} with GPU ${ACCEL}..."

echo "$@"

gcloud compute instances create ${NAME} \
--project broad-ml4cvd \
--zone us-central1-a \
--image-project broad-ml4cvd \
--image-family dl-image \
--accelerator=type=${ACCEL},count=1 \
--maintenance-policy=TERMINATE \
--boot-disk-type=pd-ssd \
--boot-disk-size=${DISK_SIZE} \
--service-account 783282864357-compute@developer.gserviceaccount.com \
--scopes https://www.googleapis.com/auth/cloud-platform \
--machine-type ${INSTANCE_TYPE} \
--metadata startup-script-url=gs://ml4cvd/projects/jamesp/home/startup.sh \
"$@"

# Previously used the base ubuntu:
# --image-project ubuntu-os-cloud \
# --image-family ubuntu-1804-lts \

# You can choose whatever size you like for the boot disk:
# --boot-disk-size 300GB
