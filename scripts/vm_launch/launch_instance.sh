#!/usr/bin/env bash

NAME=${1:-jpp-1}
shift 1
INSTANCE_TYPE=${1:-n1-standard-1}
shift 1
DISK_SIZE=${1:-100GB}
shift 1

echo "Creating instance ${NAME} from family ml4h-image of type ${INSTANCE_TYPE}..."

echo "$@"

gcloud compute instances create ${NAME} \
--project broad-ml4cvd \
--zone us-central1-a \
--image-project broad-ml4cvd \
--image-family ml4cvd-image \
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
