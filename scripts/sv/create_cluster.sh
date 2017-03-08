#!/bin/bash

# This script deletes a Google Dataproc cluster used for running the GATK-SV pipeline.

set -eu
if [[ "$#" -ne 3 ]]; then
    echo "Please provide the project, the name of cluster, and the bucket where the data lives"
    exit 1
fi

PROJECT=$1
CLUSTER_NAME=$2
INITIALIZATION_BUCKET=$3

gcloud dataproc clusters create ${CLUSTER_NAME} \
    --zone us-central1-a \
    --master-machine-type n1-highmem-8 \
    --num-worker-local-ssds 1 \
    --worker-machine-type n1-highmem-16 \
    --master-boot-disk-size 500 \
    --worker-boot-disk-size 500 \
    --num-workers 10 \
    --image-version preview \
    --project ${PROJECT} \
    --initialization-actions gs://${INITIALIZATION_BUCKET}/initialization_scripts/dataproc_initialization_new_bucket_copy_ref_img.sh \
    --initialization-action-timeout 60m

