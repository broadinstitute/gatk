#!/bin/bash

# This script deletes a Google Dataproc cluster used for running the GATK-SV pipeline.

set -eu
if [[ "$#" -ne 3 ]]; then
    echo "Please provide the project, the name of cluster, and the bucket where the data lives"
    echo "For the name of the cluster, if running on [1] NA12878 sample, [2] CHM mix sample, "
    echo "  please indicate in the name of the cluster by including a string [1] 'na12878', [2]'chmmix', "
    echo "  so that the correct data is copied directly for you."
    exit 1
fi

PROJECT=$1
CLUSTER_NAME=$2
INITIALIZATION_BUCKET=$3

ON_NA12878=$(echo "$CLUSTER_NAME" | grep -i 'na12878' | wc -l | awk '{print $1}')
ON_CHM=$(echo "$CLUSTER_NAME" | grep -i 'chmmix' | wc -l | awk '{print $1}')
INTE_TEST=$(echo "$CLUSTER_NAME" | grep -i 'intetest' | wc -l | awk '{print $1}')

if [[ $ON_NA12878 != 1 && $ON_CHM != 1 && $INTE_TEST != 1 ]]; then
    echo "The name of the cluster to be created is:" "${CLUSTER_NAME}"
    echo "  which does not indicate in its name which sample it intends to run on"
    echo "If you intend running on [1] NA12878 sample, [2] CHM mix sample, "
    echo "  please type exactly '-na12878' or '-chmmix', "
    echo "  so that the correct data is copied directly for you."
    echo "Otherwise no data will be copied."
    read postfix
    CLUSTER_NAME+="$postfix"
    echo "Name of cluster is now:"
    echo "$CLUSTER_NAME"
fi

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
    --initialization-actions gs://${INITIALIZATION_BUCKET}/initialization_scripts/dataproc_initialization_new_bucket_copy_ref_img_custom_data.sh \
    --initialization-action-timeout 60m

