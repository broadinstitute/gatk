#!/bin/bash

# This script creates a Google Dataproc cluster used for running the GATK-SV pipeline.

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] project name (required)"
    echo -e "  [3] cluster name (required)"
    echo -e "  [4] GCS directory containing the correct reference & indices (required)"
    echo -e "  [5] GCS directory containing the BAM and its index (required), and"
    echo -e "either when a local initialization script is to be uploaded to a bucket:"
    echo -e "  [6] local initialization script"
    echo -e "  [7] GCS path to upload initialization script to (required as per Google)"
    echo -e "or when the initialization script is already in the cloud: "
    echo -e "  [6] path to initialization script on GCS"
    echo -e "\nExample when you have a local custom initialization script:"
    echo -e "  bash create_cluster.sh \\"
    echo -e "       <GATK-DIR> \\"
    echo -e "       <PROJECT-NAME> \\"
    echo -e "       <CLUSTER-NAME> \\"
    echo -e "       <BUCKET-DIR-STORING-REFERENCE> \\"
    echo -e "       <BUCKET-DIR-STORING-BAM> \\"
    echo -e "       <LOCAL-CUSTOM-INIT-ACTION-SCRIPT> \\"
    echo -e "       <BUCKET-DIR-TO-STORE-INIT-SCRIPT>"
    echo -e "A default initialization script \"default_init.sh\" is distributed with this script."
    echo -e "If you don't need any custom initialization actions, use that."
    echo -e "\nExample 2 when a initialization script lives in a GCS bucket"
    echo -e "  bash create_cluster.sh \\"
    echo -e "       <GATK-DIR> \\"
    echo -e "       <PROJECT-NAME> \\"
    echo -e "       <CLUSTER-NAME> \\"
    echo -e "       <BUCKET-DIR-STORING-REFERENCE> \\"
    echo -e "       <BUCKET-DIR-STORING-BAM> \\"
    echo -e "       <PATH-TO-INIT-SCRIPT-ON-GCS>"
    exit 1
fi

GATK_DIR=$1
PROJECT=$2
CLUSTER_NAME=$3
REF_DIR=$4
SAMP_DIR=$5

INIT_ACTION=""
if [[ "$#" -eq 6 ]]; then
    INIT_ACTION=$6
else
    LOCAL_INIT_SCRIPT=$6
    INIT_SCRIPT_BUCKET=$7
    gsutil cp $LOCAL_INIT_SCRIPT ${INIT_SCRIPT_BUCKET}
    INIT_ACTION=${INIT_SCRIPT_BUCKET}
    if [[ "$INIT_ACTION" =~ .+/$ ]]; then
        INIT_ACTION+=$(basename $6)
    else
        INIT_ACTION+="/"
        INIT_ACTION+=$(basename $6)
    fi
fi

gcloud dataproc clusters create ${CLUSTER_NAME} \
    --zone us-central1-a \
    --master-machine-type n1-highmem-8 \
    --worker-machine-type n1-highmem-16 \
    --master-boot-disk-size 500 \
    --worker-boot-disk-size 500 \
    --num-workers 10 \
    --num-worker-local-ssds 1 \
    --metadata "reference=$REF_DIR" \
    --metadata "sample=$SAMP_DIR" \
    --image-version 1.1 \
    --project ${PROJECT} \
    --initialization-actions ${INIT_ACTION} \
    --initialization-action-timeout 10m

MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"

if [[ ! $REF_DIR =~ .+/$ ]]; then
    REF_DIR+="/"
fi

${GATK_DIR}/gatk-launch ParallelCopyGCSDirectoryIntoHDFSSpark \
    --inputGCSPath "$REF_DIR" \
    --outputHDFSDirectory "$MASTER_NODE"/reference \
    -- \
    --sparkRunner GCS \
    --cluster "$CLUSTER_NAME" \

if [[ ! $SAMP_DIR =~ .+/$ ]]; then
    SAMP_DIR+="/"
fi

${GATK_DIR}/gatk-launch ParallelCopyGCSDirectoryIntoHDFSSpark \
    --inputGCSPath "$SAMP_DIR" \
    --outputHDFSDirectory "$MASTER_NODE"/data \
    -- \
    --sparkRunner GCS \
    --cluster "$CLUSTER_NAME" \
