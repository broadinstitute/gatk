#!/bin/bash

# This script creates a Google Dataproc cluster used for running the GATK-SV pipeline.

set -eu

if [[ "$#" -lt 8 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] project name (required)"
    echo -e "  [3] cluster name (required)"
    echo -e "  [4] maximum life time of cluster (required)"
    echo -e "  [5] maximum idling time before cluster self-termination begins (required)"
    echo -e "  [6] GCS directory containing the correct reference & indices (required)"
    echo -e "  [7] GCS directory containing the BAM and its index (required), or the GCS url of the BAM file, and"
    echo -e "either when a local initialization script is to be uploaded to a bucket:"
    echo -e "  [8] local initialization script"
    echo -e "  [9] GCS path to upload initialization script to (required as per Google)"
    echo -e "or when the initialization script is already in the cloud: "
    echo -e "  [8] path to initialization script on GCS"
    echo -e "\nExample when you have a local custom initialization script:"
    echo -e "  bash create_cluster.sh \\"
    echo -e "       <GATK-DIR> \\"
    echo -e "       <PROJECT-NAME> \\"
    echo -e "       <CLUSTER-NAME> \\"
    echo -e "       <CLUSTER-MAX-IDLE-MINUTES> \\"
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
    echo -e "       <CLUSTER-MAX-IDLE-MINUTES> \\"
    echo -e "       <BUCKET-DIR-STORING-REFERENCE> \\"
    echo -e "       <BUCKET-DIR-STORING-BAM> \\"
    echo -e "       <PATH-TO-INIT-SCRIPT-ON-GCS>"
    exit 1
fi

GATK_DIR=$1
PROJECT=$2
CLUSTER_NAME=$3
# dataproc beta feature for specifying a time frame when a cluster should terminate itself
MAX_LIFE_HOURS=$4
MAX_IDLE_MINUTES=$5
REF_DIR=$6
SAMP_INPUT=$7

INIT_ACTION=""
if [[ "$#" -eq 8 ]]; then
    INIT_ACTION=$8
else
    LOCAL_INIT_SCRIPT=$8
    INIT_SCRIPT_BUCKET=$9
    gsutil cp $LOCAL_INIT_SCRIPT ${INIT_SCRIPT_BUCKET}
    INIT_ACTION=${INIT_SCRIPT_BUCKET}
    if [[ "$INIT_ACTION" =~ .+/$ ]]; then
        INIT_ACTION+=$(basename $8)
    else
        INIT_ACTION+="/"
        INIT_ACTION+=$(basename $8)
    fi
fi

ZONE=${CLOUDSDK_COMPUTE_ZONE:-"us-central1-a"}
if [ "${CLOUDSDK_COMPUTE_ZONE:-""}" == "" ]; then
    echo "Using zone=${ZONE} because CLOUDSDK_COMPUTE_ZONE is not set"
else
    # note: not a typo, letting user know that CLOUDSDK_COMPUTE_ZONE
    # is the environmental variable being used
    echo "Using zone=CLOUDSDK_COMPUTE_ZONE=${ZONE}"
fi

# experimental feature: allow setting number of workers
# make 10 worker by default, but allow overload by setting env variable
NUM_SV_WORKERS=${NUM_SV_WORKERS:-10}
# make *no* preemptible workers by default, but allow overload by
# setting env variable
NUM_SV_PREEMPTIBLE_WORKERS=${NUM_SV_PREEMPTIBLE_WORKERS:-0}
SV_MASTER_BOOT_DISK_SIZE=${SV_MASTER_BOOT_DISK_SIZE:-500}
SV_WORKER_BOOT_DISK_SIZE=${SV_WORKER_BOOT_DISK_SIZE:-500}
SV_BOOT_DISK_TYPE=${SV_BOOT_DISK_TYPE:-"pd-standard"}

gcloud beta dataproc clusters create ${CLUSTER_NAME} \
    --zone ${ZONE} \
    --master-machine-type n1-highmem-8 \
    --worker-machine-type n1-highmem-16 \
    --master-boot-disk-size $SV_MASTER_BOOT_DISK_SIZE \
    --master-boot-disk-type $SV_BOOT_DISK_TYPE \
    --worker-boot-disk-size $SV_WORKER_BOOT_DISK_SIZE \
    --worker-boot-disk-type $SV_BOOT_DISK_TYPE \
    --preemptible-worker-boot-disk-size=$SV_WORKER_BOOT_DISK_SIZE \
    --preemptible-worker-boot-disk-type=$SV_BOOT_DISK_TYPE \
    --num-workers ${NUM_SV_WORKERS} \
    --num-preemptible-workers ${NUM_SV_PREEMPTIBLE_WORKERS} \
    --num-worker-local-ssds 1 \
    --metadata "reference=$REF_DIR" \
    --metadata "sample=$SAMP_INPUT" \
    --image-version 1.3 \
    --project ${PROJECT} \
    --initialization-actions ${INIT_ACTION} \
    --initialization-action-timeout 10m \
    --max-age ${MAX_LIFE_HOURS} \
    --max-idle ${MAX_IDLE_MINUTES} \
    --properties yarn:yarn.resourcemanager.am.max-attempts=5,mapred:mapreduce.map.maxattempts=8,mapred:mapreduce.reduce.maxattempts=8,spark:spark.task.maxFailures=8,spark:spark.stage.maxConsecutiveAttempts=8

MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"

if [[ ! $REF_DIR =~ .+/$ ]]; then
    REF_DIR+="/"
fi

${GATK_DIR}/gatk ParallelCopyGCSDirectoryIntoHDFSSpark \
    --input-gcs-path "$REF_DIR" \
    --output-hdfs-directory "$MASTER_NODE"/reference \
    -- \
    --spark-runner GCS \
    --cluster "$CLUSTER_NAME" \
    --project $PROJECT

if [[ $SAMP_INPUT =~ .+\.bam$ ]]; then
    SAMP_DIR=`dirname $SAMP_INPUT`
    SAMP_GLOB_ARG="--input-file-glob `basename $SAMP_INPUT`*"
else
    SAMP_DIR=$SAMP_INPUT
    SAMP_GLOB_ARG=""
fi

if [[ ! $SAMP_DIR =~ .+/$ ]]; then
    SAMP_DIR+="/"
fi

set -f
${GATK_DIR}/gatk ParallelCopyGCSDirectoryIntoHDFSSpark \
    --input-gcs-path "$SAMP_DIR" \
    $SAMP_GLOB_ARG \
    --output-hdfs-directory "$MASTER_NODE"/data \
    -- \
    --spark-runner GCS \
    --cluster "$CLUSTER_NAME" \
    --project $PROJECT
