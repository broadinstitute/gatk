#!/bin/bash

set -eu

if [[ "$#" -ne 4 ]]; then
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (HDFS,required)"
    echo -e "  [4] absolute path to the reference fasta on the cluster (2bit is assumed accompanying with same basename and extension \".2bit\") (HDFS,required)"

    exit 1
fi

GATK_DIR="$1"
CLUSTER_NAME="$2"
MASTER_NODE="hdfs://${CLUSTER_NAME}-m:8020"
PROJECT_OUTPUT_DIR="${MASTER_NODE}$3"
REF_FASTA="${MASTER_NODE}$4"
REF_TWOBIT=$(echo "${REF_FASTA}" | sed 's/.fasta$/.2bit/')

"${GATK_DIR}/gatk-launch" CallVariantsFromAlignedContigsSAMSpark \
    -I "${PROJECT_OUTPUT_DIR}/assemblies.sam" \
    --outputPath "${PROJECT_OUTPUT_DIR}/variants" \
    -R "${REF_TWOBIT}" \
    --fastaReference "${REF_FASTA}" \
    -- \
    --sparkRunner GCS \
    --cluster "${CLUSTER_NAME}" \
    --driver-memory 30G \
    --executor-memory 8G \
    --conf spark.yarn.executor.memoryOverhead=5000
