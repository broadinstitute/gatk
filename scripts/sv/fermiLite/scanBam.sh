#!/bin/bash

set -eu

if [[ "$#" -ne 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (required)"
    echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (required)"
    echo -e "  [5] absolute path to the reference fasta on the cluster (required)"
    echo -e "  [6] absolute path to the reference index image on each worker (required)"
    exit 1
fi

GATK_DIR="$1"
CLUSTER_NAME="$2"
MASTER_NODE="hdfs://${CLUSTER_NAME}-m:8020"
PROJECT_OUTPUT_DIR="${MASTER_NODE}$3"
INPUT_BAM="${MASTER_NODE}$4"
REF_FASTA="${MASTER_NODE}$5"
REF_INDEX_IMAGE="$6"
INTERVAL_KILL_LIST=$(echo "${REF_FASTA}" | sed 's/.fasta$/.kill.intervals/')
KMER_KILL_LIST=$(echo "${REF_FASTA}" | sed 's/.fasta$/.kill.kmers/')

"${GATK_DIR}/gatk-launch" FindBreakpointEvidenceSpark \
    -I "${INPUT_BAM}" \
    -O "${PROJECT_OUTPUT_DIR}/assemblies.sam" \
    --alignerIndexImage "${REF_INDEX_IMAGE}" \
    --exclusionIntervals "${INTERVAL_KILL_LIST}" \
    --kmersToIgnore "${KMER_KILL_LIST}" \
    --breakpointIntervals "${PROJECT_OUTPUT_DIR}/intervals" \
    --fastqDir "${PROJECT_OUTPUT_DIR}/fastq" \
    -- \
    --sparkRunner GCS \
    --cluster "${CLUSTER_NAME}" \
    --num-executors 20 \
    --driver-memory 30G \
    --executor-memory 30G \
    --conf spark.yarn.executor.memoryOverhead=5000
