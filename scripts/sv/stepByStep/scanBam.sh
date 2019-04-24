#!/bin/bash

set -eu

if [[ "$#" -ne 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (required)"
    echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (required)"
    echo -e "  [5] absolute path to the 2bit reference on the cluster (required)"
    echo -e "  [6] absolute path to the reference index image on each worker (required)"
    exit 1
fi

GATK_DIR="$1"
CLUSTER_NAME="$2"
MASTER_NODE="hdfs://${CLUSTER_NAME}-m:8020"
PROJECT_OUTPUT_DIR="${MASTER_NODE}$3"
INPUT_BAM="${MASTER_NODE}$4"
REF_TWOBIT="${MASTER_NODE}$5"
REF_INDEX_IMAGE="$6"
INTERVAL_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.intervals/')
KMER_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.kmers/')
ALTS_KILL_LIST=$(echo "${REF_TWOBIT}" | sed 's/.2bit$/.kill.alts/')

"${GATK_DIR}/gatk" FindBreakpointEvidenceSpark \
    -I "${INPUT_BAM}" \
    -O "${PROJECT_OUTPUT_DIR}/assemblies.sam" \
    --aligner-index-image "${REF_INDEX_IMAGE}" \
    --exclusion-intervals "${INTERVAL_KILL_LIST}" \
    --kmers-to-ignore "${KMER_KILL_LIST}" \
    --cross-contigs-to-ignore "${ALTS_KILL_LIST}" \
    --breakpoint-intervals "${PROJECT_OUTPUT_DIR}/intervals" \
    --fastq-dir "${PROJECT_OUTPUT_DIR}/fastq" \
    -- \
    --spark-runner GCS \
    --cluster "${CLUSTER_NAME}" \
    --num-executors 20 \
    --driver-memory 30G \
    --executor-memory 30G \
    --conf spark.executor.memoryOverhead=5000
