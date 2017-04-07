#!/bin/bash

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (HDFS,required)"
    echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (HDFS,required)"
    echo -e "  [5] absolute path to the reference fasta on the cluster (2bit is assumed accompanying with same basename and extension \".2bit\", skiplist with extension \".kill.intervals\") (HDFS,required)"
    echo -e "  [6] absolute path to the reference index image on each worker node's local file system (required)"
    echo -e "Example:"
    echo -e " bash svDiscover.sh \\"
    echo -e "      ~/GATK/gatk \\"
    echo -e "      my-test-cluster \\"
    echo -e "      /test-sample \\"
    echo -e "      /data/NA12878_test.bam \\"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta.img"
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
REF_TWOBIT=$(echo "${REF_FASTA}" | sed 's/.fasta$/.2bit/')

"${GATK_DIR}/gatk-launch" StructuralVariationDiscoveryPipelineSpark \
    -I "${INPUT_BAM}" \
    -O "${PROJECT_OUTPUT_DIR}/variants/inv_del_ins.vcf" \
    -R "${REF_TWOBIT}" \
    --fastaReference "${REF_FASTA}" \
    --alignerIndexImage "${REF_INDEX_IMAGE}" \
    --exclusionIntervals "${INTERVAL_KILL_LIST}" \
    --kmersToIgnore "${KMER_KILL_LIST}" \
    --breakpointIntervals "${PROJECT_OUTPUT_DIR}/intervals" \
    --fastqDir "${PROJECT_OUTPUT_DIR}/fastq" \
    --contigSAMFile "${PROJECT_OUTPUT_DIR}/assemblies.sam" \
    -- \
    --sparkRunner GCS \
    --cluster "${CLUSTER_NAME}" \
    --num-executors 20 \
    --driver-memory 30G \
    --executor-memory 30G \
    --conf spark.yarn.executor.memoryOverhead=5000
