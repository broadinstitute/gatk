#!/bin/bash

set -eu

# either have the required variables set in the running environment, or provide them 
#  in command line
if [[ -z ${GATK_DIR+x} || -z ${CLUSTER_NAME+x} || -z ${PROJECT_OUTPUT_DIR+x} || -z ${INPUT_BAM+x} || -z ${REF_FASTA+x} || -z ${INTERVAL_KILL_LIST+x} || -z ${KMER_KILL_LIST+x} ]]; then 
    if [[ "$#" -ne 5 ]]; then
        echo -e "Please provide:"
        echo -e "  [1] local directory of GATK build (required)"
        echo -e "  [2] cluster name (required)"
        echo -e "  [3] output directory on the cluster (required)"
        echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (required)"
        echo -e "  [5] absolute path to the reference fasta on the cluster (required)"
        exit 1
    fi
    GATK_DIR=$1
    CLUSTER_NAME=$2

    MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"
    PROJECT_OUTPUT_DIR="$MASTER_NODE"/$3
    INPUT_BAM="$MASTER_NODE"/$4
    REF_FASTA="$MASTER_NODE"/$5
    INTERVAL_KILL_LIST=$(echo $REF_FASTA | sed 's/.fasta$/.kill.intervals/')
    KMER_KILL_LIST=$(echo $REF_FASTA | sed 's/.fasta$/.kill.kmers/')
fi

cd "$GATK_DIR" 

./gatk-launch FindBreakpointEvidenceSpark \
    -R "$REF_FASTA" \
    -I "$INPUT_BAM" \
    -O "$PROJECT_OUTPUT_DIR"/fastq \
    --exclusionIntervals "$INTERVAL_KILL_LIST" \
    --kmersToIgnore "$KMER_KILL_LIST" \
    --kmerIntervals "$PROJECT_OUTPUT_DIR"/kmerIntervals \
    --breakpointEvidenceDir "$PROJECT_OUTPUT_DIR"/evidence \
    --breakpointIntervals "$PROJECT_OUTPUT_DIR"/intervals \
    --qnameIntervalsMapped "$PROJECT_OUTPUT_DIR"/qnameIntervalsMapped \
    --qnameIntervalsForAssembly "$PROJECT_OUTPUT_DIR"/qnameIntervalsForAssembly\
    --maxFASTQSize 10000000 \
    -- \
    --sparkRunner GCS \
    --cluster "$CLUSTER_NAME" \
    --num-executors 20 \
    --driver-memory 30G \
    --executor-memory 30G \
    --conf spark.yarn.executor.memoryOverhead=5000
