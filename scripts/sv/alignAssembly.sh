#!/bin/bash

set -eu

# either have the required variables set in the running environment, or provide them 
#  in command line
if [[ -z ${GATK_DIR+x} || -z ${CLUSTER_NAME+x} || -z ${PROJECT_OUTPUT_DIR+x} || -z ${REF_INDEX_IMG+x} ]]; then 
    if [[ "$#" -ne 4 ]]; then
        echo -e "Please provide:"
        echo -e "  [1] local directory of GATK build (required)"
        echo -e "  [2] cluster name (required)"
        echo -e "  [3] output directory on the cluster (required)"
        echo -e "  [4] absolute path to the reference index image on each worker node's local file system (required)"
    exit 1
    fi
    GATK_DIR=$1
    CLUSTER_NAME=$2
    OUTPUT_DIR=$3
    REF_INDEX_IMG=$4
    MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"
    PROJECT_OUTPUT_DIR="$MASTER_NODE"/"$OUTPUT_DIR"
fi

cd "$GATK_DIR" 

./gatk-launch AlignAssembledContigsSpark \
    --bwamemIndexImage "$REF_INDEX_IMG" \
    --inputFile "$PROJECT_OUTPUT_DIR"/assembly_0 \
    -O "$PROJECT_OUTPUT_DIR"/aligned_assemblies \
    -- \
    --sparkRunner GCS \
    --cluster "$CLUSTER_NAME" \
    --driver-memory 30G \
    --executor-memory 30G \
    --num-executors 20 \
    --conf spark.yarn.executor.memoryOverhead=16000
