#!/bin/bash

set -eu

if [[ -z ${GATK_DIR+x} || -z ${CLUSTER_NAME+x} || -z ${MASTER_NODE+x} || -z ${PROJECT_OUTPUT_DIR+x} ]]; then 
    if [[ "$#" -ne 3 ]]; then
        echo "Please provide: local directory of GATK build, cluster name, output dir on the cluster"
    exit 1
    fi
    GATK_DIR=$1
    CLUSTER_NAME=$2
    OUTPUT_DIR=$3
    MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"
    PROJECT_OUTPUT_DIR="$MASTER_NODE"/"$OUTPUT_DIR"
fi

REFERENCE_IMG_LOCATION="$MASTER_NODE"/reference/Homo_sapiens_assembly19.fasta.64.img

echo "Assuming reference image: " "$REFERENCE_IMG_LOCATION"

cd "$GATK_DIR" 

./gatk-launch AlignAssembledContigsSpark \
    --bwamemIndexImage "$REFERENCE_IMG_LOCATION" \
    --inputFile "$PROJECT_OUTPUT_DIR"/assembly_0 \
    -O "$PROJECT_OUTPUT_DIR"/aligned_assemblies \
    -- \
    --sparkRunner GCS \
    --cluster "$CLUSTER_NAME" \
    --driver-memory 30G \
    --executor-memory 30G \
    --num-executors 20 \
    --conf spark.yarn.executor.memoryOverhead=16000
