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

REFERENCE_LOCATION="$MASTER_NODE"/reference/Homo_sapiens_assembly19.fasta
TWOBIT_REFERENCE_LOCATION="$MASTER_NODE"/reference/Homo_sapiens_assembly19.2bit
echo "Assuming reference: " "$REFERENCE_LOCATION"
echo "Assuming 2-bit reference: " "$TWOBIT_REFERENCE_LOCATION"

cd "$GATK_DIR" 

./gatk-launch CallVariantsFromAlignedContigsSpark \
    --inputAlignments "$PROJECT_OUTPUT_DIR"/aligned_assemblies \
    --inputAssemblies "$PROJECT_OUTPUT_DIR"/assembly_0 \
    --outputPath "$PROJECT_OUTPUT_DIR"/variants \
    -R "$TWOBIT_REFERENCE_LOCATION" \
    -fastaReference "$REFERENCE_LOCATION" \
    -- \
    --sparkRunner GCS \
    --cluster "$CLUSTER_NAME" \
    --driver-memory 30G \
    --executor-memory 8G \
    --conf spark.yarn.executor.memoryOverhead=5000
