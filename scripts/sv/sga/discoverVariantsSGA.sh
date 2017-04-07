#!/bin/bash

set -eu

# either have the required variables set in the running environment, or provide them 
#  in command line
if [[ -z ${GATK_DIR+x} || -z ${CLUSTER_NAME+x} || -z ${PROJECT_OUTPUT_DIR+x} || -z ${REF_FASTA+x} || -z ${REF_2BIT+x} ]]; then 
    if [[ "$#" -ne 4 ]]; then
        echo -e "Please provide:"
        echo -e "  [1] local directory of GATK build (required)"
        echo -e "  [2] cluster name (required)"
        echo -e "  [3] output directory on the cluster (required)"
        echo -e "  [4] absolute path to the reference fasta on the cluster (required)"
    exit 1
    fi
    GATK_DIR=$1
    CLUSTER_NAME=$2
    MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"
    PROJECT_OUTPUT_DIR="$MASTER_NODE"/$3
    REF_FASTA="$MASTER_NODE"/$4
    REF_2BIT=$(echo $REF_FASTA | sed 's/.fasta$/.2bit/')
fi

cd "$GATK_DIR" 

./gatk-launch DiscoverVariantsFromAlignedSGAContigsSpark \
    --inputAlignments "$PROJECT_OUTPUT_DIR"/aligned_assemblies \
    --inputAssemblies "$PROJECT_OUTPUT_DIR"/assembly_0 \
    --outputPath "$PROJECT_OUTPUT_DIR"/variants \
    -R "$REF_2BIT" \
    -fastaReference "$REF_FASTA" \
    -- \
    --sparkRunner GCS \
    --cluster "$CLUSTER_NAME" \
    --driver-memory 30G \
    --executor-memory 8G \
    --conf spark.yarn.executor.memoryOverhead=5000
