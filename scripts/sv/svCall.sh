#!/bin/bash

set -eu

if [[ "$#" -lt 3 ]]; then
    echo "Please provide: local directory of GATK build, cluster name, output dir on the cluster"
    echo "                reference location, "
    exit 1
else
    echo "assuming default reference, sample and reference image locations (see below)"
fi

GATK_DIR=$1
CLUSTER_NAME=$2
MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"
OUTPUT_DIR=$3

PROJECT_OUTPUT_DIR="$MASTER_NODE"/"$OUTPUT_DIR"

export GATK_DIR
export CLUSTER_NAME
export MASTER_NODE
export PROJECT_OUTPUT_DIR


./scanBam.sh

./assembly.sh

./alignAssembly.sh

./callVariants.sh
