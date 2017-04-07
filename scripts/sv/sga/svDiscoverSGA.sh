#!/bin/bash

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] absolute path to the output directory on the cluster (HDFS,required)"
    echo -e "  [4] absolute path to the directory containing interleaved FASTQ files on the cluster (index assumed accompanying the bam) (HDFS,required)"
    echo -e "  [5] absolute path to the reference fasta on the cluster (2bit is assumed accompanying with same basename and extension \".2bit\", skiplist with extension \".kill.intervals\") (HDFS,required)"
    echo -e "  [6] absolute path to the reference index image on each worker node's local file system (required)"
    echo -e "Example:"
    echo -e " bash svDiscover.sh \\"
    echo -e "      ~/GATK/gatk \\"
    echo -e "      my-test-cluster \\"
    echo -e "      /test-sample \\"
    echo -e "      /test-sample/fastq \\"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta.img"
    exit 1
fi

GATK_DIR=$1
CLUSTER_NAME=$2

MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"

PROJECT_OUTPUT_DIR="$MASTER_NODE"/$3
FASTQ_DIR="$MASTER_NODE"/$4
REF_FASTA="$MASTER_NODE"/$5
REF_2BIT=$(echo $REF_FASTA | sed 's/.fasta$/.2bit/')
REF_INDEX_IMG=$6

export GATK_DIR
export CLUSTER_NAME
export PROJECT_OUTPUT_DIR
export REF_FASTA
export REF_2BIT
export REF_INDEX_IMG

baseDir=$(dirname -- "$0")
${baseDir}/assemblySGA.sh
${baseDir}/alignAssembly.sh
${baseDir}/discoverVariantsSGA.sh
