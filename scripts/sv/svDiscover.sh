#!/bin/bash

set -eu

if [[ "$#" -lt 6 ]]; then
    echo -e "Please provide:"
    echo -e "  [1] local directory of GATK build (required)"
    echo -e "  [2] cluster name (required)"
    echo -e "  [3] output directory on the cluster (HDFS,required)"
    echo -e "  [4] absolute path to the BAM on the cluster (index assumed accompanying the bam) (HDFS,required)"
    echo -e "  [5] absolute path to the reference fasta on the cluster (2bit is assumed accompanying with same basename and extension \".2bit\", skiplist with extension \".kill.intervals\") (HDFS,required)"
    echo -e "  [6] absolute path to the reference index image on each worker node's local file system (required)"
    echo -e "Example:"
    echo -e " bash svDiscover.sh \\"
    echo -e "      ~/GATK/gatk \\"
    echo -e "      my-test-cluster \\"
    echo -e "      test-sample \\"
    echo -e "      /data/NA12878_test.bam \\"
    echo -e "      /reference/Homo_sapiens_assembly38.fasta"
    exit 1
fi

GATK_DIR=$1
CLUSTER_NAME=$2

MASTER_NODE="hdfs://""$CLUSTER_NAME""-m:8020"

PROJECT_OUTPUT_DIR="$MASTER_NODE"/$3
INPUT_BAM="$MASTER_NODE"/$4
REF_FASTA="$MASTER_NODE"/$5
REF_2BIT=$(echo $REF_FASTA | sed 's/.fasta$/.2bit/')
INTERVAL_KILL_LIST=$(echo $REF_FASTA | sed 's/.fasta$/.kill.intervals/')
KMER_KILL_LIST=$(echo $REF_FASTA | sed 's/.fasta$/.kill.kmers/')
ALT_MOL_FILE=$(echo $REF_FASTA | sed 's/.fasta$/.mol')
REF_INDEX_IMG=$6

export GATK_DIR
export CLUSTER_NAME
export MASTER_NODE
export PROJECT_OUTPUT_DIR
export INPUT_BAM
export REF_FASTA
export REF_2BIT
export REF_INDEX_IMG
export INTERVAL_KILL_LIST
export KMER_KILL_LIST
export ALT_MOL_FILE

./scanBam.sh

./assembly.sh

./alignAssembly.sh

./discoverVariants.sh
