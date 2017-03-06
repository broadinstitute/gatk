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

cd "$GATK_DIR" 

./gatk-launch FindBreakpointEvidenceSpark \
    -R "$MASTER_NODE"/reference/Homo_sapiens_assembly19.fasta \
    -I "$MASTER_NODE"/data/NA12878_PCR-_30X.bam \
    -O "$PROJECT_OUTPUT_DIR"/fastq \
    --exclusionIntervals "$MASTER_NODE"/reference/GRCh37.kill.intervals \
    --kmersToIgnore "$MASTER_NODE"/reference/Homo_sapiens_assembly38.dups \
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
