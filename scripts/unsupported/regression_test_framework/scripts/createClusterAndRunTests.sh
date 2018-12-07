#!/usr/bin/env bash

###############################################################################
# Set up
###############################################################################

INPUT_BAM="${1}"
INTERVAL_LIST_FILE="${2}"
OUTPUT_BAM_LOC="${3}"

# Env variables
export DATAPROC_CLUSTER_NAME=cluster-$USER
export CLOUDSDK_CORE_PROJECT=broad-gatk-collab

# Create a dataproc cluster with 7 worker nodes, 200GB of disk, 8 vcores
gcloud beta dataproc clusters create ${DATAPROC_CLUSTER_NAME}-seven --max-age=1h --zone us-central1-b --master-machine-type n1-standard-8 --master-boot-disk-size 500 --num-workers 7 --worker-machine-type n1-highmem-8 --worker-boot-disk-size 2000 --image-version 1.2 --project broad-dsde-dev

#gcloud dataproc clusters create ${DATAPROC_CLUSTER_NAME}-fourteen --subnet default --zone us-central1-b --master-machine-type n1-standard-8 --master-boot-disk-size 500 --num-workers 14 --worker-machine-type n1-highmem-4 --worker-boot-disk-size 2000 --image-version 1.2 --project broad-dsde-dev --max-age 1d

###############################################################################
# Benchmark performance
###############################################################################

haplotype-call-small() {
      ref_fasta="gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta"

    LOG=logs/${CLASS}_$(date +%Y%m%d_%H%M%S).log
    RESULTS_CSV=results/run.csv
    /Users/emeryj/hellbender/gatk/gatk HaplotypeCallerSpark \
      -I $1 \
      --interval-padding 500 \
      -L $2 \
      -R $ref_fasta \
      -O $3 \
      -- \
      --spark-runner GCS --cluster ${DATAPROC_CLUSTER_NAME}-seven \
      --project broad-dsde-dev \
      --num-executors 7 --executor-cores 8 --executor-memory 4g \
      --conf spark.yarn.executor.memoryOverhead=600
    RC=$?
    DURATION_SEC=$(grep 'Time taken' $LOG | grep -Eo "[0-9]+")
    echo "$1,$RC,$DURATION_SEC" >> $RESULTS_CSV
}
#
#for bam in \
#    gs://disq-tom-testdata/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg002_gr37_X.bam
#do
    haplotype-call-small $INPUT_BAM $INTERVAL_LIST_FILE $OUTPUT_BAM_LOC
#done