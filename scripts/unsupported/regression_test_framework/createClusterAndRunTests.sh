#!/usr/bin/env bash

###############################################################################
# Set up
###############################################################################

INPUT_BAM="${1}"

# Env variables
export DATAPROC_CLUSTER_NAME=cluster-$USER
export CLOUDSDK_CORE_PROJECT=broad-gatk-collab


# Run a dataproc cluster with 5 worker nodes, 200GB of disk, 8 vcores
gcloud dataproc clusters create ${DATAPROC_CLUSTER_NAME}-seven --subnet default --zone us-central1-b --master-machine-type n1-standard-8 --master-boot-disk-size 500 --num-workers 7 --worker-machine-type n1-highmem-8 --worker-boot-disk-size 2000 --image-version 1.2 --project broad-dsde-dev --max-age 1d

gcloud dataproc clusters create ${DATAPROC_CLUSTER_NAME}-fourteen --subnet default --zone us-central1-b --master-machine-type n1-standard-8 --master-boot-disk-size 500 --num-workers 14 --worker-machine-type n1-highmem-4 --worker-boot-disk-size 2000 --image-version 1.2 --project broad-dsde-dev --max-age 1d

###############################################################################
# Benchmark performance
###############################################################################

haplotype-caller-small() {
    CLASS=$1
    ARGS=$2
    LOG=logs/${CLASS}_$(date +%Y%m%d_%H%M%S).log
    RESULTS_CSV=results/run.csv
    gcloud dataproc jobs submit spark --cluster $DATAPROC_CLUSTER_NAME \
     --project $CLOUDSDK_CORE_PROJECT \
     --properties spark.executor.instances=5,spark.executor.cores=8,spark.executor.memory=4g,spark.driver.memory=4g,spark.dynamicAllocation.enabled=false \
     --jars gs://disq-tom-testdata/jars/disq-benchmarks-0.0.1-SNAPSHOT.jar \
     --class $CLASS \
     -- \
     $ARGS 2>&1 | tee /dev/tty > $LOG

     ./gatk HaplotypeCallerSpark \
    -I gs://my-gcs-bucket/path/to/input.bam \
    -O gs://my-gcs-bucket/path/to/output.bam \
    -R
    -- \
    --spark-runner GCS --cluster $DATAPROC_CLUSTER_NAME \
    --num-executors 7 --executor-cores 8 --executor-memory 52g \
    --conf spark.yarn.executor.memoryOverhead=600
    RC=$?
    DURATION_SEC=$(grep 'Time taken' $LOG | grep -Eo "[0-9]+")
    echo "$CLASS,$ARGS,$RC,$DURATION_SEC" >> $RESULTS_CSV
}

count-reads-sparkbam() {
    ARGS=$1
    gcloud dataproc jobs submit spark --cluster $DATAPROC_CLUSTER_NAME \
     --project $CLOUDSDK_CORE_PROJECT \
     --properties spark.executor.instances=5,spark.executor.cores=8,spark.executor.memory=4g,spark.driver.memory=4g,spark.dynamicAllocation.enabled=false \
     --jars gs://disq-tom-testdata/jars/$GOOGLE_CLOUD_NIO_JAR \
     --jar gs://disq-tom-testdata/jars/$CLI_JAR \
     -- \
     count-reads $ARGS
}

MID_SIZED_BAM=gs://disq-tom-testdata/1000genomes/ftp/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam

count-reads com.tom_e_white.disq.benchmarks.DisqCountReads "$MID_SIZED_BAM yarn false 134217728"
count-reads com.tom_e_white.disq.benchmarks.DisqCountReads "$MID_SIZED_BAM yarn true 134217728"
count-reads com.tom_e_white.disq.benchmarks.HadoopBamCountReads "$MID_SIZED_BAM yarn"

LARGE_SIZED_BAM=gs://disq-tom-testdata/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam

count-reads com.tom_e_white.disq.benchmarks.DisqCountReads "$LARGE_SIZED_BAM yarn false 134217728"
count-reads com.tom_e_white.disq.benchmarks.DisqCountReads "$LARGE_SIZED_BAM yarn true 134217728"
count-reads com.tom_e_white.disq.benchmarks.HadoopBamCountReads "$LARGE_SIZED_BAM yarn"
count-reads-sparkbam $LARGE_SIZED_BAM

# Before running against files in HDFS, ssh into the master node and run
# gsutil cp gs://disq-tom-testdata/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam - | hadoop fs -put - /tmp/sorted_final_merged.bam

count-reads com.tom_e_white.disq.benchmarks.DisqCountReads "hdfs://$DATAPROC_CLUSTER_NAME-m/tmp/sorted_final_merged.bam yarn false 134217728"
count-reads com.tom_e_white.disq.benchmarks.HadoopBamCountReads "hdfs://$DATAPROC_CLUSTER_NAME-m/tmp/sorted_final_merged.bam yarn"

###############################################################################
# Benchmark accuracy
###############################################################################

check-bam() {
    CLASS=$1
    ARGS=$2
    LOG=logs/${CLASS}_$(date +%Y%m%d_%H%M%S).log
    gcloud dataproc jobs submit spark --cluster $DATAPROC_CLUSTER_NAME \
     --project $CLOUDSDK_CORE_PROJECT \
     --properties spark.executor.instances=5,spark.executor.cores=8,spark.executor.memory=4g,spark.driver.memory=6g,spark.dynamicAllocation.enabled=false \
     --jars gs://disq-tom-testdata/jars/disq-benchmarks-0.0.1-SNAPSHOT.jar \
     --class $CLASS \
     -- \
     $ARGS 2>&1 | tee /dev/tty > $LOG
}

# 1kg
for bam in \
    gs://disq-tom-testdata/1000genomes/ftp/data/HG00096/alignment/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam \
    gs://disq-tom-testdata/1000genomes/ftp/data/HG00096/alignment/HG00096.unmapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam \
    gs://disq-tom-testdata/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.mapped.ILLUMINA.bwa.GBR.exome.20120522.bam \
    gs://disq-tom-testdata/1000genomes/ftp/data/HG00096/exome_alignment/HG00096.unmapped.ILLUMINA.bwa.GBR.exome.20120522.bam
do
    check-bam com.tom_e_white.disq.benchmarks.DisqCheckBam "$bam yarn false 134217728"
done

# giab
for bam in \
    gs://disq-tom-testdata/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg002_gr37_X.bam \
    gs://disq-tom-testdata/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/MtSinai_blasr_bam_GRCh37/hg002_gr37_Y.bam \
    gs://disq-tom-testdata/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam
do
    check-bam com.tom_e_white.disq.benchmarks.DisqCheckBam "$bam yarn false 134217728"
done