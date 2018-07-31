#!/usr/bin/env bash

# Copy genome data to HDFS on a GCS cluster.

${GATK_HOME:-../..}/gatk ParallelCopyGCSDirectoryIntoHDFSSpark \
    --input-gcs-path gs://broad-spark-eval-test-data/genome/ \
    --output-hdfs-directory hdfs://${GCS_CLUSTER}-m:8020/user/$USER/q4_spark_eval \
    -- \
    --spark-runner GCS \
    --cluster $GCS_CLUSTER