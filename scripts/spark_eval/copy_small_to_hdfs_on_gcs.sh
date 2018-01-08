#!/usr/bin/env bash

# Copy small data to HDFS on a GCS cluster.

${GATK_HOME:-../..}/gatk ParallelCopyGCSDirectoryIntoHDFSSpark \
    --input-gcs-path gs://broad-spark-eval-test-data/small/ \
    --output-hdfs-directory hdfs://${GCS_CLUSTER}-m:8020/user/$USER/small_spark_eval \
    -- \
    --spark-runner GCS \
    --cluster $GCS_CLUSTER