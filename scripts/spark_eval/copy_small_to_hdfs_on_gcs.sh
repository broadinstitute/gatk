#!/usr/bin/env bash

# Copy small data to HDFS on a GCS cluster.

${GATK_HOME:-../..}/gatk-launch ParallelCopyGCSDirectoryIntoHDFSSpark \
    --inputGCSPath gs://broad-spark-eval-test-data/small/ \
    --outputHDFSDirectory hdfs://${GCS_CLUSTER}-m:8020/user/$USER/small_spark_eval \
    -apiKey $API_KEY \
    -- \
    --sparkRunner GCS \
    --cluster $GCS_CLUSTER