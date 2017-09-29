#!/usr/bin/env bash

# Run BQSR on genome data on a GCS Dataproc cluster. Data is in HDFS.

. utils.sh

time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded --shardedOutput true -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit --knownSites hdfs://${GCS_CLUSTER}-m:8020/user/$USER/q4_spark_eval/dbsnp_138.b37.vcf --joinStrategy OVERLAPS_PARTITIONER" 10 8 42g 4g
