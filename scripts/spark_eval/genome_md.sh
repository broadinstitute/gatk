#!/usr/bin/env bash

# Run Mark Duplicates on genome data on a Spark cluster.

. utils.sh

time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -O hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded --shardedOutput true -XL NC_007605" 64 1 4g 4g
