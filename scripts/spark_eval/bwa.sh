#!/usr/bin/env bash

# Run BwaSpark on genome data in HDFS.

. utils.sh

time_gatk "BwaSpark -I hdfs:///user/$USER/q4_spark_eval/G15512.HCC1954.RG.bam -O hdfs:///user/$USER/q4_spark_eval/out/bwa-sharded -R hdfs://${HDFS_HOST_PORT}/user/$USER/q4_spark_eval/human_g1k_v37.fasta --disableSequenceDictionaryValidation true --shardedOutput true" 20 4 46g 8g
