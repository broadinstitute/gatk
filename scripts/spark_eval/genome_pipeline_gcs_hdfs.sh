#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on genome data on a GCS Dataproc cluster. Data is in HDFS.

. utils.sh

time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878-no-NC_007605.bam -O hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded --shardedOutput true" 128 1 4g 4g
time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded --shardedOutput true -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit --knownSites hdfs://${GCS_CLUSTER}-m:8020/user/$USER/q4_spark_eval/dbsnp_138.b37.vcf" 10 8 42g 4g
time_gatk "HaplotypeCallerSpark -I hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit -O hdfs://${GCS_CLUSTER}-m:8020/user/$USER/q4_spark_eval/out/WGS-G94982-NA12878.vcf -pairHMM AVX_LOGLESS_CACHING -maxReadsPerAlignmentStart 10" 30 1 12g 8g
