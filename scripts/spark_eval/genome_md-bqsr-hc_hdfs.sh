#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on genome data in HDFS.

. utils.sh

time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878-no-NC_007605.bam -O hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded --sharded-output true" 256 1 4g 4g
time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded --sharded-output true -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit --known-sites hdfs://${HDFS_HOST_PORT}/user/$USER/q4_spark_eval/dbsnp_138.b37.vcf" 20 8 42g 4g
time_gatk "HaplotypeCallerSpark -I hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit -O hdfs://${HDFS_HOST_PORT}/user/$USER/q4_spark_eval/out/WGS-G94982-NA12878.vcf -pairHMM AVX_LOGLESS_CACHING --max-reads-per-alignment-start 10" 60 1 12g 8g
