#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on genome data in HDFS.

. utils.sh

time_gatk "ReadsPipelineSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878-no-NC_007605.bam -O hdfs://${HDFS_HOST_PORT}/user/$USER/q4_spark_eval/out/WGS-G94982-NA12878.vcf -R hdfs://${HDFS_HOST_PORT}/user/$USER/q4_spark_eval/human_g1k_v37.fasta --known-sites hdfs://${HDFS_HOST_PORT}/user/$USER/q4_spark_eval/dbsnp_138.b37.vcf.gz -pairHMM AVX_LOGLESS_CACHING --maxReadsPerAlignmentStart 10" 40 7 20g 8g
