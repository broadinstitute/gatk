#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on exome data in HDFS.

. utils.sh

time_gatk "ReadsPipelineSpark -I hdfs:///user/$USER/exome_spark_eval/NA12878.ga2.exome.maq.raw.bam -O hdfs://${HDFS_HOST_PORT}/user/$USER/exome_spark_eval/out/NA12878.ga2.exome.maq.raw.vcf -R hdfs:///user/$USER/exome_spark_eval/Homo_sapiens_assembly18.2bit --known-sites hdfs://${HDFS_HOST_PORT}/user/$USER/exome_spark_eval/dbsnp_138.hg18.vcf -pairHMM AVX_LOGLESS_CACHING --max-reads-per-alignment-start 10" 20 7 28g 4g

# Notes
# 20 executors - 2 per node (this is run on a 10 node cluster of n1-standard-16, each with 16 cores, 60g)
# 7 cores each - total 140, gives 20 left - 2 per node
# 28g each - total 560, gives 40g left - 4g per node