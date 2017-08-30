#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on exome data in HDFS.

. utils.sh

time_gatk "ReadsPipelineSpark -I hdfs:///user/$USER/exome_spark_eval/NA12878.ga2.exome.maq.raw.bam -O hdfs://${HDFS_HOST_PORT}/user/$USER/exome_spark_eval/out/NA12878.ga2.exome.maq.raw.vcf -R hdfs:///user/$USER/exome_spark_eval/Homo_sapiens_assembly18.2bit --knownSites hdfs://${HDFS_HOST_PORT}/user/$USER/exome_spark_eval/dbsnp_138.hg18.vcf -pairHMM AVX_LOGLESS_CACHING -maxReadsPerAlignmentStart 10" 8 8 32g 4g