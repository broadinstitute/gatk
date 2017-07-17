#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on exome data on a GCS Dataproc cluster. Data is in HDFS.

. utils.sh

time_gatk "ReadsPipelineSpark -I hdfs:///user/$USER/exome_spark_eval/NA12878.ga2.exome.maq.raw.bam -O hdfs://${GCS_CLUSTER}-m:8020/user/$USER/exome_spark_eval/out/NA12878.ga2.exome.maq.raw.vcf -R hdfs:///user/$USER/exome_spark_eval/Homo_sapiens_assembly18.2bit --knownSites hdfs://${GCS_CLUSTER}-m:8020/user/$USER/exome_spark_eval/dbsnp_138.hg18.vcf -pairHMM AVX_LOGLESS_CACHING -maxReadsPerAlignmentStart 10" 4 8 32g 4g