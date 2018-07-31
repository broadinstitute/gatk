#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on exome data in HDFS.

. utils.sh

time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/exome_spark_eval/NA12878.ga2.exome.maq.raw.bam -O hdfs:///user/$USER/exome_spark_eval/out/markdups-sharded --sharded-output true" 96 1 4g 4g
time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/exome_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/exome_spark_eval/out/bqsr-sharded --sharded-output true -R hdfs:///user/$USER/exome_spark_eval/Homo_sapiens_assembly18.2bit --known-sites hdfs://${HDFS_HOST_PORT}/user/$USER/exome_spark_eval/dbsnp_138.hg18.vcf" 8 8 32g 4g
time_gatk "HaplotypeCallerSpark -I hdfs:///user/$USER/exome_spark_eval/out/bqsr-sharded -R hdfs:///user/$USER/exome_spark_eval/Homo_sapiens_assembly18.2bit -O hdfs://${HDFS_HOST_PORT}/user/$USER/exome_spark_eval/out/NA12878.ga2.exome.maq.raw.vcf -pairHMM AVX_LOGLESS_CACHING --max-reads-per-alignment-start 10" 64 1 6g 4g