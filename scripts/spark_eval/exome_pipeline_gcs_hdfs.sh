#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on exome data on a GCS Dataproc cluster. Data is in HDFS.

. utils.sh

time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/exome_spark_eval/NA12878.ga2.exome.maq.raw.bam -O hdfs:///user/$USER/exome_spark_eval/out/markdups-sharded --shardedOutput true" 48 1 4g 4g
# Note need to include host/port on following line. Need to investigate how to remove.
time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/exome_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/exome_spark_eval/out/bqsr-sharded --shardedOutput true -R hdfs:///user/$USER/exome_spark_eval/Homo_sapiens_assembly18.2bit --knownSites hdfs://${GCS_CLUSTER}-m:8020/user/$USER/exome_spark_eval/dbsnp_138.hg18.vcf" 4 8 32g 4g
time_gatk "HaplotypeCallerSpark -I hdfs:///user/$USER/exome_spark_eval/out/bqsr-sharded -R hdfs:///user/$USER/exome_spark_eval/Homo_sapiens_assembly18.2bit -O hdfs://${GCS_CLUSTER}-m:8020/user/$USER/exome_spark_eval/out/NA12878.ga2.exome.maq.raw.vcf -pairHMM AVX_LOGLESS_CACHING -maxReadsPerAlignmentStart 10" 32 1 6g 4g