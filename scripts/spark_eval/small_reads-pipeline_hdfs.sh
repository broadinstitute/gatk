#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on small data in HDFS.

. utils.sh

time_gatk "ReadsPipelineSpark -I hdfs:///user/$USER/small_spark_eval/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam -O hdfs://${HDFS_HOST_PORT}/user/$USER/small_spark_eval/out/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.vcf -R hdfs:///user/$USER/small_spark_eval/human_g1k_v37.20.21.2bit --known-sites hdfs://${HDFS_HOST_PORT}/user/$USER/small_spark_eval/dbsnp_138.b37.20.21.vcf -pairHMM AVX_LOGLESS_CACHING --max-reads-per-alignment-start 10" 1 8 32g 4g
