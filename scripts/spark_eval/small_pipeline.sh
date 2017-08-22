#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on small data on a Spark cluster.

. utils.sh

time_gatk "MarkDuplicatesSpark -I hdfs:///user/$USER/small_spark_eval/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam -O hdfs:///user/$USER/small_spark_eval/out/markdups-sharded --shardedOutput true" 8 1 4g 4g
time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/small_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/small_spark_eval/out/bqsr-sharded --shardedOutput true -R hdfs:///user/$USER/small_spark_eval/human_g1k_v37.20.21.2bit --knownSites hdfs:///user/$USER/small_spark_eval/dbsnp_138.b37.20.21.vcf -L 20:10000000-10100000" 1 8 32g 4g
time_gatk "HaplotypeCallerSpark -I hdfs:///user/$USER/small_spark_eval/out/bqsr-sharded -R hdfs:///user/$USER/small_spark_eval/human_g1k_v37.20.21.2bit -O hdfs:///user/$USER/small_spark_eval/out/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.vcf -pairHMM AVX_LOGLESS_CACHING" 8 1 4g 4g
