#!/usr/bin/env bash

# Run BQSR on genome data on a Spark cluster.

. utils.sh

# Small BAM (6.6GB), but full reference and known sites VCF
#time_gatk "BQSRPipelineSpark -I bam/CEUTrio.HiSeq.WGS.b37.ch1.1m-65m.NA12878.bam -O hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded-6gb-bam --shardedOutput true -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit --knownSites hdfs:///user/$USER/q4_spark_eval/dbsnp_138.b37.vcf --joinStrategy OVERLAPS_PARTITIONER" 4 8 32g 4g

time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/q4_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/q4_spark_eval/out/bqsr-sharded --shardedOutput true -R hdfs:///user/$USER/q4_spark_eval/human_g1k_v37.2bit --knownSites hdfs:///user/$USER/q4_spark_eval/dbsnp_138.b37.vcf --joinStrategy OVERLAPS_PARTITIONER" 4 8 48g 4g
