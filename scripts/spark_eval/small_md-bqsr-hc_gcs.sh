#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on small data in GCS.

. utils.sh

time_gatk "MarkDuplicatesSpark -I gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam -O hdfs:///user/$USER/small_spark_eval/out/markdups-sharded --sharded-output true" 8 1 4g 4g
time_gatk "BQSRPipelineSpark -I hdfs:///user/$USER/small_spark_eval/out/markdups-sharded -O hdfs:///user/$USER/small_spark_eval/out/bqsr-sharded --sharded-output true -R gs://hellbender/test/resources/large/human_g1k_v37.20.21.2bit --known-sites gs://hellbender/test/resources/large/dbsnp_138.b37.20.21.vcf -L 20:10000000-10100000" 1 8 32g 4g
time_gatk "HaplotypeCallerSpark -I hdfs:///user/$USER/small_spark_eval/out/bqsr-sharded -R gs://hellbender/test/resources/large/human_g1k_v37.20.21.2bit -O gs://broad-spark-eval-test-data/small/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.vcf -pairHMM AVX_LOGLESS_CACHING" 8 1 4g 4g
