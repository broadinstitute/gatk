#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on small data on a GCS Dataproc cluster. Data is in GCS.

# TODO: change output to a GCS bucket when writing works again

. utils.sh

time_gatk "ReadsPipelineSpark -I gs://hellbender/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam -O hdfs://${GCS_CLUSTER}-m:8020/user/$USER/small_spark_eval/out/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.vcf -R gs://hellbender/test/resources/large/human_g1k_v37.20.21.2bit --knownSites gs://hellbender/test/resources/large/dbsnp_138.b37.20.21.vcf -pairHMM AVX_LOGLESS_CACHING -maxReadsPerAlignmentStart 10" 1 8 32g 4g
