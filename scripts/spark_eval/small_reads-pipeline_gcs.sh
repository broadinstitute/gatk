#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on small data in GCS.

. utils.sh

time_gatk "ReadsPipelineSpark -I gs://broad-spark-eval-test-data/small/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam -O gs://broad-spark-eval-test-data/small/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.vcf -R gs://broad-spark-eval-test-data/small/human_g1k_v37.20.21.2bit --known-sites gs://broad-spark-eval-test-data/small/dbsnp_138.b37.20.21.vcf -pairHMM AVX_LOGLESS_CACHING --max-reads-per-alignment-start 10" 1 8 32g 4g
