#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on exome data in GCS.

. utils.sh

time_gatk "ReadsPipelineSpark -I gs://broad-spark-eval-test-data/exome/NA12878.ga2.exome.maq.raw.bam -O gs://broad-spark-eval-test-data/exome/NA12878.ga2.exome.maq.raw.vcf -R gs://broad-spark-eval-test-data/exome/Homo_sapiens_assembly18.2bit --known-sites gs://broad-spark-eval-test-data/exome/dbsnp_138.hg18.vcf -pairHMM AVX_LOGLESS_CACHING --max-reads-per-alignment-start 10" 8 8 32g 4g