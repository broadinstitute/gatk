#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on genome data in GCS.

. utils.sh

time_gatk "ReadsPipelineSpark -I gs://broad-spark-eval-test-data/genome/WGS-G94982-NA12878-no-NC_007605.bam -O gs://broad-spark-eval-test-data/genome/out/WGS-G94982-NA12878.vcf -R gs://broad-spark-eval-test-data/genome/human_g1k_v37.fasta --known-sites gs://broad-spark-eval-test-data/genome/dbsnp_138.b37.vcf.gz -pairHMM AVX_LOGLESS_CACHING --maxReadsPerAlignmentStart 10" 40 7 20g 8g
