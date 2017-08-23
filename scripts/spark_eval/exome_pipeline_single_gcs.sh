#!/usr/bin/env bash

# Run the pipeline (ReadsPipelineSpark) on exome data on a GCS Dataproc cluster. Data is in GCS.

. utils.sh

time_gatk "ReadsPipelineSpark -I gs://gatk-tom-testdata-exome/NA12878.ga2.exome.maq.raw.bam -O gs://gatk-tom-testdata-exome/NA12878.ga2.exome.maq.raw.vcf -R gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.2bit --knownSites gs://gatk-tom-testdata-exome/dbsnp_138.hg18.vcf -pairHMM AVX_LOGLESS_CACHING -maxReadsPerAlignmentStart 10" 4 8 32g 4g