#!/usr/bin/env bash

# Download all required data for genomes and store in HDFS.

TARGET_DIR=${1:-q4_spark_eval}

hadoop fs -stat $TARGET_DIR > /dev/null 2>&1
if [ $? -eq 0 ]; then
  echo "$TARGET_DIR already exists. Delete it and try again."
  exit 1
fi

set -e
set -x

# Create data directory in HDFS
hadoop fs -mkdir -p $TARGET_DIR

# Download WGS BAM
#gsutil cp gs://hellbender/q4_spark_eval/WGS-G94982-NA12878.bam - | hadoop fs -put - $TARGET_DIR/WGS-G94982-NA12878.bam
#gsutil cp gs://hellbender/q4_spark_eval/WGS-G94982-NA12878.bai - | hadoop fs -put - $TARGET_DIR/WGS-G94982-NA12878.bai
# BAM with NC_007605 reads removed since this contig is not in the reference
gsutil cp gs://gatk-tom-testdata/WGS-G94982-NA12878-no-NC_007605.bam - | hadoop fs -put - /user/tom/q4_spark_eval/WGS-G94982-NA12878-no-NC_007605.bam

# Download reference (b37)
gsutil cp gs://gatk-tom-testdata/human_g1k_v37.2bit - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.2bit
gsutil cp gs://gatk-tom-testdata/human_g1k_v37.dict - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.dict
gsutil cp gs://gatk-tom-testdata/human_g1k_v37.fasta - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.fasta
gsutil cp gs://gatk-tom-testdata/human_g1k_v37.fasta.fai - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.fasta.fai

# Download known sites VCF (b37)
gsutil cp gs://gatk-tom-testdata/dbsnp_138.b37.vcf - | hadoop fs -put - $TARGET_DIR/dbsnp_138.b37.vcf

# Download exome intervals
gsutil cp gs://hellbender/q4_spark_eval/Broad.human.exome.b37.interval_list Broad.human.exome.b37.interval_list

# List data
hadoop fs -ls -h $TARGET_DIR
