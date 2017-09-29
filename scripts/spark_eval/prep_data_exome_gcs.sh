#!/usr/bin/env bash

# Download all required data for exomes and store in HDFS.

TARGET_DIR=${1:-exome_spark_eval}

hadoop fs -stat $TARGET_DIR > /dev/null 2>&1
if [ $? -eq 0 ]; then
  echo "$TARGET_DIR already exists. Delete it and try again."
  exit 1
fi

set -e
set -x

# Create data directory in HDFS
hadoop fs -mkdir -p $TARGET_DIR

# Download exome BAM
gsutil cp gs://gatk-tom-testdata-exome/NA12878.ga2.exome.maq.raw.bam - | hadoop fs -put - $TARGET_DIR/NA12878.ga2.exome.maq.raw.bam

# Download reference (hg18)
gsutil cp gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.2bit - | hadoop fs -put - $TARGET_DIR/Homo_sapiens_assembly18.2bit
gsutil cp gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.dict - | hadoop fs -put - $TARGET_DIR/Homo_sapiens_assembly18.dict
gsutil cp gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.fasta.fai - | hadoop fs -put - $TARGET_DIR/Homo_sapiens_assembly18.fasta.fai
gsutil cp gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.fasta - | hadoop fs -put - $TARGET_DIR/Homo_sapiens_assembly18.fasta

# Download known sites VCF (hg18)
gsutil cp gs://gatk-tom-testdata-exome/dbsnp_138.hg18.vcf - | hadoop fs -put - $TARGET_DIR/dbsnp_138.hg18.vcf

# List data
hadoop fs -ls -h $TARGET_DIR
