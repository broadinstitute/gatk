#!/usr/bin/env bash

# Download all required data for genomes and store in HDFS. Use this for non-GCS clusters.

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
gsutil gs://broad-spark-eval-test-data/genome/WGS-G94982-NA12878-no-NC_007605.bam - | hadoop fs -put - $TARGET_DIR/WGS-G94982-NA12878-no-NC_007605.bam

# Download reference (b37) (ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/)
gsutil gs://broad-spark-eval-test-data/genome/human_g1k_v37.2bit - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.2bit
gsutil gs://broad-spark-eval-test-data/genome/human_g1k_v37.dict - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.dict
gsutil gs://broad-spark-eval-test-data/genome/human_g1k_v37.fasta - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.fasta
gsutil gs://broad-spark-eval-test-data/genome/human_g1k_v37.fasta.fai - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.fasta.fai

# (Code for generating 2bit)
#hadoop fs -get $TARGET_DIR/human_g1k_v37.fasta
#curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
#chmod +x faToTwoBit
#./faToTwoBit human_g1k_v37.fasta human_g1k_v37.2bit

# Download known sites VCF (b37) (ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/)
gsutil gs://broad-spark-eval-test-data/genome/dbsnp_138.b37.vcf - | hadoop fs -put - $TARGET_DIR/dbsnp_138.b37.vcf

# Download exome intervals
gsutil cp gs://broad-spark-eval-test-data/genome/Broad.human.exome.b37.interval_list Broad.human.exome.b37.interval_list

# List data
hadoop fs -ls -h $TARGET_DIR
