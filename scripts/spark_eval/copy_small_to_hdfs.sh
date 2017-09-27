#!/usr/bin/env bash

# Download all required data for small BAM and store in HDFS. Use this for non-GCS clusters.

TARGET_DIR=${1:-small_spark_eval}

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
gsutil cp gs://broad-spark-eval-test-data/small/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam - | hadoop fs -put - $TARGET_DIR/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam
gsutil cp gs://broad-spark-eval-test-data/small/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam.bai - | hadoop fs -put - $TARGET_DIR/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam.bai

# Download reference
gsutil cp gs://broad-spark-eval-test-data/small/human_g1k_v37.20.21.2bit - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.20.21.2bit
gsutil cp gs://broad-spark-eval-test-data/small/human_g1k_v37.20.21.dict - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.20.21.dict
gsutil cp gs://broad-spark-eval-test-data/small/human_g1k_v37.20.21.fasta.fai - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.20.21.fasta.fai
gsutil cp gs://broad-spark-eval-test-data/small/human_g1k_v37.20.21.fasta - | hadoop fs -put - $TARGET_DIR/human_g1k_v37.20.21.fasta

# Download known sites VCF
gsutil cp gs://broad-spark-eval-test-data/small/dbsnp_138.b37.20.21.vcf - | hadoop fs -put - $TARGET_DIR/dbsnp_138.b37.20.21.vcf

# List data
hadoop fs -ls -h $TARGET_DIR
