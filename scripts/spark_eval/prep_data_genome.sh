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
gsutil cp gs://hellbender/q4_spark_eval/WGS-G94982-NA12878.bam - | hadoop fs -put - $TARGET_DIR/WGS-G94982-NA12878.bam
gsutil cp gs://hellbender/q4_spark_eval/WGS-G94982-NA12878.bai - | hadoop fs -put - $TARGET_DIR/WGS-G94982-NA12878.bai

# Download reference (b37)
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.dict.gz | gunzip | hadoop fs -put - $TARGET_DIR/human_g1k_v37.dict
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.gz | gunzip | hadoop fs -put - $TARGET_DIR/human_g1k_v37.fasta
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37.fasta.fai.gz | gunzip | hadoop fs -put - $TARGET_DIR/human_g1k_v37.fasta.fai

# Generate 2bit
hadoop fs -get $TARGET_DIR/human_g1k_v37.fasta
curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x faToTwoBit
./faToTwoBit human_g1k_v37.fasta human_g1k_v37.2bit
hadoop fs -put human_g1k_v37.2bit $TARGET_DIR/human_g1k_v37.2bit
rm human_g1k_v37.* faToTwoBit

# Download known sites VCF (b37)
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz | gunzip | hadoop fs -put - $TARGET_DIR/dbsnp_138.b37.vcf

# Download exome intervals
gsutil cp gs://hellbender/q4_spark_eval/Broad.human.exome.b37.interval_list Broad.human.exome.b37.interval_list

# List data
hadoop fs -ls -h $TARGET_DIR
