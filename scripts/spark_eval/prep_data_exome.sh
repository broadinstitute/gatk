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
curl ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20101201_cg_NA12878/NA12878.ga2.exome.maq.raw.bam | hadoop fs -put - $TARGET_DIR/NA12878.ga2.exome.maq.raw.bam

# Download reference (hg18)
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg18/Homo_sapiens_assembly18.dict.gz | gunzip | hadoop fs -put - $TARGET_DIR/Homo_sapiens_assembly18.dict
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg18/Homo_sapiens_assembly18.fasta.fai.gz | gunzip | hadoop fs -put - $TARGET_DIR/Homo_sapiens_assembly18.fasta.fai
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg18/Homo_sapiens_assembly18.fasta.gz | gunzip | hadoop fs -put - $TARGET_DIR/Homo_sapiens_assembly18.fasta

# Generate 2bit
hadoop fs -get $TARGET_DIR/Homo_sapiens_assembly18.fasta
curl -O http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
chmod +x faToTwoBit
./faToTwoBit Homo_sapiens_assembly18.fasta Homo_sapiens_assembly18.2bit
hadoop fs -put Homo_sapiens_assembly18.2bit $TARGET_DIR/Homo_sapiens_assembly18.2bit
rm Homo_sapiens_assembly18.* faToTwoBit

# Download known sites VCF (hg18)
curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg18/dbsnp_138.hg18.vcf.gz | gunzip | hadoop fs -put - $TARGET_DIR/dbsnp_138.hg18.vcf

# List data
hadoop fs -ls -h $TARGET_DIR
