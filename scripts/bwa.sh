#!/bin/bash

# Run bwa on YARN, loading files from HDFS and writing output to HDFS

# Sample usage:
#
# curl -s http://hypervolu.me/%7Eerik/genomes/E.coli_K12_MG1655.fa > E.coli_K12_MG1655.fa
# fastq-dump --split-files SRR1770413
# hadoop fs -mkdir bwa-test
# hadoop fs -put E.coli_K12_MG1655.fa bwa-test/E.coli_K12_MG1655.fa
# hadoop fs -put SRR1770413_1.fastq bwa-test/SRR1770413_1.fastq
# hadoop fs -put SRR1770413_2.fastq bwa-test/SRR1770413_2.fastq
# hadoop fs -chmod 777 bwa-test
# hadoop fs -rm /user/$USER/bwa-test/k12.sam
# yarn jar /opt/cloudera/parcels/CDH-5.6.0-1.cdh5.6.0.p0.45/lib/hadoop-yarn/hadoop-yarn-applications-distributedshell.jar \
# --jar /opt/cloudera/parcels/CDH-5.6.0-1.cdh5.6.0.p0.45/lib/hadoop-yarn/hadoop-yarn-applications-distributedshell.jar \
# --num_containers 1 --master_memory 512 --master_vcores 2 --container_memory 2048 --container_vcores 1 \
# --shell_script bwa.sh \
# --shell_args "/user/$USER/bwa-test/E.coli_K12_MG1655.fa /user/$USER/bwa-test/SRR1770413_1.fastq /user/$USER/bwa-test/SRR1770413_2.fastq /user/$USER/bwa-test/k12.sam"

set -x
pwd

if [ $# -eq 0 ]; then
  echo "Usage: $0 [ref.fa] [seq1.fq] [seq2.fq] [out.sam]";
  exit 1
fi

ref=$1
seq1=$2
seq2=$2
out=$4

# Install BWA
if [ ! -e /tmp/bwa ]; then
  cd /tmp
  git clone https://github.com/lh3/bwa
  (cd bwa && make)
  cd -
fi
BWA=/tmp/bwa/bwa

# Localize files from HDFS
localref=$(basename $ref)
localseq1=$(basename $seq1)
localseq2=$(basename $seq2)
hadoop fs -get $ref $localref
hadoop fs -get $seq1 $localseq1
hadoop fs -get $seq2 $localseq2

# Run BWA
$BWA index $localref
$BWA mem -t 1 -R '@RG\tID:K12\tSM:K12' \
  $localref $localseq1 $localseq2 | hadoop fs -put - $out