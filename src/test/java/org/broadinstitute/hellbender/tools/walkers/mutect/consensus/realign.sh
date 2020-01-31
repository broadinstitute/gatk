#!/bin/bash

if [ "$#" -ne 3 ]; then
  echo "this script must 3 arguments: <grouped_bam> <base_name> <out_dir>"
  exit 1
fi

export consensus_bam=$1
export base_name=$2
export out_dir=$3

export bwa_path="/Users/tsato/Downloads/bwa-0.7.17/bwa"
export picard_jar="/Users/tsato/workspace/picard/build/libs/picard.jar"
export reference="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
export fastq1_name="${out_dir}${base_name}".1.fastq
export fastq2_name="${out_dir}${base_name}".2.fastq

java -jar ${picard_jar} SamToFastq \
INPUT="${consensus_bam}" \
FASTQ="${fastq1_name}" \
SECOND_END_FASTQ="${fastq2_name}" \
UNPAIRED_FASTQ="${out_dir}${base_name}".unpaired.fastq

# Align fastq with bwa
${bwa_path} mem -K 100000000 -t 16 \
${reference} \
"${fastq1_name}" \
"${fastq2_name}" > "${out_dir}${base_name}".aligned.sam

# Then MergeBamAlignment

