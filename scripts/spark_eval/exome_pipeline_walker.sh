#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on exome data using walkers.

. utils.sh

WORKDIR=$HOME

if [ ! -f $WORKDIR/dbsnp_138.hg18.vcf ]; then
  gsutil cp gs://gatk-tom-testdata-exome/dbsnp_138.hg18.vcf $WORKDIR/dbsnp_138.hg18.vcf
fi

if [ ! -f $WORKDIR/dbsnp_138.hg18.vcf.idx ]; then
  curl ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg18/dbsnp_138.hg18.vcf.idx.gz | gunzip > $WORKDIR/dbsnp_138.hg18.vcf.idx
fi

if [ ! -f $WORKDIR/Homo_sapiens_assembly18.dict ]; then
  gsutil cp gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.dict $WORKDIR/Homo_sapiens_assembly18.dict
fi

if [ ! -f $WORKDIR/Homo_sapiens_assembly18.fasta ]; then
  gsutil cp gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.fasta $WORKDIR/Homo_sapiens_assembly18.fasta
fi

if [ ! -f $WORKDIR/Homo_sapiens_assembly18.fasta.fai ]; then
  gsutil cp gs://gatk-tom-testdata-exome/Homo_sapiens_assembly18.fasta.fai $WORKDIR/Homo_sapiens_assembly18.fasta.fai
fi

if [ ! -f $WORKDIR/NA12878.ga2.exome.maq.raw.bam ]; then
  gsutil cp gs://gatk-tom-testdata-exome/NA12878.ga2.exome.maq.raw.bam $WORKDIR/NA12878.ga2.exome.maq.raw.bam
fi

time_gatk_walker "MarkDuplicates -I $WORKDIR/NA12878.ga2.exome.maq.raw.bam -O $WORKDIR/NA12878.ga2.exome.maq.raw.md.bam --METRICS_FILE $WORKDIR/NA12878.ga2.exome.maq.raw.metrics.txt --VALIDATION_STRINGENCY LENIENT"
time_gatk_walker "BaseRecalibrator -I $WORKDIR/NA12878.ga2.exome.maq.raw.md.bam -O $WORKDIR/NA12878.ga2.exome.maq.raw.recal_data.table -R $WORKDIR/Homo_sapiens_assembly18.fasta --knownSites $WORKDIR/dbsnp_138.hg18.vcf"
time_gatk_walker "ApplyBQSR -I $WORKDIR/NA12878.ga2.exome.maq.raw.md.bam -O $WORKDIR/NA12878.ga2.exome.maq.raw.bqsr.bam -bqsr $WORKDIR/NA12878.ga2.exome.maq.raw.recal_data.table"
time_gatk_walker "HaplotypeCaller -I $WORKDIR/NA12878.ga2.exome.maq.raw.bqsr.bam -R $WORKDIR/Homo_sapiens_assembly18.fasta -O $WORKDIR/NA12878.ga2.exome.maq.raw.vcf -pairHMM AVX_LOGLESS_CACHING -maxReadsPerAlignmentStart 10"