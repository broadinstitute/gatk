#!/bin/bash

set -eu

mkdir -p ../postprocess

cat ../Master/GATKVCFHeader.txt \
	../Master/GATK_primaryContigs_var_no_warn_Duplicates.txt \
	> ../postprocess/GATK_primaryContigs_var_no_warn_Duplicates.MASTER.vcf

grep -v 'BND' \
	../Feature/GATK_primaryContigs_var_no_warn_Duplicates.txt \
	> temp.GATK_primaryContigs_var_no_warn_Duplicates.bdnLess.txt
cat ../Feature/GATKVCFHeader.txt \
	temp.GATK_primaryContigs_var_no_warn_Duplicates.bdnLess.txt \
	> ../postprocess/GATK_primaryContigs_var_no_warn_Duplicates.FEATURE.bndLess.vcf

awk 'match($8, /END=[0-9]+/){endpos=substr($8,RSTART+4,RLENGTH-4); print $1","$2","endpos}' \
	../Master/GATK_primaryContigs_var_no_warn_Duplicates.txt \
	> temp.bed
awk 'match($8, /END=[0-9]+/){endpos=substr($8,RSTART+4,RLENGTH-4); print $1","$2","endpos}' \
	temp.GATK_primaryContigs_var_no_warn_Duplicates.bdnLess.txt \
	>> temp.bed

cat temp.bed | tr ',' '\t' | gsort -V > ../postprocess/duplicateRecords.bed

rm -f temp*
