#!/bin/bash

set -eu

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "BEGIN TO PERFORM SIMPLE COMPARISON BETWEEN MASTER AND FEATURE"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'

if [[ "$GATK_VCF_MASTER" == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$GATK_VCF_MASTER")
    pattern=".gz"
    GATKVCF_m=${COMPRESSED//$pattern/}
    bgzip -c -d "$GATK_VCF_MASTER" > "$GATKVCF_m"
elif [[ "$GATK_VCF_MASTER" == *.vcf  ]]; then
    GATKVCF_m="$GATK_VCF_MASTER"
fi

if [[ "$GATK_VCF_FEATURE" == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$GATK_VCF_FEATURE")
    pattern=".gz"
    GATKVCF_f=${COMPRESSED//$pattern/}
    bgzip -c -d "$GATK_VCF_FEATURE" > "$GATKVCF_f"
elif [[ "$GATK_VCF_FEATURE" == *.vcf  ]]; then
    GATKVCF_f="$GATK_VCF_FEATURE"
fi

if [[ "$GATK_CPX_DERIVED_ONE_SEG_VCF_FEATURE" == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$GATK_CPX_DERIVED_ONE_SEG_VCF_FEATURE")
    pattern=".gz"
    GATKVCF_c_1=${COMPRESSED//$pattern/}
    bgzip -c -d "$GATK_CPX_DERIVED_ONE_SEG_VCF_FEATURE" > "$GATKVCF_c_1"
elif [[ "$GATK_CPX_DERIVED_ONE_SEG_VCF_FEATURE" == *.vcf  ]]; then
    GATKVCF_c_1="$GATK_CPX_DERIVED_ONE_SEG_VCF_FEATURE"
fi


if [[ "$GATK_CPX_DERIVED_MULTI_SEG_VCF_FEATURE" == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$GATK_CPX_DERIVED_MULTI_SEG_VCF_FEATURE")
    pattern=".gz"
    GATKVCF_c_2=${COMPRESSED//$pattern/}
    bgzip -c -d "$GATK_CPX_DERIVED_MULTI_SEG_VCF_FEATURE" > "$GATKVCF_c_2"
elif [[ "$GATK_CPX_DERIVED_MULTI_SEG_VCF_FEATURE" == *.vcf  ]]; then
    GATKVCF_c_2="$GATK_CPX_DERIVED_MULTI_SEG_VCF_FEATURE"
fi

cd "$ANALYSIS_DIR_MASTER_VS_FEATURE"

###################
## collect the IDs then do some accounting
cat "$ANALYSIS_DIR_MASTER""Deletion/gatkIDsWithMatchingPacBio.txt" \
    "$ANALYSIS_DIR_MASTER""InsDupRPL/gatkIDsWithMatchingPacBio.txt" \
    | sort | uniq > "master.gatkIDsWithMatchingPacBio.txt"

cat "$ANALYSIS_DIR_MASTER""Deletion/gatkIDsNoMatchingPacBio.txt" \
    "$ANALYSIS_DIR_MASTER""InsDupRPL/gatkIDsNoMatchingPacBio.txt" \
    | sort | uniq > "master.gatkIDsNoMatchingPacBio.txt"

cat "$ANALYSIS_DIR_FEATURE""Deletion/gatkIDsWithMatchingPacBio.txt" \
    "$ANALYSIS_DIR_FEATURE""InsDupRPL/gatkIDsWithMatchingPacBio.txt" \
    | sort | uniq > "feature.gatkIDsWithMatchingPacBio.txt"

cat "$ANALYSIS_DIR_FEATURE""Deletion/gatkIDsNoMatchingPacBio.txt" \
    "$ANALYSIS_DIR_FEATURE""InsDupRPL/gatkIDsNoMatchingPacBio.txt" \
    | sort | uniq > "feature.gatkIDsNoMatchingPacBio.txt"


echo "Number of variants \"validated\" by PacBio haploid callsets from master:"
wc -l "master.gatkIDsWithMatchingPacBio.txt" | awk '{print $1}'
echo "Number of variants \"validated\" by PacBio haploid callsets from feature:"
wc -l "feature.gatkIDsWithMatchingPacBio.txt" | awk '{print $1}'
echo "The number of variants shared by them is"
comm -12 \
    "master.gatkIDsWithMatchingPacBio.txt" \
    "feature.gatkIDsWithMatchingPacBio.txt" \
    > "shared.gatkIDsWithMatchingPacBio.txt"
wc -l "shared.gatkIDsWithMatchingPacBio.txt" | awk '{print $1}'


echo "Number of variants NOT \"validated\" by PacBio haploid callsets from master:"
wc -l "master.gatkIDsNoMatchingPacBio.txt" | awk '{print $1}'
echo "Number of variants NOT \"validated\" by PacBio haploid callsets from feature:"
wc -l "feature.gatkIDsNoMatchingPacBio.txt" | awk '{print $1}'
echo "The number of variants shared by them is"
comm -12 \
    "master.gatkIDsNoMatchingPacBio.txt" \
    "feature.gatkIDsNoMatchingPacBio.txt" \
    > "shared.gatkIDsNoMatchingPacBio.txt"
wc -l "shared.gatkIDsNoMatchingPacBio.txt" | awk '{print $1}'

###################
# now do some filtering on featuer callset: MQ_min >= 20
MQ_FT=20
LEN_FT=30
echo
echo "Now do some filtering on featuer callset: MQ_min >= $MQ_FT && ALN_min >= $LEN_FT"

cat "$GATKVCF_f" "$GATKVCF_c_1" "$GATKVCF_c_2" | grep -v '^#' > temp.feature.var.txt
awk '{print $8}' temp.feature.var.txt > temp.info.txt
grep -Eo "MAPPING_QUALITIES=[0-9]+(,[0-9]+){0,}" temp.info.txt | \
    awk -F '=' '{print $2}' | tr , ' ' | \
    awk -v filter="$MQ_FT" '{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i; if(m>=filter) print "PASS"; else print "FAIL"}' > \
    temp.mqFilter.txt
grep -Eo "MAX_ALIGN_LENGTH=[0-9]+" temp.info.txt | awk -F '=' '{print $2}' > temp.maxAlignLen.txt
awk -v filter="$LEN_FT" '{if ($1 >= filter) print "PASS"; else print "FAIL";}' temp.maxAlignLen.txt > \
    temp.alnLenFilter.txt

paste <(awk 'BEGIN {OFS="	"}; {print $1, $2, $3}' temp.feature.var.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.mqFilter.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.alnLenFilter.txt) \
      | grep -vF 'FAIL' | awk '{print $3}' | sort | uniq > temp.pass.gatkIDs.txt

comm -12 "feature.gatkIDsWithMatchingPacBio.txt" temp.pass.gatkIDs.txt > "feature.mq20_aln30.gatkIDsWithMatchingPacBio.txt"
echo "Number of variants \"validated\" by PacBio haploid callsets from feature after filtering with MQ $MQ_FT and contig alignment length $LEN_FT:"
wc -l "feature.mq20_aln30.gatkIDsWithMatchingPacBio.txt" | awk '{print $1}'

comm -12 "feature.gatkIDsNoMatchingPacBio.txt" temp.pass.gatkIDs.txt > "feature.mq20_aln30.gatkIDsNoMatchingPacBio.txt"
echo "Number of variants NOT \"validated\" by PacBio haploid callsets from feature after filtering with MQ $MQ_FT and contig alignment length $LEN_FT:"
wc -l "feature.mq20_aln30.gatkIDsNoMatchingPacBio.txt" | awk '{print $1}'

echo

###################
## extract assembly contigs that triggered only calls in either mater or feature but not both
echo "Now extract FN (validated in master but absent in feature) and FP (appear in feature but unvalidated by PacBio)"

# master and feature only contigs
parallel_grep 1 "master.gatkIDsWithMatchingPacBio.txt" "$GATKVCF_m" | \
    grep -Eo 'CTG_NAMES=((,)?asm[0-9]{6,6}:tig[0-9]{5,5})+' | awk -F '=' '{print $2}' | tr ',' '\n' | \
    sort | uniq > "master.ctgNamesWithMatchingPacBio.txt"
parallel_grep 1 "feature.gatkIDsWithMatchingPacBio.txt" temp.feature.var.txt | \
    grep -Eo 'CTG_NAMES=((,)?asm[0-9]{6,6}:tig[0-9]{5,5})+' | awk -F '=' '{print $2}' | tr ',' '\n' | \
    sort | uniq > "feature.ctgNamesWithMatchingPacBio.txt"
comm -23 "master.ctgNamesWithMatchingPacBio.txt" \
    "feature.ctgNamesWithMatchingPacBio.txt" \
    > "masterOnly.ctgNamesWithMatchingPacBio.txt"

parallel_grep 1 "master.gatkIDsNoMatchingPacBio.txt" "$GATKVCF_m" | \
    grep -Eo 'CTG_NAMES=((,)?asm[0-9]{6,6}:tig[0-9]{5,5})+' | awk -F '=' '{print $2}' | tr ',' '\n' | \
    sort | uniq > "master.ctgNamesNoMatchingPacBio.txt"
parallel_grep 1 "feature.gatkIDsNoMatchingPacBio.txt" temp.feature.var.txt | \
    grep -Eo 'CTG_NAMES=((,)?asm[0-9]{6,6}:tig[0-9]{5,5})+' | awk -F '=' '{print $2}' | tr ',' '\n' | \
    sort | uniq > "feature.ctgNamesNoMatchingPacBio.txt"
comm -13 "master.ctgNamesNoMatchingPacBio.txt" \
    "feature.ctgNamesNoMatchingPacBio.txt" \
    > "featureOnly.ctgNamesNoMatchingPacBio.txt"

## FN
parallel_grep 1 "masterOnly.ctgNamesWithMatchingPacBio.txt" "$GATKVCF_m" | \
    awk 'match($8, /END=[0-9]+/){ending=substr($8,RSTART+4,RLENGTH-4); print $1","$2","ending","$3}' \
    > temp.FN.bed
sed -e $'s/,/\t/g' temp.FN.bed | gsort -V | uniq > suspectsFN.bed

## FP
parallel_grep 1 "featureOnly.ctgNamesNOMatchingPacBio.txt" temp.feature.var.txt | \
    awk 'match($8, /END=[0-9]+/){ending=substr($8,RSTART+4,RLENGTH-4); print $1","$2","ending","$3;next;}{print $1","$2","$2","$3}' \
    > temp.FP.bed
sed -e $'s/,/\t/g' temp.FP.bed | gsort -V | uniq > suspectsFP.bed
echo "BED file storing records:"
echo "  1. matched by PacBio exist only in master, and "
echo "  1. not matched by PacBio exist only in feature"
echo "is output to"
realpath suspectsFN.bed
realpath suspectsFP.bed

echo "Among FN suspects"
grep -c -F 'DEL' suspectsFN.bed
echo "  are deletions, and "
grep -c -F 'INS' suspectsFN.bed
echo "  are insertions"
echo "Among FP suspects"
grep -c -F 'DEL' suspectsFP.bed
echo "  are deletions, and "
grep -c -F 'INS' suspectsFP.bed
echo "  are insertions"

########## clean up
rm -f temp*
if [[ "$GATK_VCF_MASTER" == *.vcf.gz ]]; then
    rm -f "$GATKVCF_m"
fi
if [[ "$GATK_VCF_FEATURE" == *.vcf.gz ]]; then
    rm -f "$GATKVCF_f"
fi
if [[ "$GATK_CPX_DERIVED_ONE_SEG_VCF_FEATURE" == *.vcf.gz ]]; then
    rm -f "$GATKVCF_c_1"
fi
if [[ "$GATK_CPX_DERIVED_MULTI_SEG_VCF_FEATURE" == *.vcf.gz ]]; then
    rm -f "$GATKVCF_c_2"
fi

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE SIMPLE COMPARISON BETWEEN MASTER AND FEATURE"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
