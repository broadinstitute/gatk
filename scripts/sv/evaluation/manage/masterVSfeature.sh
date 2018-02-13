#!/bin/bash

set -eu

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "BEGIN TO PERFORM SIMPLE COMPARISON BETWEEN MASTER AND FEATURE"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'

echo "Diverts different types of variants in VCF files"
echo "  and collects & saves variants sizes by type"


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

if [[ "$GATK_CPX_VCF_FEATURE" == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$GATK_CPX_VCF_FEATURE")
    pattern=".gz"
    GATKVCF_c=${COMPRESSED//$pattern/}
    bgzip -c -d "$GATK_CPX_VCF_FEATURE" > "$GATKVCF_c"
elif [[ "$GATK_CPX_VCF_FEATURE" == *.vcf  ]]; then
    GATKVCF_c="$GATK_CPX_VCF_FEATURE"
fi

mkdir -p "../masterVSfeature" && cd "../masterVSfeature/"

###################
## collect the IDs then do some accounting
cat "$ANALYSIS_DIR_MASTER""Deletion/gatkIDsWithMatchingPacBio.txt" \
    "$ANALYSIS_DIR_MASTER""InsDupRPL/gatkIDsWithMatchingPacBio.txt" \
    | sort > "master.gatkIDsWithMatchingPacBio.txt"

cat "$ANALYSIS_DIR_MASTER""Deletion/gatkIDsNoMatchingPacBio.txt" \
    "$ANALYSIS_DIR_MASTER""InsDupRPL/gatkIDsNoMatchingPacBio.txt" \
    | sort > "master.gatkIDsNoMatchingPacBio.txt"

cat "$ANALYSIS_DIR_FEATURE""Deletion/gatkIDsWithMatchingPacBio.txt" \
    "$ANALYSIS_DIR_FEATURE""InsDupRPL/gatkIDsWithMatchingPacBio.txt" \
    "$ANALYSIS_DIR_FEATURE""CPX/gatkIDsWithMatchingPacBio.txt" \
    | sort > "feature.gatkIDsWithMatchingPacBio.txt"

cat "$ANALYSIS_DIR_FEATURE""Deletion/gatkIDsNoMatchingPacBio.txt" \
    "$ANALYSIS_DIR_FEATURE""InsDupRPL/gatkIDsNoMatchingPacBio.txt" \
    "$ANALYSIS_DIR_FEATURE""CPX/gatkIDsNoMatchingPacBio.txt" \
    | sort > "feature.gatkIDsNoMatchingPacBio.txt"

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
wc -l "master.gatkIDsNOMatchingPacBio.txt" | awk '{print $1}'
echo "Number of variants NOT \"validated\" by PacBio haploid callsets from feature:"
wc -l "feature.gatkIDsNOMatchingPacBio.txt" | awk '{print $1}'
echo "The number of variants shared by them is"
comm -12 \
    "master.gatkIDsNOMatchingPacBio.txt" \
    "feature.gatkIDsNOMatchingPacBio.txt" \
    > "shared.gatkIDsNOMatchingPacBio.txt"
wc -l "shared.gatkIDsNoMatchingPacBio.txt" | awk '{print $1}'


###################
## extract assembly contigs that triggered only calls in either mater or feature but not both
comm -23 \
    "master.gatkIDsWithMatchingPacBio.txt" \
    "feature.gatkIDsWithMatchingPacBio.txt" \
    | sort > "temp.masterOnly.gatkIDsWithMatchingPacBio.txt"
grep -f "temp.masterOnly.gatkIDsWithMatchingPacBio.txt" "$GATKVCF_m" | \
    grep -Eo 'CTG_NAMES=asm[0-9]{6,6}:tig[0-9]{5,5}(,asm[0-9]{6,6}:tig[0-9]{5,5}){0,}' | \
    tr ',' '\n' | sort | uniq > masterOnly.ctgNamesWithMatchingPacBio.txt

comm -13 \
    "master.gatkIDsWithMatchingPacBio.txt" \
    "feature.gatkIDsWithMatchingPacBio.txt" \
    | sort > "temp.featureOnly.gatkIDsWithMatchingPacBio.txt"
grep -f "temp.featureOnly.gatkIDsWithMatchingPacBio.txt" "$GATKVCF_f" | \
    grep -Eo 'CTG_NAMES=asm[0-9]{6,6}:tig[0-9]{5,5}(,asm[0-9]{6,6}:tig[0-9]{5,5}){0,}' | \
    tr ',' '\n' | sort | uniq > temp.featureOnly.ctgNamesWithMatchingPacBio.1.txt
grep -f "temp.featureOnly.gatkIDsWithMatchingPacBio.txt" "$GATKVCF_c" | \
    grep -Eo 'CTG_NAMES=asm[0-9]{6,6}:tig[0-9]{5,5}(,asm[0-9]{6,6}:tig[0-9]{5,5}){0,}' | \
    tr ',' '\n' | sort | uniq > temp.featureOnly.ctgNamesWithMatchingPacBio.2.txt
cat temp.featureOnly.ctgNamesWithMatchingPacBio.1.txt \
    temp.featureOnly.ctgNamesWithMatchingPacBio.2.txt | \
    sort | uniq > featureOnly.ctgNamesWithMatchingPacBio.txt

## then extract BED files based on the assembly names
grep -f masterOnly.ctgNamesWithMatchingPacBio.txt "$GATKVCF_m" | \
    awk 'match($8, /END=[0-9]+/){ending=substr($8,RSTART+4,RLENGTH-4); print $1",$2,"ending",MASTER"}' \
    > temp.singletons.bed
grep -f featureOnly.ctgNamesWithMatchingPacBio.txt "$GATKVCF_f" | \
    awk 'match($8, /END=[0-9]+/){ending=substr($8,RSTART+4,RLENGTH-4); print $1",$2,"ending",FEATURE"}' \
    >> temp.singletons.bed
grep -f featureOnly.ctgNamesWithMatchingPacBio.txt "$GATKVCF_c" | \
    awk 'match($8, /END=[0-9]+/){ending=substr($8,RSTART+4,RLENGTH-4); print $1",$2,"ending",FEATURE"}' \
    >> temp.singletons.bed
cat temp.singletons.bed | tr ',' '\t' | gsort -V > singletons.bed

########## clean up
rm -f temp*
if [[ $1 == *.vcf.gz ]]; then
    rm -f "$GATKVCF_m"
fi
if [[ $2 == *.vcf.gz ]]; then
    rm -f "$GATKVCF_f"
fi
if [[ $3 == *.vcf.gz ]]; then
    rm -f "$GATKVCF_c"
fi

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE SIMPLE COMPARISON BETWEEN MASTER AND FEATURE"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'