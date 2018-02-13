#!/bin/bash

set -eu

if [[ "$#" -lt 9 ]]; then
    echo "Please provide:"
    echo "[1] absolute path to GATK insertion VCF"
    echo "[2] absolute path to GATK duplication VCF"
    echo "[3] absolute path to Manta insertion VCF"
    echo "[4] absolute path to Manta duplication VCF"
    echo "[5] absolute path to the directory where analysis outputs are to be written to"
    echo "[6] absolute path to GATK active region interval bed file"
    echo "[7] validataion call set tool name (\"Manta\"|\"Pacbio\")"
    echo "[8] absolute path to file listing uniq GATK IDs with overlapping Manta calls, and"
    echo "[9] absolute path to file listing uniq Manta IDs with overlapping GATK calls, or"
    echo "[8] absolute path to file listing uniq PacBio IDs with overlapping GATK calls, and"
    echo "[9] absolute path to file listing uniq PacBio IDs with overlapping Manta calls, and"
    echo "[10] absolute path to file listing uniq GATK IDs with overlapping Manta calls, and"
    echo "[11] absolute path to file listing uniq GATK IDs with overlapping PacBio calls"

    exit 1
fi

GATK_INS_VCF=$1
GATK_DUP_VCF=$2

MANTA_INS_VCF=$3
MANTA_DUP_VCF=$4

OUTPUT_DIR_ROOT=$5
if [[ ! $5 =~ .+/$ ]]; then
    OUTPUT_DIR_ROOT+="/"
fi
mkdir -p "$OUTPUT_DIR_ROOT""Results"
OUTPUT_DIR="$OUTPUT_DIR_ROOT""Results/"

GATK_ACTIVE_REGIONS=$6

echo "#################################################"
echo "Call sets sensitivity and specificity checks"
echo "#################################################"
echo

########## get bed file of Manta records without matching GATK records

if [[ $7 == "Manta" ]]; then

    ## get gatk ids un-validated by Manta ("false positive")
    grep -v '^#' "$GATK_INS_VCF" > temp.txt
    grep -v '^#' "$GATK_DUP_VCF" >> temp.txt
    grep -vf "$8" temp.txt | \
        awk 'BEGIN {OFS="	"}; match($8, /END=[0-9]+/){print $1, $2, substr($8,RSTART+4,RLENGTH-4), $3}' > \
        temp.gatk_notMatched.bed
    gsort -V -k1,1 -k2,3n temp.gatk_notMatched.bed | \
        uniq > "$OUTPUT_DIR""gatk_notMatched.bed"
    echo "Number of GATK records without matching Manta calls:"
    wc -l "$OUTPUT_DIR""gatk_notMatched.bed" | awk '{print $1}'
    echo `grep -F 'INS' "$OUTPUT_DIR""gatk_notMatched.bed" | wc -l` "are called as insertions"
    echo `grep -F 'DUP' "$OUTPUT_DIR""gatk_notMatched.bed" | wc -l` "are called as duplication"
    echo

    ## get manta ids not caught by GATK ("false negative")
    grep -v '^#' "$MANTA_INS_VCF" > temp.txt
    grep -v '^#' "$MANTA_DUP_VCF" >> temp.txt
    grep -vf "$9" temp.txt | \
        awk 'BEGIN {OFS="	"}; match($8, /END=[0-9]+/){print $1, $2, substr($8,RSTART+4,RLENGTH-4), $3}' > \
        temp.manta_notMatched.bed
    gsort -V -k1,1 -k2,3n temp.manta_notMatched.bed | \
        uniq > "$OUTPUT_DIR""manta_notMatched.bed"
    echo "Number of Manta records without matching GATK calls:"
    wc -l "$OUTPUT_DIR""manta_notMatched.bed" | awk '{print $1}'
    echo `grep -F 'INS' "$OUTPUT_DIR""manta_notMatched.bed" | wc -l` "are called as insertions"
    echo `grep -F 'DUP' "$OUTPUT_DIR""manta_notMatched.bed" | wc -l` "are called as duplication"
    echo

    ## active region sensitivity
    echo "Overlapping the unmatched Manta calls with GATK-SV pipeline local assembly intervals"
    bedtools window \
        -a "$OUTPUT_DIR""manta_notMatched.bed" \
        -b "$GATK_ACTIVE_REGIONS" \
        -w 50 | 
        awk '{print $NF}' | \
        sort | uniq | wc -l
    echo "If this number is close to the number of unmatched Manta calls, it is pointing to a possibility"
    echo " that the dicovery stage in the GATK-SV pipeline is having an overly agressive filtering step."
    echo

    rm -f temp*
elif [[ $7 == "Pacbio" ]]; then
    
    ## validated Manta calls not caught by GATK
    echo "Number of Manta calls validated by PacBio but not caught in the GATK call set:"
    comm -13 "$8" "$9" |
        awk -F ':' 'BEGIN{OFS="	"} {print $1, $2-1, $2}' > \
        "$OUTPUT_DIR""mantaRecordsValidatedByPacbioUncaughtByGATK.bed"
    wc -l "$OUTPUT_DIR""mantaRecordsValidatedByPacbioUncaughtByGATK.bed" | awk '{print $1}'
    echo 

    ## active region sensitivity
    echo "Overlapping the unmatched, Pacbio-validated Manta calls with GATK-SV pipeline local assembly intervals"
    bedtools window \
        -a "$OUTPUT_DIR""mantaRecordsValidatedByPacbioUncaughtByGATK.bed" \
        -b "$GATK_ACTIVE_REGIONS" \
        -w 50 |
        awk '{print $NF}' | \
        sort | uniq | wc -l | awk '{print $1}'
    echo "If this number is close to the number of unmatched, validated Manta calls, it is pointing to a possibility"
    echo " that the dicovery stage in the GATK-SV pipeline is having an overly agressive filtering step."
    echo

    ## Manta and GATK calls un-validated by PacBio
    echo "Number of calls shared by GATK and Manta but un-validated by PacBio:"
    comm -23 "${10}" "${11}" > \
        "$OUTPUT_DIR""recordsSharedByGATKandManta.unvalidatedByPacbio.gatkIDs.txt"
    wc -l "$OUTPUT_DIR""recordsSharedByGATKandManta.unvalidatedByPacbio.gatkIDs.txt" | awk '{print $1}'
    echo

else
    echo -e "Un-recogonazied validate call set tool name:" "$7"
    exit 1
fi