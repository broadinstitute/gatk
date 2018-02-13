#!/bin/bash

set -eu

if [[ "$#" -ne 4 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to GATK VCF"
	echo "[2] absolute path to PacBio CHM1 VCF"
	echo "[3] absolute path to PacBio CHM13 VCF"
	echo "[4] absolute path to the directory where analysis outputs are to be written to"
	exit 1
fi

INPUTVCF=$1
if [[ -z $INPUTVCF ]]; then
	echo "Seems like there's no CPX vcf avaible. Quit."
	exit 1
fi
if [[ $1 == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$1")
    pattern=".gz"
    GATKVCF=${COMPRESSED//$pattern/}
    bgzip -c -d "$1" > "$GATKVCF"
elif [[ $1 == *.vcf  ]]; then
    GATKVCF=$1
fi

PACBIO_VCF_CHM1=$2
PACBIO_VCF_CHM13=$3

OUTPUT_DIR=$4
if [[ ! $4 =~ .+/$ ]]; then
	OUTPUT_DIR+="/"
fi

# MQ_FT=20
# LEN_FT=30
# if [[ "$#" -eq 4 ]]; then
#   MQ_FT=$3;
#   LEN_FT=$4;
# fi

echo "#################################################"
echo "CPX VCF parsing on:"
echo "  $GATKVCF"
echo "#################################################"
echo

VARIANTS_TXT="temp.GATK_primaryContigs_cpx.txt"
grep -v '^#'  "$GATKVCF" > "$VARIANTS_TXT"

awk '{print $8}' "$VARIANTS_TXT" > temp.info.txt

grep -Eo "END=[0-9]+" temp.info.txt | awk -F '=' '{print $2}' > temp.end.txt

# grep -Eo "MAPPING_QUALITIES=[0-9]+(,[0-9]+){0,}" temp.info.txt | \
#     awk -F '=' '{print $2}' | tr , ' ' | \
#     awk -v filter="$MQ_FT" '{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i; if(m>=filter) print "PASS"; else print "FAIL"}' > \
#     temp.mqFilter.txt

# grep -Eo "MAX_ALIGN_LENGTH=[0-9]+" temp.info.txt || true | awk -F '=' '{print $2}' > temp.maxAlignLen.txt
# awk -v filter="$LEN_FT" '{if ($1 >= filter) print "PASS"; else print "FAIL";}' temp.maxAlignLen.txt > \
#     temp.alnLenFilter.txt

paste <(awk 'BEGIN {OFS="	"}; {print $1, $2}' "$VARIANTS_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.end.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $3}' "$VARIANTS_TXT") \
      > "$OUTPUT_DIR""gatk.cpx.all.bed"
      # <(awk 'BEGIN {OFS="	"}; {print $1}' temp.maxAlignLen.txt) \
      # <(awk 'BEGIN {OFS=" "}; {print $1}' temp.mqFilter.txt) \
      # <(awk 'BEGIN {OFS=" "}; {print $1}' temp.alnLenFilter.txt) \
      # > temp.bed

# awk 'BEGIN {OFS="	"}; {s=$4";"$5";"$6";"$7";"$8; print $1, $2, $3, s}' temp.bed \
#     > "$OUTPUT_DIR""gatk.cleanDel.all.bed"

# grep -vF 'FAIL' "$OUTPUT_DIR""gatk.cleanDel.all.bed" \
#     > "$OUTPUT_DIR""gatk.cleanDel.bed"


bedtools intersect \
	-a "$OUTPUT_DIR""gatk.cpx.all.bed" \
	-b "$PACBIO_VCF_CHM1" \
	-c | awk '{if ($NF!=0) print $4}' | \
	sort | uniq > temp.match1.txt

bedtools intersect \
	-a "$OUTPUT_DIR""gatk.cpx.all.bed" \
	-b "$PACBIO_VCF_CHM13" \
	-c | awk '{if ($NF!=0) print $4}' | \
	sort | uniq > temp.match13.txt

cat temp.match1.txt temp.match13.txt | \
	sort | uniq > "$OUTPUT_DIR""gatkIDsWithMatchingPacBio.txt"


awk '{print $4}' "$OUTPUT_DIR""gatk.cpx.all.bed" | 
	sort | uniq > temp.gatkIDs.txt

comm -23 temp.gatkIDs.txt "$OUTPUT_DIR""gatkIDsWithMatchingPacBio.txt" \
	> "$OUTPUT_DIR""gatkIDsNoMatchingPacBio.txt"


echo "Total number of variants from GATK to be tested:"
wc -l "$VARIANTS_TXT" | awk '{print $1}'
echo "Total overlaps on CHM1"
wc -l temp.match1.txt | awk '{print $1}'
echo "Total overlaps on CHM13"
wc -l temp.match13.txt | awk '{print $1}'
echo "among these, number of uniq GATK IDs"
wc -l "$OUTPUT_DIR""gatkIDsWithMatchingPacBio.txt" | awk '{print $1}'


########## clean up
rm -f temp*
if [[ $1 == *.vcf.gz ]]; then
    rm -f "$GATKVCF"
fi

echo "#################################################"
echo "Done parsing GATK VCF:"
echo "  $INPUTVCF"
echo "#################################################"
echo
