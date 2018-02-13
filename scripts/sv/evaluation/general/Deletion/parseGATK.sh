#!/bin/bash

set -eu

if [[ "$#" -lt 2 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to GATK VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
  echo "[3] (optional) mapping quality filter"
  echo "[4] (optional) alignment length filter"
	exit 1
fi

INPUTVCF=$1

OUTPUT_DIR=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR+="/"
fi

MQ_FT=20
LEN_FT=30
if [[ "$#" -eq 4 ]]; then
  MQ_FT=$3;
  LEN_FT=$4;
fi

echo "#################################################"
echo "GATK VCF parsing on:"
echo "  $INPUTVCF"
echo "#################################################"
echo

VARIANTS_TXT="temp.GATK_primaryContigs_cleanDel.txt"
grep -v '^#'  "$INPUTVCF" > "$VARIANTS_TXT"

awk '{print $8}' "$VARIANTS_TXT" > temp.info.txt
grep -Eo "MAX_ALIGN_LENGTH=[0-9]+" temp.info.txt | awk -F '=' '{print $2}' > temp.maxAlignLen.txt
grep -Eo "END=[0-9]+" temp.info.txt | awk -F '=' '{print $2}' > temp.end.txt
awk -F';' 'match($0, /HOMLEN=[0-9]+/){print substr($0,RSTART,RLENGTH);next};{print "HOMLEN=0"}' temp.info.txt | \
    awk -F '=' '{print $2}' > temp.homLen.txt
grep -Eo "MAPPING_QUALITIES=[0-9]+(,[0-9]+){0,}" temp.info.txt | \
    awk -F '=' '{print $2}' | tr , ' ' | \
    awk -v filter="$MQ_FT" '{m=$1;for(i=1;i<=NF;i++)if($i>m)m=$i; if(m>=filter) print "PASS"; else print "FAIL"}' > \
    temp.mqFilter.txt
awk -v filter="$LEN_FT" '{if ($1 >= filter) print "PASS"; else print "FAIL";}' temp.maxAlignLen.txt > \
    temp.alnLenFilter.txt

paste <(awk 'BEGIN {OFS="	"}; {print $1, $2}' "$VARIANTS_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.end.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $3}' "$VARIANTS_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.maxAlignLen.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.homLen.txt) \
      <(awk 'BEGIN {OFS=" "}; {print $1}' temp.mqFilter.txt) \
      <(awk 'BEGIN {OFS=" "}; {print $1}' temp.alnLenFilter.txt) \
      > temp.bed

awk 'BEGIN {OFS="	"}; {s=$4";"$5";"$6";"$7";"$8; print $1, $2, $3, s}' temp.bed \
    > "$OUTPUT_DIR""gatk.cleanDel.all.bed"

grep -vF 'FAIL' "$OUTPUT_DIR""gatk.cleanDel.all.bed" \
    > "$OUTPUT_DIR""gatk.cleanDel.bed"

rm -f temp*

echo "#################################################"
echo "Done parsing GATK VCF:"
echo "  $INPUTVCF"
echo "#################################################"
echo