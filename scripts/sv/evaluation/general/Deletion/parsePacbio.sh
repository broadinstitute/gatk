#!/bin/bash

set -eu

if [[ "$#" -ne 3 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to PacBio VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
  echo "[3] which halploid sample it is (\"0\"|\"1\"|\"13\") with 0 for mixture"
	exit 1
fi

INPUTVCF=$1

OUTPUT_DIR=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR+="/"
fi

POSTFIX=$3

echo "#################################################"
echo "PacBio VCF parsing on:"
echo "  $INPUTVCF"
echo "#################################################"
echo

VARIANTS_TXT="temp.PacBio_primaryContigs_cleanDel_""$POSTFIX"".txt"
grep -v '^#' "$INPUTVCF" > "$VARIANTS_TXT"

awk '{print $5}' "$VARIANTS_TXT" | grep -Eo '[A-Z]+' > temp.type.txt
grep -Eo "	END=[0-9]+" "$VARIANTS_TXT" | awk -F '=' '{print $2}' > temp.end.txt
grep -Eo "CONTIG_DEPTH=[0-9]+" "$VARIANTS_TXT" | awk -F '=' '{print $2}' > temp.contigDepth.txt
grep -Eo "CONTIG_SUPPORT=[0-9]+" "$VARIANTS_TXT" | awk -F '=' '{print $2}' > temp.contigSupport.txt
awk 'match($8, /REPEAT_TYPE=[A-Za-z0-9]+(_[A-Za-z0-9]+)?;/){print substr($8,RSTART,RLENGTH-1);next;}{print "REPEAT_TYPE=NA"}' \
    "$VARIANTS_TXT" | awk -F '=' '{print $2}' > temp.repeatType.txt

paste <(awk 'BEGIN {OFS="	"}; {print $1, $2}' "$VARIANTS_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.end.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $6}' "$VARIANTS_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.contigDepth.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.contigSupport.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.repeatType.txt) \
      <(awk 'BEGIN {OFS=" "}; {print $1}' temp.type.txt) \
      > temp.bed

awk 'BEGIN {OFS="	"}; {s="PB_"$8";"$4";"$5";"$6";"$7; print $1, $2, $3, s}' temp.bed \
    > "$OUTPUT_DIR""pacbio.cleanDel_""$POSTFIX"".bed"

rm -f temp*

echo "#################################################"
echo "Done parsing PacBio VCF:"
echo "  $INPUTVCF"
echo "#################################################"
echo