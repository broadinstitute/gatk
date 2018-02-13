#!/bin/bash

set -eu

if [[ "$#" -ne 2 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to Manta VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
	exit 1
fi

SORT="sort"
if [[ "$OSTYPE" == "darwin"* ]]; then
	command -v gsort >/dev/null 2>&1 || { echo >&2 "I require gsort for sorting output bed file but it's not installed.  Aborting."; exit 1; }
	SORT="gsort"
elif [[ "$OSTYPE" == "linux-gnu" ]]; then
        SORT="sort"
else
	echo "unsupported OS, quit"
	exit 1
fi

INPUTVCF=$1

OUTPUT_DIR=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR+="/"
fi

echo "#################################################"
echo "Manta VCF parsing on:"
echo "  $INPUTVCF"
echo "#################################################"
echo

VARIANTS_TXT="temp.Manta_PASS_PRECISE_nonBND_primaryContigs_cleanDel.txt"
grep -v '^#'  "$INPUTVCF" > "$VARIANTS_TXT"

########## non-symbolic deletions (those with explicit REF allele sequence in the VCF)
echo "First deal with Manta clean deletions that have explicit REF allele (i.e. non-symbolic ALT allele)."
NONSYMB_TXT="temp.Manta_PASS_PRECISE_nonBND_primaryContigs_cleanDel.nonSymb.txt"
grep -v '<DEL>' "$VARIANTS_TXT" > "$NONSYMB_TXT"
awk '{print $10}' "$NONSYMB_TXT" | \
    awk -F':' 'BEGIN {OFS="	"}; {print $1, $3}' > temp.gtgq.txt
awk '{print $8}'  "$NONSYMB_TXT" | \
    awk -F';' 'BEGIN {OFS="	"}; match($0, /HOMLEN=[0-9]+/){print $1, substr($0,RSTART,RLENGTH);next;}{print $1, "HOMLEN=0"}' > \
    temp.endHomLen.txt

# BED4 file generation
paste <(awk 'BEGIN {OFS="	"}; {print $1, $2}' "$NONSYMB_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.endHomLen.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $3, $6}' "$NONSYMB_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1, $2}' temp.gtgq.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $2}' temp.endHomLen.txt) \
      > temp.bed

sed -i.bak 's/END=//' temp.bed && sed -i.bak 's/HOMLEN=//' temp.bed
awk 'BEGIN {OFS="	"}; {s=$4";"$5";"$6";"$7";"$8; print $1, $2, $3, s}' temp.bed \
    > "$OUTPUT_DIR""manta.nonSymbDel.bed"

echo
########## symbolic deletions (those with ALT allele sequence in the VCF as <DEL>)
echo "Then deal with Manta clean deletions that have symbolic ALT allele."

SYMB_TXT="temp.Manta_PASS_PRECISE_nonBND_primaryContigs_cleanDel.symb.txt"
grep -F '<DEL>' "$VARIANTS_TXT" > "$SYMB_TXT"

awk '{print $10}' "$SYMB_TXT" | \
    awk -F':' 'BEGIN {OFS="	"}; {print $1, $3}' > temp.gtgq.txt
awk '{print $8}'  "$SYMB_TXT" | \
    awk -F';' 'BEGIN {OFS="	"}; match($0, /HOMLEN=[0-9]+/){print $1, substr($0,RSTART,RLENGTH);next;}{print $1, "HOMLEN=0"}' \
    > temp.endHomLen.txt

# BED4 file generation
paste <(awk 'BEGIN {OFS="	"}; {print $1, $2}' "$SYMB_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.endHomLen.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $3, $6}' "$SYMB_TXT") \
      <(awk 'BEGIN {OFS="	"}; {print $1, $2}' temp.gtgq.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $2}' temp.endHomLen.txt) \
      > temp.bed

sed -i.bak 's/END=//' temp.bed && sed -i.bak 's/HOMLEN=//' temp.bed
awk 'BEGIN {OFS="	"}; {s=$4";"$5";"$6";"$7";"$8; print $1, $2, $3, s}' temp.bed \
    > "$OUTPUT_DIR""manta.symbDel.bed"
cat "$OUTPUT_DIR""manta.nonSymbDel.bed" \
    "$OUTPUT_DIR""manta.symbDel.bed" | \
    $SORT \
    > "$OUTPUT_DIR""manta.cleanDel.bed"

rm -f temp*

echo "#################################################"
echo "Done parsing Manta VCF:"
echo "  $INPUTVCF"
echo "#################################################"
echo
