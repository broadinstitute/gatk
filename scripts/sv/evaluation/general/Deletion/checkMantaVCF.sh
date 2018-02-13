#!/bin/bash

set -eu

if [[ "$#" -lt 2 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to Manta VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
	exit 1
fi

UNIQ="uniq"
if [[ "$OSTYPE" == "darwin"* ]]; then
	command -v guniq >/dev/null 2>&1 || { echo >&2 "I require guniq for sorting output bed file but it's not installed.  Aborting."; exit 1; }
        UNIQ="guniq"
elif [[ "$OSTYPE" == "linux-gnu" ]]; then
        UNIQ="uniq"
else
        echo "unsupported OS, quit"
        exit 1
fi

MANTAVCF=$1

echo "#################################################"
echo "Checking against several observations made while looking at the Manta deletion records on:"
echo "  $MANTAVCF"
echo "#################################################"
echo

VARIANTS_TXT="temp.Manta103_PASS_PRECISE_nonBND_primaryContigs_cleanDel.txt"
grep -v '^#'  "$MANTAVCF" > "$VARIANTS_TXT"

echo "check assumptions on observations made during eye-ball review:"

echo "  * the annotation CIPOS/CIEND is always presented with the format \"CIPOS=0,[0-9]+\" or \"CIEND=0,[0-9]+\""
STATUS=$(grep -Eo 'CI(POS|END)=[0-9]+,[0-9]+' "$VARIANTS_TXT" | awk -F'=' '{print $2}' | awk -F',' '{print $1}' | sort -n | uniq)
if [[ $STATUS != 0 ]]; then
	echo "CIPOS and CIEND doesn't always have the format \"CIPOS=0,[0-9]+\" or \"CIEND=0,[0-9]+\"";
	exit 1;
fi
echo "  * some deletions are annotated with CIPOS, but not all are annotated with CIEND"
echo "    those having CIPOS but not CIEND all have concrete REF allele and single base, i.e. the anchor base, ALT allele"
SEQS=$(grep -F 'CIPOS=' "$VARIANTS_TXT" | grep -vF 'CIEND=' | awk '{print $5}' | sort | uniq | xargs echo)
if [[ $SEQS != "A C G T" ]]; then 
	echo "Assertion that deletions having CIPOS but not CIEND all have single base ALT allele is false";
	exit 1;
fi
echo "  * in the mean time, deletions having CIPOS all have HOMESEQ and HOMLEN annotation"
N1=$(grep -F 'CIPOS=' -c "$VARIANTS_TXT")
N2=$(grep -F 'CIPOS=' "$VARIANTS_TXT" | grep -F 'HOMLEN' -c)
if [[ "$N1" != "$N2" ]]; then
	echo "CIPOS doesn't always have accompanying HOMLEN annotation";
	exit 1;
fi
echo "  * and the HOMLEN value always equals the second integer in the CIPOS annotation"
grep -Eo 'CIPOS=0,[0-9]+' "$VARIANTS_TXT" | awk -F',' '{print $2}' > temp.txt
grep -F 'CIPOS=' "$VARIANTS_TXT" | grep -Eo 'HOMLEN=[0-9]+' | awk -F'=' '{print $2}' > temp2.txt
STATUS=$(diff temp.txt temp2.txt | wc | awk '{print $1}')
if [[ $STATUS != 0 ]]; then
	echo "CIPOS doesn't always have the same length as indicated in HOMLEN";
	exit 1;
fi
echo "  * similarly, when both CIPOS and CIEND are present, they are the same"
grep -E 'CIEND=0,[0-9]+' "$VARIANTS_TXT" | grep -Eo 'CIPOS=0,[0-9]+' | awk -F'=' '{print $2}' > temp.txt
grep -Eo 'CIEND=0,[0-9]+' "$VARIANTS_TXT" | awk -F'=' '{print $2}' > temp2.txt
STATUS=$(diff temp.txt temp2.txt | wc | awk '{print $1}')
if [[ $STATUS != 0 ]]; then
	echo "CIPOS doesn't always have the same length as indicated in CIEND";
	exit 1;
fi
echo

rm -f temp*

echo "#################################################"
echo "Done checking on"
echo "  $MANTAVCF"
echo "#################################################"
echo