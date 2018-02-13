#!/bin/bash

set -eu

if [[ "$#" -lt 2 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to PacBio VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
	exit 1
fi

PACBIOVCF=$1

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

OUTPUT_DIR_ROOT=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR_ROOT+="/"
fi

echo "#################################################"
echo "Checking against several observations made while looking at the PacBio deletion records on:"
echo "  $PACBIOVCF"
echo "#################################################"
echo

VARIANTS_TXT="temp.PACBIO_primaryContigs_cleanDel.txt"
grep -v '^#'  "$PACBIOVCF" > "$VARIANTS_TXT"
NUM_VAR=$(wc -l "$VARIANTS_TXT" | awk '{print $1}')

echo "First check that all PacBio records have CONTIG_DEPTH annotation"
N2=$(grep -E 'CONTIG_DEPTH=[0-9]+' -c "$VARIANTS_TXT")
if [[ "$NUM_VAR" != "$N2" ]]; then
	echo "Not all PacBio records have CONTIG_DEPTH annotation";
	exit 1;
fi

echo "Then check that all PacBio records have CONTIG_SUPPORT annotation"
N2=$(grep -E 'CONTIG_SUPPORT=[0-9]+' -c "$VARIANTS_TXT")
if [[ "$NUM_VAR" != "$N2" ]]; then
	echo "Not all PacBio records have CONTIG_SUPPORT annotation";
	exit 1;
fi

echo "Then check that all PacBio records have SEQ annotation"
N2=$(grep -vE 'SEQ=[a-zA-Z]+' -c "$VARIANTS_TXT" || true)
if [[ "$N2" != 0 ]]; then
	echo "Some PacBio records don't have SEQ annotation"
	exit 1
fi

echo "Then check that all PacBio records' SEQ annotation has the same length as SVLEN"
grep -oE 'SVLEN=[0-9]+' "$VARIANTS_TXT" | awk -F '=' '{print $2}' > temp.svlen.txt
grep -oE 'SEQ=[a-zA-Z]+' "$VARIANTS_TXT" | awk -F '=' '{print length($2)}'  > temp.seq.txt
paste <(awk 'BEGIN {OFS="	"}; {print $1}' temp.svlen.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.seq.txt) \
      > temp.txt
N=$(awk '{print $1-$2}' temp.txt | sort -n | uniq)
if [[ $N != 0 ]]; then
	echo "Some PacBio deletion records have their SEQ annotation's length not equaling their SVLEN"
	exit 1
fi
echo "Then check that all PacBio records' END-POS+1 annotation's value has the same value as SVLEN"
awk 'match($8, /END=[0-9]+/){end=substr($8,RSTART+4,RLENGTH-4);end+=0;print end-$2+1}' "$VARIANTS_TXT" > temp.dist.txt
paste <(awk 'BEGIN {OFS="	"}; {print $1}' temp.svlen.txt) \
      <(awk 'BEGIN {OFS="	"}; {print $1}' temp.dist.txt) \
      > temp.txt
N=$(awk '{print $1-$2}' temp.txt | sort -n | uniq)
if [[ $N != 0 ]]; then
	echo "Some PacBio deletion records have their END-POS+1 annotation's value not equaling their SVLEN"
	exit 1
fi

echo "And the REPEAT_TYPE annotation takes the following values: "
echo "  (REPEAT_TYPE=NA, if printed, is added post-hoc to show that this annotation is missing for some records)"
awk 'match($8, /REPEAT_TYPE=[A-Za-z0-9]+(_[A-Za-z0-9]+)?;/){print substr($8,RSTART,RLENGTH-1);next;}{print "REPEAT_TYPE=NA"}' "$VARIANTS_TXT" | \
    sort | uniq

rm -f temp*

echo "#################################################"
echo "Done checking on"
echo "  $PACBIOVCF"
echo "#################################################"
echo