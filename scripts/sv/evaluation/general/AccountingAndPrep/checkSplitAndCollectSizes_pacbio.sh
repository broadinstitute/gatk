#!/bin/bash

set -eu

if [[ "$#" -ne 4 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to PacBio VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
    echo "[3] which halploid sample it is (\"0\"|\"1\"|\"13\") with 0 for mixture"
	echo "[4] reference version (\"19\"|\"38\")"
	exit 1
fi

if [[ $1 == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$1")
    pattern=".gz"
    PACBIOVCF=${COMPRESSED//$pattern/}
    bgzip -c -d "$1" > "$PACBIOVCF"
elif [[ $1 == *.vcf  ]]; then
    PACBIOVCF=$1
fi

OUTPUT_DIR_ROOT=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR_ROOT+="/"
fi


POSTFIX=$3


PRIMARY_CONTIGS_PATTERN=""
if [[ $4 == "19" ]]; then
	PRIMARY_CONTIGS_PATTERN="^([0-9]{1,2}|X|Y)	"
elif [[ $4 == "38" ]]; then
	PRIMARY_CONTIGS_PATTERN="^chr([0-9]{1,2}|X|Y)	"
else
	echo "reference version must be either 19 or 38"
    exit 1
fi

echo "#################################################"
echo "Checking and splitting PacBio VCF:"
echo "  $1"
echo "#################################################"
echo

grep '^#' "$PACBIOVCF" > "$OUTPUT_DIR_ROOT""PacBioVCFHeader.txt"

########## check some observations made when eye-balling and extraction for primary contigs
PB_VARIANTS_OF_PRIMARY_INTEREST="$OUTPUT_DIR_ROOT""PacBio_primaryContigs_var_""$POSTFIX"".txt"
echo "Total number of variants:"
grep -v '^#' -c "$PACBIOVCF"
echo

echo "check that all variants PASSes filter and have SVLEN"
FILTER=$(grep -v '^#' "$PACBIOVCF"| awk '{print $7}' | sort | uniq)
if [[ $FILTER == "PASS" ]]; then echo yes; else echo no; fi
ALL_HAS_LEN=$(grep -v '^#' "$PACBIOVCF"| awk 'match($8, /SVLEN=(-?)[0-9]+/){print 1;next;}{print 0}' | sort -n | uniq)
if [[ $ALL_HAS_LEN == 1 ]]; then echo yes; else echo no; fi
echo

echo "type of variants PacBio callset contains:"
grep -v '^#' "$PACBIOVCF" | \
    awk 'match($8, /SVTYPE=[a-zA-Z]+/){print substr($8, RSTART, RLENGTH); next;}{print "Unknown"}' | \
    sort | uniq
echo

grep -v '^#' "$PACBIOVCF" | grep -E "$PRIMARY_CONTIGS_PATTERN" > "$PB_VARIANTS_OF_PRIMARY_INTEREST"
echo "Total number of variants on the \"primary\" contigs:"
wc -l "$PB_VARIANTS_OF_PRIMARY_INTEREST" | awk '{print $1}'
echo

########## INV 
echo "num of inversions"
grep '\t<INV>\t' "$PB_VARIANTS_OF_PRIMARY_INTEREST" > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | \
    awk -F'=' '{print $2}' | \
    sort -n > "$OUTPUT_DIR_ROOT""pacbioInversionSizes_""$POSTFIX"".txt"
mkdir -p "$OUTPUT_DIR_ROOT""Inversion"
mv tempVar.txt "$OUTPUT_DIR_ROOT""Inversion/PacBio_primaryContigs_inversions_""$POSTFIX"".txt"
echo

########## DEL (Huddelston call set doesn't contain scarred deletions)
echo "num of deletions"
grep '\t<DEL>\t' "$PB_VARIANTS_OF_PRIMARY_INTEREST" > \
    tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | \
    awk -F'=' '{print $2}' | \
    sort -n > "$OUTPUT_DIR_ROOT""pacbioDeletionSizes_""$POSTFIX"".txt"
mkdir -p "$OUTPUT_DIR_ROOT""Deletion"
mv tempVar.txt "$OUTPUT_DIR_ROOT""Deletion/PacBio_primaryContigs_cleanDel_""$POSTFIX"".txt"
echo

########## INS 
echo "num of insertions"
grep '\t<INS>\t' "$PB_VARIANTS_OF_PRIMARY_INTEREST" > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | \
    awk -F'=' '{print $2}' | \
    sort -n > "$OUTPUT_DIR_ROOT""pacbioInsertionSizes_""$POSTFIX"".txt"
mkdir -p "$OUTPUT_DIR_ROOT""InsDupRPL"
mv tempVar.txt "$OUTPUT_DIR_ROOT""InsDupRPL/PacBio_primaryContigs_ins_""$POSTFIX"".txt"
echo

########## clean up
rm -f temp*
if [[ $1 == *.vcf.gz ]]; then
    rm -f "$PACBIOVCF"
fi

echo "#################################################"
echo "Done for PacBio VCF"
echo "  $1"
echo "#################################################"
echo
