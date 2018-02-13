#!/bin/bash

set -eu

if [[ "$#" -ne 3 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to Manta VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
	echo "[3] reference version (\"19\"|\"38\")"
	exit 1
fi

if [[ $1 == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$1")
    pattern=".gz"
    MANTAVCF=${COMPRESSED//$pattern/}
    bgzip -c -d "$1" > "$MANTAVCF"
elif [[ $1 == *.vcf  ]]; then
    MANTAVCF=$1
fi


OUTPUT_DIR_ROOT=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR_ROOT+="/"
fi

PRIMARY_CONTIGS_PATTERN=""
if [[ $3 == "19" ]]; then
	PRIMARY_CONTIGS_PATTERN="^([0-9]{1,2}|X|Y)	"
elif [[ $3 == "38" ]]; then
	PRIMARY_CONTIGS_PATTERN="^chr([0-9]{1,2}|X|Y)	"
else
	echo "reference version must be either 19 or 38"
	exit 1
fi

echo "#################################################"
echo "Checking and splitting Manta VCF:"
echo "  $1"
echo "#################################################"
echo

grep '^#' "$MANTAVCF" > "$OUTPUT_DIR_ROOT""MantaVCFHeader.txt"
grep -v '^#' "$MANTAVCF" > tempMantaVar.txt

########## size accounting on the raw and extraction for primary contigs
echo "num of variants"
wc -l tempMantaVar.txt | awk '{print $1}'

echo "num of filter-passing variants"
grep '\tPASS\t' tempMantaVar.txt > temp2.txt
wc -l temp2.txt | awk '{print $1}'

echo "num of filter-passing precise variants"
grep -v 'IMPRECISE' temp2.txt > tempMantaVar.txt
wc -l tempMantaVar.txt | awk '{print $1}'

echo "num of filter-passing, precise, non-BND variants"
grep -v '\tMantaBND' tempMantaVar.txt > temp2.txt
wc -l temp2.txt | awk '{print $1}'

echo "num of filter-passing, precise, non-BND variants on the primary contigs "
echo "  (\"primary\" defined as chr1-22 and chrX and chrY)"
grep -E "$PRIMARY_CONTIGS_PATTERN" temp2.txt > temp3.txt
wc -l temp3.txt | awk '{print $1}'
echo

MANTA_VARIANTS_OF_PRIMARY_INTEREST="$OUTPUT_DIR_ROOT""Manta_PASS_PRECISE_nonBND_primaryContigs.txt"
MANTA_STRANGE_RECORDS="$OUTPUT_DIR_ROOT""Manta_PASS_PRECISE_nonBND_primaryContigs_strangeSize.txt"
awk 'match($8, /SVLEN=(-)?[0-9]+/){sz=substr($8,RSTART+6,RLENGTH-6); sz+=0; sz_abs=(sz<0 ? -sz : sz); if (sz_abs<50 || sz_abs > 50000) print $0}' \
    temp3.txt > "$MANTA_STRANGE_RECORDS"
wc -l "$MANTA_STRANGE_RECORDS" | awk '{prin $1}'
echo "Any Manta records that have strange length values (<50, or >50K) saved in"
echo "  $MANTA_STRANGE_RECORDS"

awk '{print $3}' "$MANTA_STRANGE_RECORDS" > temp.mantaStrangeOnes.IDs.txt
grep -vf temp.mantaStrangeOnes.IDs.txt temp3.txt > "$MANTA_VARIANTS_OF_PRIMARY_INTEREST"
echo "num of filter-passing, precise, non-BND variants on the primary contigs with size [50, 50k]"
wc -l "$MANTA_VARIANTS_OF_PRIMARY_INTEREST" | awk '{print $1}'
echo

########## INV 
echo "num of inversions"
grep -E '\tMantaINV:([0-9]+(:)?){6,6}\t' \
    "$MANTA_VARIANTS_OF_PRIMARY_INTEREST" > \
    tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | \
    awk -F '=' '{print $2}' > \
    "$OUTPUT_DIR_ROOT""mantaExactInversionSizes.txt"
mkdir -p "$OUTPUT_DIR_ROOT""Inversion"
mv tempVar.txt \
    "$OUTPUT_DIR_ROOT""Inversion/Manta_PASS_PRECISE_nonBND_primaryContigs_inversions.txt"
echo 

########## clean DEL
echo "num of clean deletions:"
echo "    those having ID of MantaDEL and having neither a CIGAR of format \"1M[0-9]+I[0-9]+D\" nor having annotation of SVINSLEN/SVINSSEQ"
grep -E '\tMantaDEL:([0-9]+(:)?){6,6}\t' \
    "$MANTA_VARIANTS_OF_PRIMARY_INTEREST" | \
    grep -vE 'CIGAR=1M[0-9]+I[0-9]+D' | \
    grep -vF 'SVINSLEN' \
    > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=-[0-9]+' tempVar.txt | \
    awk -F '=-' '{print $2}' > \
    "$OUTPUT_DIR_ROOT""mantaExactCleanDeletionSizes.txt"
mkdir -p "$OUTPUT_DIR_ROOT""Deletion"
mv tempVar.txt \
    "$OUTPUT_DIR_ROOT""Deletion/Manta_PASS_PRECISE_nonBND_primaryContigs_cleanDel.txt"
echo 

########## manta long-range substitutions comes in two ways: those with CIGAR of format "1M[0-9]+I[0-9]+D" and those DEL without CIGAR but with SVINSSEQ & SVINSLEN
echo "num of long-range substitutions"
grep -E '(CIGAR=1M[0-9]+I[0-9]+D|	MantaDEL:([0-9]+(:)?){6,6}	.+SVINSLEN=)' \
    "$MANTA_VARIANTS_OF_PRIMARY_INTEREST" > \
    tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
echo "num of long-range substitutions called as deletions"
grep -E '\tMantaDEL:([0-9]+(:)?){6,6}\t' -c tempVar.txt
echo "num of long-range substitutions called as insertions"
grep -E '\tMantaINS:([0-9]+(:)?){6,6}\t' -c tempVar.txt
mkdir -p "$OUTPUT_DIR_ROOT""InsDupRPL"
mv tempVar.txt \
    "$OUTPUT_DIR_ROOT""InsDupRPL/Manta_PASS_PRECISE_nonBND_primaryContigs_longRangeSubstitution.txt"
echo 

########## INS 
echo "num of insertions with inferred size"
grep -E '\tMantaINS:([0-9]+(:)?){6,6}\t' \
    "$MANTA_VARIANTS_OF_PRIMARY_INTEREST" | \
    grep -vE 'CIGAR=1M[0-9]+I[0-9]+D' | \
    grep -F 'SVLEN=' > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | awk -F '=' '{print $2}' > \
    "$OUTPUT_DIR_ROOT""mantaExactCleanInsertionSizes.txt"
mkdir -p "$OUTPUT_DIR_ROOT""InsDupRPL"
mv tempVar.txt \
    "$OUTPUT_DIR_ROOT""InsDupRPL/Manta_PASS_PRECISE_nonBND_primaryContigs_cleanIns.txt"
echo 

########## DUP:TANDEM 
echo "num of tandem duplications"
grep -E '\tMantaDUP:TANDEM:([0-9]+(:)?){6,6}\t' "$MANTA_VARIANTS_OF_PRIMARY_INTEREST" | \
    grep -vE 'CIGAR=1M[0-9]+I[0-9]+D' | \
    grep -F 'SVLEN=' > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
awk 'match($8,/SVLEN=[0-9]+/){lenone=substr($8, RSTART+6, RLENGTH-6);lenone+=0}match($8,/SVINSLEN=[0-9]+/){lentwo=substr($8, RSTART+9, RLENGTH-9);lentwo+=0;next}{lentwo=0}{print lenone+lentwo}' tempVar.txt \
    > "$OUTPUT_DIR_ROOT""mantaExactCleanTandupSizes.txt"
mkdir -p "$OUTPUT_DIR_ROOT""InsDupRPL"
mv tempVar.txt \
    "$OUTPUT_DIR_ROOT""InsDupRPL/Manta_PASS_PRECISE_nonBND_primaryContigs_tandup.txt"
echo 

########## unknown size INS (those having ID of MantaINS and (LEFT|RIGHT)_SVINSSEQ annotation, hence no SVLEN annotation)
echo "num of insertion with unknown sizes"
grep -E '\tMantaINS:([0-9]+(:)?){6,6}\t' \
    "$MANTA_VARIANTS_OF_PRIMARY_INTEREST" | \
    grep -E '(LEFT|RIGHT)_SVINSSEQ' > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
mkdir -p "$OUTPUT_DIR_ROOT""UnknownSizedInsertion"
mv tempVar.txt \
    "$OUTPUT_DIR_ROOT""UnknownSizedInsertion/Manta_PASS_PRECISE_nonBND_primaryContigs_unknownSizedInsertion.txt"
echo 

########## clean up
rm -f temp*
if [[ $1 == *.vcf.gz ]]; then
    rm -f "$MANTAVCF"
fi


echo "#################################################"
echo "Done for Manta VCF"
echo "  $1"
echo "#################################################"
echo
