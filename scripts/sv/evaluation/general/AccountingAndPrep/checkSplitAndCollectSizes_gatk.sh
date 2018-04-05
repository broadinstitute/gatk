#!/bin/bash

set -eu

if [[ "$#" -ne 5 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to GATK VCF"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
    echo "[3] absolute path to GATK VCF containing simple variants re-interpreted from 1-seg CPX calls"
    echo "[4] absolute path to GATK VCF containing simple variants re-interpreted from multi-seg CPX calls"
	echo "[5] reference version (\"19\"|\"38\")"
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

OUTPUT_DIR=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR+="/"
fi

mkdir -p "$OUTPUT_DIR"


GATK_one_seg_VCF=""
if [[ $3 == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$3")
    pattern=".gz"
    GATK_one_seg_VCF=${COMPRESSED//$pattern/}
    bgzip -c -d "$3" > "$GATK_one_seg_VCF"
elif [[ $3 == *.vcf  ]]; then
    GATK_one_seg_VCF=$3
fi

GATK_multi_seg_VCF=""
if [[ $4 == *.vcf.gz ]]; then
    COMPRESSED=$(basename "$4")
    pattern=".gz"
    GATK_multi_seg_VCF=${COMPRESSED//$pattern/}
    bgzip -c -d "$4" > "$GATK_multi_seg_VCF"
elif [[ $4 == *.vcf  ]]; then
    GATK_multi_seg_VCF=$4
fi


PRIMARY_CONTIGS_PATTERN=""
if [[ $5 == "19" ]]; then
	PRIMARY_CONTIGS_PATTERN="^([0-9]{1,2}|X|Y)	"
elif [[ $5 == "38" ]]; then
	PRIMARY_CONTIGS_PATTERN="^chr([0-9]{1,2}|X|Y)	"
else
	echo "reference version must be either 19 or 38"
    exit 1
fi

echo "#################################################"
echo "Checking and splitting GATK VCF:"
echo "  $1"
echo "Number of variants in file"
grep -v '^#' -c "$GATKVCF"
echo "#################################################"
echo

grep '^#' "$GATKVCF" > "$OUTPUT_DIR""GATKVCFHeader.txt"
grep -v '^#' "$GATKVCF" > tempGATKRawVariants.txt

########## extraction for primary contigs
GATK_PRIME_VAR="$OUTPUT_DIR""GATK_primaryContigs_var.txt"

echo "##########"
echo "Extracting variants on primary contigs to"
echo "  $GATK_PRIME_VAR"
grep -E "$PRIMARY_CONTIGS_PATTERN" tempGATKRawVariants.txt > "$GATK_PRIME_VAR"
echo "Number of variants in" "$GATK_PRIME_VAR"
wc -l "$GATK_PRIME_VAR" | awk '{print $1}'
echo

########## check for GATK vcf for strange length (too small or large) values
GATK_STRANGE_RECORDS="$OUTPUT_DIR""WARN_gatk_strange_records.txt"
GATK_PRIME_VAR_NO_WARN="$OUTPUT_DIR""GATK_primaryContigs_var_no_warn.txt"

echo "##########"
echo "Checking on length annotation of GATK VCF"
echo "Any GATK records that have strange length values (<50, or >50K) saved in"
echo "  $GATK_STRANGE_RECORDS"
awk 'match($8, /SVLEN=(-)?[0-9]+/){sz=substr($8,RSTART+6,RLENGTH-6); sz+=0; sz_abs=(sz<0 ? -sz : sz); if (sz_abs<50 || sz_abs > 50000) print $0}' \
    "$GATK_PRIME_VAR" | \
    grep -vF 'BND' > "$GATK_STRANGE_RECORDS"

awk '{print $3}' "$GATK_STRANGE_RECORDS" > temp.gatkStrangeOnes.IDs.txt

grep -vf temp.gatkStrangeOnes.IDs.txt "$GATK_PRIME_VAR" > "$GATK_PRIME_VAR_NO_WARN"
echo "Number of variants in" "$GATK_PRIME_VAR_NO_WARN"
wc -l "$GATK_PRIME_VAR_NO_WARN" | awk '{print $1}'
echo

########## check for duplicate records
GATK_PRIME_VAR_NO_WARN_DUPLICATED_RECORDS="$OUTPUT_DIR""GATK_primaryContigs_var_no_warn_Duplicates.txt"
GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS="$OUTPUT_DIR""GATK_primaryContigs_var_no_warn_noDuplicates.txt"

echo "##########"
echo "There might be records in GATK's VCF that are \"duplicated\"."
echo " That is, there are GATK records with the same POS and END (hence the same ID) fields."
echo " We split the GATK variants between the \"duplicated\" and the \"non-duplicated\" ones,"
echo " and save the \"duplicated\" variants to $GATK_PRIME_VAR_NO_WARN_DUPLICATED_RECORDS"
echo " Please manually review the \"duplicated\" entries, to verify the assertion that "
echo " they differ, very slightly, in their INSSEQ and/or DUPLICATED_SEQUENCE annotations." 

awk '{print $3}' "$GATK_PRIME_VAR_NO_WARN" | uniq -d > temp.duplicatedGATKIDs.ins.txt

grep -vf temp.duplicatedGATKIDs.ins.txt \
    "$GATK_PRIME_VAR_NO_WARN" \
    > "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS"
echo "Number of variants in" "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS"
wc -l "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" | awk '{print $1}'
echo

parallel_grep 10 temp.duplicatedGATKIDs.ins.txt \
    "$GATK_PRIME_VAR_NO_WARN" \
    "$GATK_PRIME_VAR_NO_WARN_DUPLICATED_RECORDS"
echo "Number of variants in" "$GATK_PRIME_VAR_NO_WARN_DUPLICATED_RECORDS"
wc -l "$GATK_PRIME_VAR_NO_WARN_DUPLICATED_RECORDS" | awk '{print $1}'
echo

echo "#################################################"
echo "Done checking on"
echo "  $1"

if [[ ! -z "${GATK_one_seg_VCF}" ]] && [[ ! -z "${GATK_multi_seg_VCF}" ]]; then

    echo "Now merging $GATKVCF with $GATK_one_seg_VCF and $GATK_multi_seg_VCF"
    echo

    grep -v '^#' "$GATK_one_seg_VCF" > temp.oneseg.txt
    grep -v '^#' "$GATK_multi_seg_VCF" > temp.multiseg.txt

    cp "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" temp.major.txt
    cat temp.major.txt temp.oneseg.txt temp.multiseg.txt | \
       gsort -V \
       > temp.merged.txt
    mv temp.merged.txt "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS"

    echo "Done merging"
fi
echo "#################################################"
echo

echo "#################################################"
echo "Diverts different types of variants in GATK into different VCF files"
echo "  and collects & saves variants sizes by type"
echo "#################################################"
echo
mkdir -p "$OUTPUT_DIR""Accounting"


########## INV 
echo "Number of inversions"
grep '\t<INV>\t' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" > tempVar.txt || true
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | \
    awk -F'=' '{print $2}' | \
    sort -n > "$OUTPUT_DIR""Accounting/gatkInversionSizes.txt"
mkdir -p "$OUTPUT_DIR""Inversion"
mv tempVar.txt "$OUTPUT_DIR""Inversion/GATK_primaryContigs_inversions.txt"
echo

########## imprecise DEL
echo "Number of imprecise deletions"
grep '\t<DEL>\t' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" | \
    grep -F ';IMPRECISE;' > tempVar.txt || true
wc -l tempVar.txt | awk '{print $1}'
mkdir -p "$OUTPUT_DIR""ImpreciseDeletion"
mv tempVar.txt "$OUTPUT_DIR""ImpreciseDeletion/GATK_primaryContigs_impreciseDel.txt"
echo

########## scarred DEL
echo "Number of scarred deletions"
grep '\t<DEL>\t' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" | \
    grep -vF ';IMPRECISE;' | \
    grep -F 'INSSEQ=' > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
mkdir -p "$OUTPUT_DIR""InsDupRPL"
mv tempVar.txt "$OUTPUT_DIR""InsDupRPL/GATK_primaryContigs_scarredDel.txt"
echo

########## clean DEL
echo "Number of clean deletions"
grep '\t<DEL>\t' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" | \
    grep -vF 'INSSEQ=' | \
    grep -vF ';IMPRECISE;' > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=-[0-9]+' tempVar.txt | \
    awk -F'=-' '{print $2}' | \
    sort -n > "$OUTPUT_DIR""Accounting/gatkCleanDeletionSizes.txt"
mkdir -p "$OUTPUT_DIR""Deletion"
mv tempVar.txt "$OUTPUT_DIR""Deletion/GATK_primaryContigs_cleanDel.txt"
echo

########## INS 
echo "Number of insertions"
grep '\t<INS>\t' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | \
    awk -F'=' '{print $2}' | \
    sort -n > "$OUTPUT_DIR""Accounting/gatkInsertionSizes.txt"
mkdir -p "$OUTPUT_DIR""InsDupRPL"
mv tempVar.txt "$OUTPUT_DIR""InsDupRPL/GATK_primaryContigs_ins.txt"
echo

########## DUP:TANDEM 
echo "Number of tandem duplications"
grep '\t<DUP>\t' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" > tempVar.txt
wc -l tempVar.txt | awk '{print $1}'
grep -oE 'SVLEN=[0-9]+' tempVar.txt | \
    awk -F'=' '{print $2}' | \
    sort -n > "$OUTPUT_DIR""Accounting/gatkTandupSizes.txt"
mkdir -p "$OUTPUT_DIR""InsDupRPL"
mv tempVar.txt "$OUTPUT_DIR""InsDupRPL/GATK_primaryContigs_dup.txt"
echo


########## BND
echo "Number of BND records (counting both MATE's)"
NUM=$(grep -c '\tBND_chr' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" || true)
if [[ $NUM != 0 ]]; then
    grep '\tBND_chr' "$GATK_PRIME_VAR_NO_WARN_UNIQUE_RECORDS" > tempVar.txt
    wc -l tempVar.txt | awk '{print $1}'
    mkdir -p "$OUTPUT_DIR""BND"
    mv tempVar.txt "$OUTPUT_DIR""BND/GATK_primaryContigs_bnd.txt"
else 
    echo "0"
fi
echo 

########## clean up
rm -f temp*
if [[ $1 == *.vcf.gz ]]; then
    rm -f "$GATKVCF"
fi

echo "#################################################"
echo "Done for GATK VCF"
echo "  $1"
echo "#################################################"
echo
