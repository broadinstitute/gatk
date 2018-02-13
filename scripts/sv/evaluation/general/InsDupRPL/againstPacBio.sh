#!/bin/bash

set -eu

command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools for interval manipulation but it's not installed. Aborting."; exit 1; }

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source "$DIR""/utilities.sh"

if [[ "$#" -lt 9 ]]; then
    echo "Please provide:"
    echo "[1] name of call set: GATK or MANTA"
    echo "[2] absolute path to insertion VCF"
    echo "[3] absolute path to tandem duplication VCF"
    echo "[4] absolute path to scarred deletion (aka RPL) VCF"
    echo "[5] absolute path to CHM1 PacBio insertion VCF"
    echo "[6] absolute path to CHM13 PacBio insertion VCF"
    echo "[7] absolute path to CHM1 PacBio deletion VCF"
    echo "[8] absolute path to CHM13 PacBio deletion VCF"
    echo "[9] absolute path to the directory where analysis outputs are to be written to"
    exit 1
fi

NAME=$1

TEST_INS_VCF=$2
TEST_DUP_VCF=$3
TEST_RPL_VCF=$4

PABIO_CHM_1_INS_VCF=$5
PABIO_CHM_13_INS_VCF=$6
PABIO_CHM_1_DEL_VCF=$7
PABIO_CHM_13_DEL_VCF=$8

OUTPUT_DIR=$9
if [[ ! $9 =~ .+/$ ]]; then
    OUTPUT_DIR+="/"
fi

IDNAMES=$1"IDS"
COMPARISON=$1"vsPacbio"
echo "#################################################"
echo "$NAME VS PacBio"
echo

ins_against_ins "$TEST_INS_VCF" "$PABIO_CHM_1_INS_VCF" "CONTIG_DEPTH" "$COMPARISON" 1
dup_against_ins "$TEST_DUP_VCF" "$PABIO_CHM_1_INS_VCF" "CONTIG_DEPTH" "$COMPARISON" 1
rpl_against_ins "$TEST_RPL_VCF" "$PABIO_CHM_1_INS_VCF" "CONTIG_DEPTH" "$COMPARISON" 1
rpl_against_del "$TEST_RPL_VCF" "$PABIO_CHM_1_DEL_VCF" "CONTIG_DEPTH" "$COMPARISON" 1

ins_against_ins "$TEST_INS_VCF" "$PABIO_CHM_13_INS_VCF" "CONTIG_DEPTH" "$COMPARISON" 13
dup_against_ins "$TEST_DUP_VCF" "$PABIO_CHM_13_INS_VCF" "CONTIG_DEPTH" "$COMPARISON" 13
rpl_against_ins "$TEST_RPL_VCF" "$PABIO_CHM_13_INS_VCF" "CONTIG_DEPTH" "$COMPARISON" 13
rpl_against_del "$TEST_RPL_VCF" "$PABIO_CHM_13_DEL_VCF" "CONTIG_DEPTH" "$COMPARISON" 13


grep -v '^#' "$TEST_INS_VCF" | awk '{print $3}' > temp."$IDNAMES".ins.txt
grep -v '^#' "$TEST_DUP_VCF" | awk '{print $3}' > temp."$IDNAMES".dup.txt
grep -v '^#' "$TEST_RPL_VCF" | awk '{print $3}' > temp."$IDNAMES".rpl.txt
cat temp."$IDNAMES".ins.txt temp."$IDNAMES".dup.txt temp."$IDNAMES".rpl.txt | sort > temp."$IDNAMES".txt
echo "Total number of variants from $NAME to be tested:"
wc -l temp."$IDNAMES".txt | awk '{print $1}'
echo "Total overlaps on CHM1"
cat "$OUTPUT_DIR""insIns.""$COMPARISON""_1.window_w50.IDs.txt" \
    "$OUTPUT_DIR""dupIns.""$COMPARISON""_1.window_l0r50.IDs.txt" \
    "$OUTPUT_DIR""rplIns.""$COMPARISON""_1.window_w50.IDs.txt" \
    "$OUTPUT_DIR""rplDel.""$COMPARISON""_1.window_w50.IDs.txt" \
    > temp_1.txt
wc -l temp_1.txt | awk '{print $1}'

echo "Total overlaps on CHM13"
cat "$OUTPUT_DIR""insIns.""$COMPARISON""_13.window_w50.IDs.txt" \
    "$OUTPUT_DIR""dupIns.""$COMPARISON""_13.window_l0r50.IDs.txt" \
    "$OUTPUT_DIR""rplIns.""$COMPARISON""_13.window_w50.IDs.txt" \
    "$OUTPUT_DIR""rplDel.""$COMPARISON""_13.window_w50.IDs.txt" \
    > temp_13.txt
wc -l temp_13.txt | awk '{print $1}'


caseID=$(echo "$NAME"|tr '[:upper:]' '[:lower:]')
caseID+="IDs"
cat temp_1.txt temp_13.txt > temp.txt

echo "among these, number of uniq $NAME IDs"
awk '{print $1}' temp.txt | sort | uniq > "$OUTPUT_DIR$caseID""WithMatchingPacBio.txt"
wc -l "$OUTPUT_DIR$caseID""WithMatchingPacBio.txt" | awk '{print $1}'
comm -23 temp."$IDNAMES".txt "$OUTPUT_DIR$caseID""WithMatchingPacBio.txt" \
    > "$OUTPUT_DIR$caseID""NoMatchingPacBio.txt"

echo "and number of uniq PacBio IDs"
awk '{print $2}' temp.txt | sort | uniq > "$OUTPUT_DIR""combined.$COMPARISON.uniqPacbioIDs.txt"
wc -l "$OUTPUT_DIR""combined.$COMPARISON.uniqPacbioIDs.txt" | awk '{print $1}'
echo

rm -f temp*
echo "#################################################"
echo "Done for $NAME"
echo
