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
caseID=$(echo "$NAME"|tr '[:upper:]' '[:lower:]')

TEST_INS_VAR=$2
TEST_DUP_VAR=$3
TEST_RPL_VAR=$4

PABIO_CHM_1_INS_VAR=$5
PABIO_CHM_13_INS_VAR=$6
PABIO_CHM_1_DEL_VAR=$7
PABIO_CHM_13_DEL_VAR=$8

OUTPUT_DIR=$9
if [[ ! $9 =~ .+/$ ]]; then
    OUTPUT_DIR+="/"
fi

IDNAMES=$1"IDS"
COMPARISON=$1"vsPacBio"
echo "#################################################"
echo "$NAME VS PacBio"
echo

# run simple intersection analysis which outputs variant id and matching PacBio POS (possible many)
ins_against_ins "$TEST_INS_VAR" "$PABIO_CHM_1_INS_VAR" "$COMPARISON" 1 "$caseID"
dup_against_ins "$TEST_DUP_VAR" "$PABIO_CHM_1_INS_VAR" "$COMPARISON" 1 "$caseID"
rpl_against_ins "$TEST_RPL_VAR" "$PABIO_CHM_1_INS_VAR" "$COMPARISON" 1 "$caseID"
rpl_against_del "$TEST_RPL_VAR" "$PABIO_CHM_1_DEL_VAR" "$COMPARISON" 1 "$caseID"

ins_against_ins "$TEST_INS_VAR" "$PABIO_CHM_13_INS_VAR" "$COMPARISON" 13 "$caseID"
dup_against_ins "$TEST_DUP_VAR" "$PABIO_CHM_13_INS_VAR" "$COMPARISON" 13 "$caseID"
rpl_against_ins "$TEST_RPL_VAR" "$PABIO_CHM_13_INS_VAR" "$COMPARISON" 13 "$caseID"
rpl_against_del "$TEST_RPL_VAR" "$PABIO_CHM_13_DEL_VAR" "$COMPARISON" 13 "$caseID"

# get IDS of test variant
awk '{print $3}' "$TEST_INS_VAR" > temp."$IDNAMES".ins.txt
awk '{print $3}' "$TEST_DUP_VAR" > temp."$IDNAMES".dup.txt
awk '{print $3}' "$TEST_RPL_VAR" > temp."$IDNAMES".rpl.txt
cat temp."$IDNAMES".ins.txt \
    temp."$IDNAMES".dup.txt \
    temp."$IDNAMES".rpl.txt | \
    sort \
    > temp.allTestVars.txt
echo "Total number of variants from $NAME to be tested:"
wc -l temp.allTestVars.txt | awk '{print $1}'

# get total overlaps from either cell lines
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

cat temp_1.txt temp_13.txt > temp.txt

# count number of uniq variants from test set and from PacBio calls
caseID+="IDs"
awk '{print $1}' temp.txt | \
    sort | \
    uniq \
    > "$OUTPUT_DIR$caseID""WithMatchingPacBio.txt"
comm -23 \
    temp.allTestVars.txt \
    "$OUTPUT_DIR$caseID""WithMatchingPacBio.txt" \
    > "$OUTPUT_DIR$caseID""NoMatchingPacBio.txt"
echo "  among these, number of uniq $NAME IDs"
wc -l "$OUTPUT_DIR$caseID""WithMatchingPacBio.txt" | awk '{print $1}'
echo "  and the number of variants without matching PacBio calls"
wc -l "$OUTPUT_DIR$caseID""NoMatchingPacBio.txt" | awk '{print $1}'

echo "  and the number of uniq PacBio call locations"
awk '{print $2}' temp.txt | \
    sort | \
    uniq \
    > "$OUTPUT_DIR""combined.$COMPARISON.uniqPacBioIDs.txt"
wc -l "$OUTPUT_DIR""combined.$COMPARISON.uniqPacBioIDs.txt" | awk '{print $1}'
echo

rm -f temp*
echo "#################################################"
echo "Done for $NAME"
echo
