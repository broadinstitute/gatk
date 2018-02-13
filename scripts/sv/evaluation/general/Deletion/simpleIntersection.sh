#!/bin/bash

set -eu

if [[ "$#" -lt 5 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to parsed GATK BED file"
	echo "[2] absolute path to parsed Manta BED file"
	echo "[3] absolute path to parsed PacBio BED file for CHM1"
	echo "[4] absolute path to parsed PacBio BED file for CHM13"
	echo "[5] absolute path to the directory where analysis outputs are to be written to"
	exit 1
fi

GATKBED=$1
MANTABED=$2

PACBIOBED_1=$3
PACBIOBED_13=$4

OUTPUT_DIR=$5
if [[ ! $5 =~ .+/$ ]]; then
	OUTPUT_DIR+="/"
fi


FRACTION=0.8
if [[ "$#" -eq 6 ]]; then
	FRACTION=$6
fi

POSTFIX=$(echo "${FRACTION/./}")
POSTFIX="f""$POSTFIX"
POSTFIX+="r"

echo "#################################################"
echo "Run simple intersection with reciprocal overlapping percentage:"
echo "  $FRACTION"
echo "with"
bedtools --version
echo "#################################################"
echo 

NUM_GATK=$(grep -v '^#' -c "$GATKBED")
NUM_MANTA=$(grep -v '^#' -c "$MANTABED")
echo "Number of GATK variants:" "$NUM_GATK"
echo "Number of Manta variants:" "$NUM_MANTA"
NUM_PACBIO=$(grep -v '^#' -c "$PACBIOBED_1")
echo "Number of PacBio variants for CHM1:" "$NUM_PACBIO"
NUM_PACBIO=$(grep -v '^#' -c "$PACBIOBED_13")
echo "Number of PacBio variants for CHM13:" "$NUM_PACBIO"
echo 

echo "######################"
echo "GATK vs Manta"

bedtools intersect \
    -a "$GATKBED" \
    -b "$MANTABED" \
    -f "$FRACTION" \
    -r \
    -wao | \
    grep -F 'Manta' > "$OUTPUT_DIR""cleanDel.gatkVSmanta.intersection_$POSTFIX.txt"
N=$(wc -l "$OUTPUT_DIR""cleanDel.gatkVSmanta.intersection_$POSTFIX.txt" | awk '{print $1}')
echo "Number of GATK records with overlapping Manta clean deletion call:"
echo "$N"
echo "Number of GATK records without overlapping Manta clean deletion call: "
expr "$NUM_GATK" - "$N"
echo

echo "######################"
echo "GATK vs Pacbio"

bedtools intersect \
    -a "$GATKBED" \
    -b "$PACBIOBED_1" \
    -f "$FRACTION" \
    -r \
    -wao | \
    grep -F 'PB' > "$OUTPUT_DIR""cleanDel.gatkVSpacbio.1.intersection_$POSTFIX.txt"
N=$(wc -l "$OUTPUT_DIR""cleanDel.gatkVSpacbio.1.intersection_$POSTFIX.txt" | awk '{print $1}')
echo "Number of GATK records with overlapping PacBio clean deletion call from CHM1: "
echo "$N"

bedtools intersect \
    -a "$GATKBED" \
    -b "$PACBIOBED_13" \
    -f "$FRACTION" \
    -r \
    -wao | \
    grep -F 'PB' > "$OUTPUT_DIR""cleanDel.gatkVSpacbio.13.intersection_$POSTFIX.txt"
N=$(wc -l "$OUTPUT_DIR""cleanDel.gatkVSpacbio.13.intersection_$POSTFIX.txt" | awk '{print $1}')
echo "Number of GATK records with overlapping PacBio clean deletion call from CHM13: "
echo "$N"

awk '{print $4}' "$GATKBED" | awk -F ';' '{print $1}' | sort > temp.gatkIDs.txt
awk '{print $4}' "$OUTPUT_DIR""cleanDel.gatkVSpacbio.1.intersection_$POSTFIX.txt" | awk -F ';' '{print $1}' > temp.gatkIDsMatchingCHM1.txt
awk '{print $4}' "$OUTPUT_DIR""cleanDel.gatkVSpacbio.13.intersection_$POSTFIX.txt" | awk -F ';' '{print $1}' > temp.gatkIDsMatchingCHM13.txt
cat temp.gatkIDsMatchingCHM1.txt temp.gatkIDsMatchingCHM13.txt | sort | uniq > "$OUTPUT_DIR""gatkIDsWithMatchingPacBio.txt"
echo "Number of GATK records with overlapping PacBio clean deletion call from either haploid: "
N=$(wc -l "$OUTPUT_DIR""gatkIDsWithMatchingPacBio.txt" | awk '{print $1}')
echo "$N"
echo "Number of GATK records without overlapping PacBio clean deletion call: "
expr "$NUM_GATK" - "$N"
comm -3 temp.gatkIDs.txt "$OUTPUT_DIR""gatkIDsWithMatchingPacBio.txt" > "$OUTPUT_DIR""gatkIDsNoMatchingPacBio.txt"
echo "They are saved in: $OUTPUT_DIR""gatkIDsNoMatchingPacBio.txt"
echo

echo "######################"
echo "Manta vs PacBio"

bedtools intersect \
    -a "$MANTABED" \
    -b "$PACBIOBED_1" \
    -f "$FRACTION" \
    -r \
    -wao | \
    grep -F 'PB' > "$OUTPUT_DIR""cleanDel.mantaVSpacbio.1.intersection_$POSTFIX.txt"
N=$(wc -l "$OUTPUT_DIR""cleanDel.mantaVSpacbio.1.intersection_$POSTFIX.txt" | awk '{print $1}')
echo "Number of MANTA records with overlapping PacBio clean deletion call from CHM1: "
echo "$N"

bedtools intersect \
    -a "$MANTABED" \
    -b "$PACBIOBED_13" \
    -f "$FRACTION" \
    -r \
    -wao | \
    grep -F 'PB' > "$OUTPUT_DIR""cleanDel.mantaVSpacbio.13.intersection_$POSTFIX.txt"
N=$(wc -l "$OUTPUT_DIR""cleanDel.mantaVSpacbio.13.intersection_$POSTFIX.txt" | awk '{print $1}')
echo "Number of MANTA records with overlapping PacBio clean deletion call from CHM13: "
echo "$N"

awk '{print $4}' "$MANTABED" | awk -F ';' '{print $1}' | sort > temp.mantaIDs.txt
awk '{print $4}' "$OUTPUT_DIR""cleanDel.mantaVSpacbio.1.intersection_$POSTFIX.txt" | awk -F ';' '{print $1}' > temp.mantaIDsMatchingCHM1.txt
awk '{print $4}' "$OUTPUT_DIR""cleanDel.mantaVSpacbio.13.intersection_$POSTFIX.txt" | awk -F ';' '{print $1}' > temp.mantaIDsMatchingCHM13.txt
cat temp.mantaIDsMatchingCHM1.txt temp.mantaIDsMatchingCHM13.txt | sort | uniq > "$OUTPUT_DIR""mantaIDsWithMatchingPacBio.txt"
echo "Number of MANTA records with overlapping PacBio clean deletion call from either haploid: "
N=$(wc -l "$OUTPUT_DIR""mantaIDsWithMatchingPacBio.txt" | awk '{print $1}')
echo "$N"
echo "Number of MANTA records without overlapping PacBio clean deletion call: "
expr "$NUM_MANTA" - "$N"
comm -3 temp.mantaIDs.txt "$OUTPUT_DIR""mantaIDsWithMatchingPacBio.txt" > "$OUTPUT_DIR""mantaIDsNoMatchingPacBio.txt"
echo "They are saved in: $OUTPUT_DIR""mantaIDsNoMatchingPacBio.txt"
echo


rm -f temp*

echo "#################################################"
echo "Done performing simple intersection on the bed files."
echo "#################################################"
echo