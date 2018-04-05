#!/bin/bash

set -eu

if [[ "$#" -lt 4 ]]; then
	echo "Please provide:"
	echo "[1] absolute path to the GATK intervals file"
	echo "[2] absolute path to the directory where analysis outputs are to be written to"
	echo "[3] absolute path to intersection file between GATK and Manta variants"
	echo "[4] absolute path to Manta deletion VCF, or, if available"
	echo "[4] (if available) absolute path to intersection file between GATK and PacBio variants"
	echo "[5] (if available) absolute path to intersection file between Manta and PacBio variants"
	exit 1
fi

GATK_INTERVALS=$1

OUTPUT_DIR_ROOT=$2
if [[ ! $2 =~ .+/$ ]]; then
	OUTPUT_DIR_ROOT+="/"
fi

mkdir -p "$OUTPUT_DIR_ROOT""Results"
OUTPUT_DIR="$OUTPUT_DIR_ROOT""Results/"

G_VS_M=$3
if [[ "$#" -eq 4 ]]; then
	MANTAVCF=$4
	echo "$MANTAVCF"
elif [[ "$#" -eq 5 ]]; then
	G_VS_P=$4
	M_VS_P=$5
fi

echo "#################################################"
echo "Checking GATK local assembly trigger sensitivity"
echo "#################################################"
echo 

echo "Next check if local assemblies are triggered in the GATK-SV pipeline at the breakpoints identified by Manta AND by PacBio."
echo "Simple venn diagram approach"
echo
echo "First convert interval output file from FindBreakpointEvidenceSpark to a bed file for easier manipulation."
GATK_INTERVALS_BED="$OUTPUT_DIR""gatk.fastqIntervals.bed"
awk 'BEGIN {OFS="	"} match($3, /produced/) {print $2, sprintf("assembly%06d.fastq", $1)}' \
    "$GATK_INTERVALS" | \
    sed 's/\(.*\):/\1	/' | \
    sed 's/\(.*\)-/\1	/' > \
    "$GATK_INTERVALS_BED"
echo -e "Next compute simple descriptive statistics of the interval sizes"
awk '{print $3-$2}' "$OUTPUT_DIR""gatk.fastqIntervals.bed" | \
    Rscript -e 'summary (as.numeric (readLines ("stdin")))'
echo

if [[ ! -z ${MANTAVCF+x} ]]; then
	echo "Extract the reference intervals in Manta records that are not having matching GATK records."
	grep -F 'Manta' "$G_VS_M" | \
	    awk '{print $(NF-1)}' | \
	    awk -F ';' '{print $1}' > \
	    temp.matched.mantaIDs.txt

	grep -vf temp.matched.mantaIDs.txt "$MANTAVCF" | \
	    awk 'BEGIN {OFS="	"}; match($8, /END=[0-9]+/){print $1, $2, substr($8,RSTART+4,RLENGTH-4), $3}' > \
	    "$OUTPUT_DIR""manta_notMatchedByGATK.bed"

	echo "To find if the unmatched Manta calls have no corresponding local assemblies"
	echo "  First find the Manta records that are covering 100% of those intervals."
	bedtools intersect \
	    -a "$GATK_INTERVALS_BED" \
	    -b "$OUTPUT_DIR""manta_notMatchedByGATK.bed" \
	    -f 1 -wao | \
	    grep -F 'Manta' | awk '{print $(NF-1)}' | sort | uniq > temp.manta_notMatch_coveringFullGATKIntervals.txt
	echo "  Such intervals are suspected to be large events, and the following simple descriptive statistics checks it"
	grep -f temp.manta_notMatch_coveringFullGATKIntervals.txt "$MANTAVCF" | \
	    awk 'match($8, /SVLEN=-[0-9]+/){print substr($8,RSTART+7,RLENGTH-7)}' | \
	    Rscript -e 'summary (as.numeric (readLines ("stdin")))'

	echo "Next overlap between the two bed files:"
	echo "  [1] ""$GATK_INTERVALS_BED"
	echo "  [2] $OUTPUT_DIR""manta_notMatchedByGATK.bed"
	bedtools intersect \
	    -a "$GATK_INTERVALS_BED" \
	    -b "$OUTPUT_DIR""manta_notMatchedByGATK.bed" \
	    -wao | \
	    grep -F 'Manta' | awk '{print $(NF-1)}' | sort | uniq | \
	    grep -vf temp.manta_notMatch_coveringFullGATKIntervals.txt | wc -l

	echo "If the number of unmatched Manta calls" `(wc -l "$OUTPUT_DIR""manta_notMatchedByGATK.bed" | awk '{print $1}')` "is far greater than the number above,"
	echo "this indicates that relatively few unmatched Manta calls have corresponding intervals that have local assembly triggered."
else
	awk 'BEGIN {OFS="	";}{print $5, $6, $7, $8}' "$M_VS_P" > temp1.txt
	awk 'BEGIN {OFS="	";}{print $5, $6, $7, $8}' "$G_VS_P" > temp2.txt
	grep -f temp1.txt temp2.txt | sort | uniq > "$OUTPUT_DIR""pacbioEntries.sharedByGATKandMANTA.bed"
	wc -l "$OUTPUT_DIR""pacbioEntries.sharedByGATKandMANTA.bed"
	grep -vf temp1.txt temp2.txt | sort | uniq > "$OUTPUT_DIR""pacbioEntries.onlyInGATK.bed"
	wc -l "$OUTPUT_DIR""pacbioEntries.onlyInGATK.bed"
	grep -vf temp2.txt temp1.txt | sort | uniq > "$OUTPUT_DIR""pacbioEntries.onlyInMANTA.bed"
	wc -l "$OUTPUT_DIR""pacbioEntries.onlyInMANTA.bed"
	echo

	echo "Next overlap between the two bed files:"
	echo "  [1] ""$GATK_INTERVALS_BED"
	echo "  [2] ""$OUTPUT_DIR""pacbioEntries.onlyInMANTA.bed"
	bedtools intersect \
	    -a "$GATK_INTERVALS_BED" \
	    -b "$OUTPUT_DIR""pacbioEntries.onlyInMANTA.bed" | \
	    awk '{print $4}' | grep -Eo '[0-9]+' | sort -n | uniq \
	    > "$OUTPUT_DIR""pacbioEntries.onlyInMANTA.overlappingAssemblies.txt"
	wc -l "$OUTPUT_DIR""pacbioEntries.onlyInMANTA.overlappingAssemblies.txt" | awk '{print $1}'

	echo "If the number of unmatched intervals" `(wc -l "$OUTPUT_DIR""pacbioEntries.onlyInMANTA.bed" | awk '{print $1}')` "is far greater than the number above,"
	echo "this indicates that relatively few unmatched Manta calls have corresponding intervals that have local assembly triggered."
	echo

	echo "Also check the shared calls between GATK and Manta but not validated in PacBio."
	awk '{print $4}' "$OUTPUT_DIR""cleanDel.gatkVSmanta.intersection_f08r.txt" > temp.txt
	grep -vf temp.txt "$OUTPUT_DIR""cleanDel.gatkVSpacbio.intersection_f08r.txt" | \
	    awk '{print $4}' | awk -F ';' '{print $1}' > \
	    "$OUTPUT_DIR""variants.discoveredByGATKandMANTA.missingInPacBio.gatkID.txt"
	wc -l "$OUTPUT_DIR""variants.discoveredByGATKandMANTA.missingInPacBio.gatkID.txt"
fi

rm -f temp*

echo -e "\nDone checking GATK local assembly trigger sensitivity.\n"
