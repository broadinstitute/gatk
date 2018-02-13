#!/bin/bash

set -eu

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "BEGIN TO WORK ON INSERTION CASES"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'

command -v bedtools >/dev/null 2>&1 || { echo >&2 "I require bedtools for interval manipulation but it's not installed. Aborting."; exit 1; }

if [[ ! "$ANALYSIS_DIR_MASTER" =~ .+/$ ]]; then
        ANALYSIS_DIR_MASTER+="/"
fi
if [[ ! "$ANALYSIS_DIR_FEATURE" =~ .+/$ ]]; then
        ANALYSIS_DIR_FEATURE+="/"
fi
if [[ ! "$SCRIPT_DIR" =~ .+/$ ]]; then
        SCRIPT_DIR+="/"
fi

if [[ ! "$MANTA_AND_PACBIO_PREANALYSIS_DIR" =~ .+/$ ]]; then
        MANTA_AND_PACBIO_PREANALYSIS_DIR+="/"
fi

ANALYSIS_DIR_MASTER="$ANALYSIS_DIR_MASTER""InsDupRPL/"
ANALYSIS_DIR_FEATURE="$ANALYSIS_DIR_FEATURE""InsDupRPL/"
SCRIPT_DIR="$SCRIPT_DIR""InsDupRPL/"

echo "#################################################"
echo "Overlapping insertions and tandem duplications"
echo "#################################################"
echo

bash "$SCRIPT_DIR""againstPacBio.sh" \
	 "MANTA" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/Manta_PASS_PRECISE_nonBND_primaryContigs_cleanIns.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/Manta_PASS_PRECISE_nonBND_primaryContigs_tandup.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/Manta_PASS_PRECISE_nonBND_primaryContigs_longRangeSubstitution.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/PacBio_primaryContigs_ins_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/PacBio_primaryContigs_ins_13.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion/PacBio_primaryContigs_cleanDel_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion/PacBio_primaryContigs_cleanDel_13.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/"

bash "$SCRIPT_DIR""againstPacBio.sh" \
	 "GATK" \
     "$ANALYSIS_DIR_MASTER""GATK_primaryContigs_ins.txt" \
     "$ANALYSIS_DIR_MASTER""GATK_primaryContigs_dup.txt" \
	 "$ANALYSIS_DIR_MASTER""GATK_primaryContigs_scarredDel.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/PacBio_primaryContigs_ins_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/PacBio_primaryContigs_ins_13.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion/PacBio_primaryContigs_cleanDel_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion/PacBio_primaryContigs_cleanDel_13.txt" \
     "$ANALYSIS_DIR_MASTER"


bash "$SCRIPT_DIR""againstPacBio.sh" \
	 "GATK" \
     "$ANALYSIS_DIR_FEATURE""GATK_primaryContigs_ins.txt" \
     "$ANALYSIS_DIR_FEATURE""GATK_primaryContigs_dup.txt" \
	 "$ANALYSIS_DIR_FEATURE""GATK_primaryContigs_scarredDel.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/PacBio_primaryContigs_ins_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""InsDupRPL/PacBio_primaryContigs_ins_13.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion/PacBio_primaryContigs_cleanDel_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion/PacBio_primaryContigs_cleanDel_13.txt" \
     "$ANALYSIS_DIR_FEATURE"

# bash "$SCRIPT_DIR/sensitivityAndSpecificity.sh" \
# 	 "$DIR"GATK_primaryContigs_ins_withoutDuplicates.vcf \
# 	 "$DIR"GATK_primaryContigs_dup_withoutDuplicates.vcf \
# 	 "$DIR"Manta_PASS_PRECISE_nonBND_primaryContigs_cleanIns.vcf \
# 	 "$DIR"Manta_PASS_PRECISE_nonBND_primaryContigs_tandup_below10000.vcf \
# 	 "$DIR" \
# 	 "$INTERVAL_BED_FILE" \
# 	 "Pacbio" \
# 	 "$DIR"Results/combined.GATKvsPacbio.uniqPacbioIDs.txt \
# 	 "$DIR"Results/combined.MantavsPacbio.uniqPacbioIDs.txt \
# 	 "$DIR"Results/combined.GATKvsManta.uniqGATKids.txt \
# 	 "$DIR"Results/combined.GATKvsPacbio.uniqGATKids.txt

# Rscript concordanceOnOneHits.R

# cat differences_insWithIns.txt > differences_combined.txt
# cat differences_insWithDup.txt | grep -vF 'gatkID' >> differences_combined.txt
# cat differences_dupWithIns.txt | grep -vF 'gatkID' >> differences_combined.txt
# cat differences_dupWithDup.txt | grep -vF 'gatkID' >> differences_combined.txt
# echo -e "\nDefine a quantity to describe how much two haplotypes differ as:"
# echo -e " 2*editDist/(haptype_1_len + haptype_2_len)"
# echo -e "Here we have a summary of the above defined quantity for the overlapping records"
# sed '1d' differences_combined.txt | awk '{print 2*$2/($3+$4)}' | \
#     Rscript -e 'summary (as.numeric (readLines ("stdin")))'

# Rscript naiveAssociation.R
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE WORKING ON INSERTION CASES"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'