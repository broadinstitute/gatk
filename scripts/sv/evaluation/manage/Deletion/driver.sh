#!/bin/bash

set -eu

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "BEGIN TO WORK ON DELETION CASES"
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

ANALYSIS_DIR_MASTER="$ANALYSIS_DIR_MASTER""Deletion/"
ANALYSIS_DIR_FEATURE="$ANALYSIS_DIR_FEATURE""Deletion/"
SCRIPT_DIR="$SCRIPT_DIR""Deletion"
MANTA_AND_PACBIO_PREANALYSIS_DIR="$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion/"

echo "#################################################"
echo "Here we make three BED4 files for clean deletions discoverred by the GATK records"
echo "  where the 4-th column is ID;MAX_ALIGN_LENGTH;HOMLEN"
echo "#################################################"
echo

##### master
bash "$SCRIPT_DIR/parseGATK.sh" \
     "$ANALYSIS_DIR_MASTER"GATK_primaryContigs_cleanDel.txt \
     "$ANALYSIS_DIR_MASTER"

##### feature
bash "$SCRIPT_DIR/parseGATK.sh" \
     "$ANALYSIS_DIR_FEATURE"GATK_primaryContigs_cleanDel.txt \
     "$ANALYSIS_DIR_FEATURE" \
     60 \
     50

echo "#################################################"
echo "Now we run simple intersection analsis: gatk master callset"
echo "#################################################"
echo

bash "$SCRIPT_DIR/simpleIntersection.sh" \
     "$ANALYSIS_DIR_MASTER""gatk.cleanDel.bed" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""manta.cleanDel.bed" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""pacbio.cleanDel_1.bed" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""pacbio.cleanDel_13.bed" \
     "$ANALYSIS_DIR_MASTER" \
     0.5

echo "#################################################"
echo "Now we run simple intersection analsis: gatk feature callset"
echo "#################################################"
echo

bash "$SCRIPT_DIR/simpleIntersection.sh" \
     "$ANALYSIS_DIR_FEATURE""gatk.cleanDel.bed" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""manta.cleanDel.bed" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""pacbio.cleanDel_1.bed" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR""pacbio.cleanDel_13.bed" \
     "$ANALYSIS_DIR_FEATURE" \
     0.5

# echo "#################################################"
# echo "Here we do a simple check and see if the variants"
# echo "  missing in GATK are caused by lack of assembly"
# echo "#################################################"
# echo
# if [[ -z ${INTERVAL_FILE+x} ]]; then
#      bash "$SCRIPT_DIR/checkGATKAssemblyTriggerSensitivity.sh" \
#           "$INTERVAL_FILE" \
#           "$ANALYSIS_DIR_MASTER" \
#           "$ANALYSIS_DIR_MASTER"Results/cleanDel.gatkVSmanta.intersection_f08r.txt \
#           "$ANALYSIS_DIR_MASTER"Results/cleanDel.gatkVSpacbio.intersection_f08r.txt \
#           "$ANALYSIS_DIR_MASTER"Results/cleanDel.mantaVSpacbio.intersection_f08r.txt
# fi


# Rscript naiveAssociation.R "$ANALYSIS_DIR_MASTER"

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE WORKING ON DELETION CASES"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'