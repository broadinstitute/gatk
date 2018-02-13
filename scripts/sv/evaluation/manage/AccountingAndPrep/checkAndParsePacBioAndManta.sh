#!/bin/bash

set -eu

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "BEGIN TO CHECK AND SPLIT MANTA AND PACBIO VCF FILES"
echo "WARNING: since Manta and PacBio callsets are unlikely to change between runs of analysis"
echo "  this script needs only to be run once for each version of their callsets."
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'

echo "Diverts different types of variants in VCF files"
echo "  and collects & saves variants sizes by type"
echo

####################

if [[ ! "$SCRIPT_DIR" =~ .+/$ ]]; then
        SCRIPT_DIR+="/"
fi

if [[ ! "$MANTA_AND_PACBIO_PREANALYSIS_DIR" =~ .+/$ ]]; then
        MANTA_AND_PACBIO_PREANALYSIS_DIR+="/"
fi

MANTA_GENERAL_SCRIPT="$SCRIPT_DIR""AccountingAndPrep/checkSplitAndCollectSizes_manta.sh"

PACBIO_GENERAL_SCRIPT="$SCRIPT_DIR""AccountingAndPrep/checkSplitAndCollectSizes_pacbio.sh"

####################
bash "$MANTA_GENERAL_SCRIPT" \
     "$MANTA_VCF" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR" \
     "$REF_VER"

####################
bash "$PACBIO_GENERAL_SCRIPT" \
     "$PACBIO_VCF_CHMMIX" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR" \
     0 \
     "$REF_VER"

####################
bash "$PACBIO_GENERAL_SCRIPT" \
     "$PACBIO_VCF_CHM1" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR" \
     1 \
     "$REF_VER"


####################
bash "$PACBIO_GENERAL_SCRIPT" \
     "$PACBIO_VCF_CHM13" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR" \
     13 \
     "$REF_VER"



echo "#################################################"
echo "Here we make three BED4 files for clean deletions discoverred by the two callers"
echo "  the 4-th column (the \"ID column\") is formatted, for Manta records, as ID;QUAL;GT;GQ;HOMLEN,"
echo "  whereas for PacBio records, the 4-th column is \"PB\";QUAL;CONTIG_DEPTH;CONTIG_SUPPORT;REPEAT_TYPE"
echo "#################################################"
echo

SCRIPT_DIR+="Deletion"
MANTA_AND_PACBIO_PREANALYSIS_DIR="$MANTA_AND_PACBIO_PREANALYSIS_DIR""Deletion"

bash "$SCRIPT_DIR/checkMantaVCF.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/Manta_PASS_PRECISE_nonBND_primaryContigs_cleanDel.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR"
bash "$SCRIPT_DIR/parseManta.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/Manta_PASS_PRECISE_nonBND_primaryContigs_cleanDel.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR"

bash "$SCRIPT_DIR/checkPacBioVCF.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/PacBio_primaryContigs_cleanDel_0.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR"
bash "$SCRIPT_DIR/parsePacbio.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/PacBio_primaryContigs_cleanDel_0.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR" \
     0

bash "$SCRIPT_DIR/checkPacBioVCF.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/PacBio_primaryContigs_cleanDel_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR"
bash "$SCRIPT_DIR/parsePacbio.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/PacBio_primaryContigs_cleanDel_1.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR" \
     1

bash "$SCRIPT_DIR/checkPacBioVCF.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/PacBio_primaryContigs_cleanDel_13.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR"
bash "$SCRIPT_DIR/parsePacbio.sh" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR/PacBio_primaryContigs_cleanDel_13.txt" \
     "$MANTA_AND_PACBIO_PREANALYSIS_DIR" \
     13

####################
echo "Results are saved to:"
echo "  $MANTA_AND_PACBIO_PREANALYSIS_DIR"
echo 


echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE CHECKING AND SPLITTING MANTA AND PACBIO VCF FILES"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'