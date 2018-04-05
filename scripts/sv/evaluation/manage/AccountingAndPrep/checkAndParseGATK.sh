#!/bin/bash

set -eu

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "BEGIN TO CHECK AND SPLIT GATK VCF FILE"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'

echo "Diverts different types of variants in VCF files"
echo "  and collects & saves variants sizes by type"
echo

####################
if [[ ! "$SCRIPT_DIR" =~ .+/$ ]]; then
        SCRIPT_DIR+="/"
fi

GATK_SCRIPT="$SCRIPT_DIR""AccountingAndPrep/checkSplitAndCollectSizes_gatk.sh"

##########
bash "$GATK_SCRIPT" \
     "$GATK_VCF_MASTER"  \
     "$ANALYSIS_DIR_MASTER" \
     "$GATK_CPX_DERIVED_ONE_SEG_VCF_MASTER" \
     "$GATK_CPX_DERIVED_MULTI_SEG_VCF_MASTER" \
     "$REF_VER"

bash "$GATK_SCRIPT" \
     "$GATK_VCF_FEATURE"  \
     "$ANALYSIS_DIR_FEATURE" \
     "$GATK_CPX_DERIVED_ONE_SEG_VCF_FEATURE" \
     "$GATK_CPX_DERIVED_MULTI_SEG_VCF_FEATURE" \
     "$REF_VER"

echo "#################################################"
echo " Delegating to R for producing plots"
echo 

#CURR_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#Rscript "$CURR_DIR/plotInsDelSizes.R" "$ANALYSIS_DIR_MASTER" "$MANTA_AND_PACBIO_PREANALYSIS_DIR" "$SCRIPT_DIR"
#Rscript "$CURR_DIR/plotInsDelSizes.R" "$ANALYSIS_DIR_FEATURE" "$MANTA_AND_PACBIO_PREANALYSIS_DIR" "$SCRIPT_DIR"

echo 
echo " Done producing plots"
echo "#################################################"

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE CHECKING AND SPLITTING GATK VCF FILE"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'