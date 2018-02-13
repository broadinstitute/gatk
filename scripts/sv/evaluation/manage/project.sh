#!/bin/bash

set -eu

# 1. prep step needs
#   * (optional) path to script for analyzing Manta and PacBio call sets
#   * (optional) path to Manta and PacBio callsets
#   * (optional) path to save analysis results on Manta and PacBio call sets
#   * path to analysis results on Manta and PacBio callsets
#   * path to script for analyzing GATK script
#   * path to GATK master and feature callsets
#   * path to save analysis results on master and feature callsets
# 2. deletion step needs
#   * path to analysis results on Manta and PacBio callsets
#   * path to script for analyzing deletion calls
#   * path to save analysis results on deletion calls for master and feature callsets
# 3. insertion step needs
#   * path to analysis results on Manta and PacBio callsets
#   * path to script for analyzing deletion calls
#   * path to save analysis results on insertion calls for master and feature callsets
# 4. cpx step needs
#   * path to PacBio callsets on the two CHM haploid celllines
#   * path to GATK cpx variant call sets from both master and feature
#   * path to save analysis results on cpx calls for master and feature callsets


export REF_VER=38

export SCRIPT_DIR="../general/"
source "$SCRIPT_DIR""misc_functions.sh"


export PACBIO_VCF_CHM1="/Users/shuang/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Huddleston2016/CHM1_SVs.annotated.vcf"
export PACBIO_VCF_CHM13="/Users/shuang/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Huddleston2016/CHM13_SVs.annotated.vcf"
export PACBIO_VCF_CHMMIX="/Users/shuang/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Huddleston2016/CHM1_CHM13_pseudodiploid_SVs.vcf"
export MANTA_AND_PACBIO_PREANALYSIS_DIR="/Users/shuang/Project/SV/Discovery/OtherCallsets/pseudoDiploid/PreAnalysis/"

export GATK_VCF_MASTER="/Users/shuang/Project/SV/Discovery/GATK/G94794.CHMI_CHMI3_WGS1/inv_del_ins.vcf.gz"
export GATK_VCF_FEATURE="/Users/shuang/Project/SV/Discovery/GATK/G94794.CHMI_CHMI3_WGS1/expInterpretNonComplex.vcf.gz"

export ANALYSIS_DIR_MASTER="/Users/shuang/Project/SV/Discovery/Evaluation/Analysis/Master/"
export ANALYSIS_DIR_FEATURE="/Users/shuang/Project/SV/Discovery/Evaluation/Analysis/Feature/"
rm -rf "$ANALYSIS_DIR_MASTER" && mkdir "$ANALYSIS_DIR_MASTER"
rm -rf "$ANALYSIS_DIR_FEATURE" && mkdir "$ANALYSIS_DIR_FEATURE"

export GATK_CPX_VCF_MASTER=""
export GATK_CPX_VCF_FEATURE="/Users/shuang/Project/SV/Discovery/GATK/G94794.CHMI_CHMI3_WGS1/expInterpretCpx.vcf.gz"

################### conditioanlly run pre-analysis on Manta and PacBio callsets
RUN_ANALYSIS_ON_MANTA_AND_PACIO=false
read -p "run prep analysis on Manta and PacBio callsets?" yn
case $yn in
    [Yy]*)  RUN_ANALYSIS_ON_MANTA_AND_PACIO=true
            ;;
    [Nn]*)  ;;
    [Cc]*)  exit
            ;;
    *)      echo "Please answer yes, no, or cancel."
            ;;
esac
if [[ $RUN_ANALYSIS_ON_MANTA_AND_PACIO == true ]]; then
    export MANTA_VCF="/Users/shuang/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Manta103/G94794.CHMI_CHMI3_WGS1/Manta103_G94794-CHMI_CHMI3_WGS1-cram-bam.vcf"
    rm -rf "$MANTA_AND_PACBIO_PREANALYSIS_DIR" && mkdir "$MANTA_AND_PACBIO_PREANALYSIS_DIR"
    bash AccountingAndPrep/checkAndParsePacBioAndManta.sh || on_error
fi

################### then run analysis size and divert variants by types
bash AccountingAndPrep/checkAndParseGATK.sh || on_error

################### then run analysis on deletion calls
bash Deletion/driver.sh || on_error

################### then run analysis on insertion calls
bash InsDupAndRPL/driver.sh || on_error

###################
bash CpxBND/driver.sh || on_error


###################
bash masterVSfeature.sh

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE COLLECTING INFORMATION AND SIMPLE TP/FP ANALYSIS"
echo "PLEASE REMEMBER TO:"
echo " 1. manually review the \"singletons.bed\" generated in \"masterVSfeature\" with the igv session (requires internet connection)"
echo " 2. further analyze variants that are filtered out as \"duplicated\" from both master and feature callsets at an early stage of the pipeline"
echo "ALSO KEEP IN MIND:"
echo "This set of script currently performs analysis with simple intersetion analysis for"
echo "   [DEL]: 50% reciprocal overlap with \"bedtools intersect\" on variants that pass fileter MQ=60 and ALN=50"
echo "   [INS/RPL]: simple intersection with \"bedtools window -w 50\" (i.e. padding left/right 50 bases) and NO filtering on variants"
echo "   [DUP]: simple intersection with \"bedtools window -l 0 -r 50\" (i.e. padding left/right 0/50 bases) and NO filtering on variants"
echo "   [CPX]: simple intersection with \"bedtools intersect\" and NO filtering on variants"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'