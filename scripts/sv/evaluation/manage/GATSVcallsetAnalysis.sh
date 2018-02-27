#!/bin/bash 

## If this script is to be appended after the pipeline in manage_sv_pipeline.sh:
##   1) the sample must be CHM mix
##   2) it should happen before results are copied to a GCS bucket.
##   3) analysis results should be copied to the bucket as well.

# # required--copy from a provided bucket
# pacbioVCF1
# pacbioVCF13
# pacbioVCFmix

# # required--copy from hdfs
# gsvcf

# # required--create on the fly
# pbdir
# gsvdir

# # optional--copy from hdfs
# gsvcpxf

# # optional--copy from a provided bucket
# mantaVCF
# gsvf2
# gsvcpxf2

# # optional--create on the fly
# gsvdir2

## If this script is to be run separately, then the following is an example
# ./GATSVcallsetAnalysis.sh \
# --pacbioVCF1 ~/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Huddleston2016/CHM1_SVs.annotated.vcf \
# --pacbioVCF13 ~/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Huddleston2016/CHM13_SVs.annotated.vcf \
# --pacbioVCFmix ~/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Huddleston2016/CHM1_CHM13_pseudodiploid_SVs.vcf \
# --gsvcf ~/Project/SV/Discovery/GATK/G94794.CHMI_CHMI3_WGS1/expInterpretNonComplex.vcf.gz \
# --pbdir ~/Project/SV/Discovery/OtherCallsets/pseudoDiploid/PreAnalysis \
# --gsvdir ~/Project/SV/Discovery/Evaluation/Analysis/Feature \
# --mantaVCF ~/Project/SV/Discovery/OtherCallsets/pseudoDiploid/Manta103/G94794.CHMI_CHMI3_WGS1/Manta103_G94794-CHMI_CHMI3_WGS1-cram-bam.vcf \
# --gsvf2 ~/Project/SV/Discovery/GATK/G94794.CHMI_CHMI3_WGS1/inv_del_ins.vcf.gz \
# --gsvcpxf ~/Project/SV/Discovery/GATK/G94794.CHMI_CHMI3_WGS1/expInterpretCpx.vcf.gz \
# --gsvdir2 ~/Project/SV/Discovery/Evaluation/Analysis/Master

set -eu

usage() 
{ 
  cat << EOF
Performs simple-overlap-based analysis on GATK-SV pipeline call set against PacBio callsets on benchmarking CHM samples.
  Options (mandatory): 
EOF
  cat <<EOF | column -s\& -t 
    --pacbioVCF1 & path to PacBio VCF on CHM-1 haploid sample
    --pacbioVCF13 & path to PacBio VCF on CHM-13 haploid sample
    --pacbioVCFmix & path to PacBio VCF on CHM-1 and CHM-13 mixture sample
    --gsvcf & path to GATK VCF on non-complex variants
    --pbdir & local directory to save analysis results on PacBio call sets
    --gsvdir & local directory to save analysis results on GATK variant calls
EOF
  cat << EOF
  Options (optional): 
EOF
  cat <<EOF | column -s\& -t 
    -h|--help & show this output 
    --mantaVCF & path to Manta VCF on CHM-1 and CHM-13 mixture sample
    --gsvf2 & path to another GATK VCF on non-complex variants (usually the call set from master)
    --gsvcpxf & path to GATK VCF on complex variants
    --gsvcpxf2 & path to another GATK VCF on complex variants (usually the call set from master)
    --gsvdir2 & local directory to save analysis results on the other GATK variant calls
    --masterVSfeatureDir & local directory to save analysis results on comparing TP/FP between the two GATK call sets
EOF
} 

throw_error() {
    echo "$1" >&2
    exit 1
}

#####
if [[ $# -eq 0 ]]; then
    usage
    exit 0;
fi

while [ $# -ge 1 ]; do
    case $1 in
        -h|-\?|--help)
            usage
            exit 0
            ;;
        --pacbioVCF1)
            if [ $# -ge 2 ]; then
                PACBIO_VCF_CHM1="$2"
                shift 2
            else
                throw_error "--pacbioVCF1 requires a non-empty argument"
            fi
            ;;
        --pacbioVCF1=?*)
            PACBIO_VCF_CHM1=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --pacbioVCF13)
            if [ $# -ge 2 ]; then
                PACBIO_VCF_CHM13="$2"
                shift 2
            else
                throw_error "--pacbioVCF13 requires a non-empty argument"
            fi
            ;;
        --pacbioVCF13=?*)
            PACBIO_VCF_CHM13=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --pacbioVCFmix)
            if [ $# -ge 2 ]; then
                PACBIO_VCF_CHMMIX="$2"
                shift 2
            else
                throw_error "--pacbioVCFmix requires a non-empty argument"
            fi
            ;;
        --pacbioVCFmix=?*)
            PACBIO_VCF_CHMMIX=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --gsvcf)
            if [ $# -ge 2 ]; then
                GATK_VCF_FEATURE="$2"
                shift 2
            else
                throw_error "--gsvcf requires a non-empty argument"
            fi
            ;;
        --gsvcf=?*)
            GATK_VCF_FEATURE=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --pbdir)
            if [ $# -ge 2 ]; then
                MANTA_AND_PACBIO_PREANALYSIS_DIR="$2"
                shift 2
            else
                throw_error "--pbdir requires a non-empty argument"
            fi
            ;;
        --pbdir=?*)
            MANTA_AND_PACBIO_PREANALYSIS_DIR=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --gsvdir)
            if [ $# -ge 2 ]; then
                ANALYSIS_DIR_FEATURE="$2"
                shift 2
            else
                throw_error "--gsvdir requires a non-empty argument"
            fi
            ;;
        --gsvdir=?*)
            ANALYSIS_DIR_FEATURE=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --gsvcpxf)
            if [ $# -ge 2 ]; then
                GATK_CPX_VCF_FEATURE="$2"
                shift 2
            else
                throw_error "--gsvcpxf requires a non-empty argument"
            fi
            ;;
        --gsvcpxf=?*)
            GATK_CPX_VCF_FEATURE=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --mantaVCF)
            if [ $# -ge 2 ]; then
                MANTA_VCF="$2"
                shift 2
            else
                throw_error "--mantaVCF requires a non-empty argument"
            fi
            ;;
        --mantaVCF=?*)
            MANTA_VCF=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --gsvf2)
            if [ $# -ge 2 ]; then
                GATK_VCF_MASTER="$2"
                shift 2
            else
                throw_error "--gsvf2 requires a non-empty argument"
            fi
            ;;
        --gsvf2=?*)
            GATK_VCF_MASTER=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --gsvcpxf2)
            if [ $# -ge 2 ]; then
                GATK_CPX_VCF_MASTER="$2"
                shift 2
            else
                throw_error "--gsvcpxf2 requires a non-empty argument"
            fi
            ;;
        --gsvcpxf2=?*)
            GATK_CPX_VCF_MASTER=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --gsvdir2)
            if [ $# -ge 2 ]; then
                ANALYSIS_DIR_MASTER="$2"
                shift 2
            else
                throw_error "--gsvdir2 requires a non-empty argument"
            fi
            ;;
        --gsvdir2=?*)
            ANALYSIS_DIR_MASTER=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --masterVSfeatureDir)
            if [ $# -ge 2 ]; then
                ANALYSIS_DIR_MASTER_VS_FEATURE="$2"
                shift 2
            else
                throw_error "--masterVSfeatureDir requires a non-empty argument"
            fi
            ;;
        --masterVSfeatureDir=?*)
            ANALYSIS_DIR_MASTER_VS_FEATURE=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        --)   # explicit call to end of all options
            shift
            break
            ;;
        -?*)  # unsupported option
            throw_error "Unknown option \"$1\". use --help for syntax"
            ;;
        *)  # not an option, a positional argument. break out
            break
            ;;
    esac
done

if [[ ! -z $ANALYSIS_DIR_MASTER_VS_FEATURE ]]; then
    if [[ -z $GATK_VCF_MASTER ]]; then
        throw_error "local directory to save analysis results on comparing TP/FP between the two GATK call sets is provided, but only one GATK call set is given"
    fi
fi
# ###################

# export PACBIO_VCF_CHM1
# export PACBIO_VCF_CHM13
# export PACBIO_VCF_CHMMIX
# export MANTA_AND_PACBIO_PREANALYSIS_DIR
# export GATK_VCF_MASTER
# export GATK_VCF_FEATURE
# export ANALYSIS_DIR_MASTER
# export ANALYSIS_DIR_FEATURE
# export GATK_CPX_VCF_MASTER
# export GATK_CPX_VCF_FEATURE
# export MANTA_VCF
# export ANALYSIS_DIR_MASTER_VS_FEATURE

# export REF_VER=38

# export SCRIPT_DIR="../general/"
# source "$SCRIPT_DIR""misc_functions.sh"

# mkdir -p "$ANALYSIS_DIR_FEATURE" "$ANALYSIS_DIR_MASTER" "$MANTA_AND_PACBIO_PREANALYSIS_DIR" "$ANALYSIS_DIR_MASTER_VS_FEATURE"

# ################### run pre-analysis on Manta and PacBio callsets
# bash AccountingAndPrep/checkAndParsePacBioAndManta.sh || on_error

# ################### then run analysis size and divert variants by types
# bash AccountingAndPrep/checkAndParseGATK.sh || on_error

# ################### then run analysis on deletion calls
# bash Deletion/driver.sh || on_error

# ################### then run analysis on insertion calls
# bash InsDupAndRPL/driver.sh || on_error

# ###################
# bash CpxBND/driver.sh || on_error


# ###################
# # optionally run analysis between GATK callsets
# if [[ ! -z $ANALYSIS_DIR_MASTER_VS_FEATURE ]]; then
#     bash masterVSfeature.sh || on_error
# fi

# echo -e '\033[0;35m#################################################\033[0m'
# echo -e '\033[0;35m#################################################\033[0m'
# echo "DONE COLLECTING INFORMATION AND SIMPLE TP/FP ANALYSIS"
# echo "PLEASE REMEMBER TO:"
# echo " 1. manually review the \"singletons.bed\" generated in \"masterVSfeature\" with the igv session (requires internet connection)"
# echo " 2. further analyze variants that are filtered out as \"duplicated\" from both master and feature callsets at an early stage of the pipeline"
# echo "ALSO KEEP IN MIND:"
# echo "This set of script currently performs analysis with simple intersetion analysis for"
# echo "   [DEL]: 50% reciprocal overlap with \"bedtools intersect\" on variants that pass fileter MQ=60 and ALN=50"
# echo "   [INS/RPL]: simple intersection with \"bedtools window -w 50\" (i.e. padding left/right 50 bases) and NO filtering on variants"
# echo "   [DUP]: simple intersection with \"bedtools window -l 0 -r 50\" (i.e. padding left/right 0/50 bases) and NO filtering on variants"
# echo "   [CPX]: simple intersection with \"bedtools intersect\" and NO filtering on variants"
# echo -e '\033[0;35m#################################################\033[0m'
# echo -e '\033[0;35m#################################################\033[0m'