#!/bin/bash

set -eu

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "BEGIN TO WORK ON COMPLEX VARIANTS"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'

if [[ ! -z "$GATK_CPX_VCF_MASTER" ]]; then
	ANALYSIS_DIR_MASTER="$ANALYSIS_DIR_MASTER/CPX"
	mkdir -p "$ANALYSIS_DIR_MASTER"
	bash "$SCRIPT_DIR""/parseCpx.sh" \
		"$GATK_CPX_VCF_MASTER" \
		"$PACBIO_VCF_CHM1" \
		"$PACBIO_VCF_CHM13" \
		"$ANALYSIS_DIR_MASTER"
fi

if [[ ! -z "$GATK_CPX_VCF_FEATURE" ]]; then
	ANALYSIS_DIR_FEATURE="$ANALYSIS_DIR_FEATURE/CPX"
	mkdir -p "$ANALYSIS_DIR_FEATURE"
	bash "$SCRIPT_DIR""/parseCpx.sh" \
		"$GATK_CPX_VCF_FEATURE" \
		"$PACBIO_VCF_CHM1" \
		"$PACBIO_VCF_CHM13" \
		"$ANALYSIS_DIR_FEATURE"
fi

echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'
echo "DONE WORKING ON COMPLEX VARIANTS"
echo -e '\033[0;35m#################################################\033[0m'
echo -e '\033[0;35m#################################################\033[0m'