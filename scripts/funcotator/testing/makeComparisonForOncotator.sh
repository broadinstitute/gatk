#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script processes files created by Funcotator and Oncotator to create 
# output files that can be diffed by Beyond Compare (or your favorite diff 
# viewer).
# It must be internally configured to point at valid Funcotator and Oncotator
# output files, and will not work for you out-of-the-box unless you are Jonn Smith.
#
# EXAMPLE:
#     ./makeComparisonForOncotator.sh
#
# AUTHOR: Jonn Smith
#
###############################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( readlink "${BASH_SOURCE[0]}" || echo "${BASH_SOURCE[0]}" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=1

################################################################################

# Change these files to point to TSV output from oncotator (onk), funcotator
# (funk), and your Gencode GTF file (GENCODE)
onk='/Users/jonn/Development/oncotator_testing/test.maf.FUNCOTATOR_VCLASS.tsv'
funk='FUNCOTATOR_OUT.maf.tsv'
 
GENCODE="/Users/jonn/Development/funcotator_dataSources_latest/gencode/hg19/gencode.v19.annotation.REORDERED.gtf"

################################################################################

CACHEDTRANSCRIPTS=cachedTranscriptInfo.txt

################################################################################

function error() {
	echo $@ 1>&2
}

function crossReferenceWithGencode() {

	fname=$1
	firstLine=true
	
	while read line ; do 
		if $firstLine ; then
			echo -e "${line}\tVARIANT_INSIDE_TX"
			firstLine=false
			continue
		fi

		variantContig=$( echo $line | awk '{print $2}' )
		variantStart=$( echo $line | awk '{print $3}' )
		variantEnd=$( echo $line | awk '{print $4}' )
		VC=$( echo $line | awk '{print $5}' )
		tx=$( echo $line | awk '{print $11}' )	
		
		if [[ "${VC}" == "IGR" ]] ; then
			echo -e "${line}\tFalse"
			error "Ignoring IGR ${variantContig}:${variantStart}"
			continue
		elif [[ ${#tx} -eq 0 ]] ; then
			echo -e "${line}\tFalse"
			error "Empty TX for ${variantContig}:${variantStart}" 
			continue
		fi

		error -n "Looking up $tx in GENCODE..." 

		local cached=false
		if [[ -f $CACHEDTRANSCRIPTS ]] ; then 
			gencodeLine=$( grep -m1 "$tx" $CACHEDTRANSCRIPTS )
			[[ ${#gencodeLine} -eq 0 ]] && gencodeLine=$( grep -m1 "$tx" $GENCODE ) && echo "$gencodeLine" >> $CACHEDTRANSCRIPTS && cached=true
		else 
			gencodeLine=$( grep -m1 "$tx" $GENCODE )
			echo "$gencodeLine" >> $CACHEDTRANSCRIPTS
		fi
		tx_gene=$( echo "${gencodeLine}" | awk '{print $18}' | sed -e 's#"##g' -e 's#;##g' )
		tx_contig=$( echo "${gencodeLine}" | awk '{print $1}')
		tx_start=$( echo "${gencodeLine}" | awk '{print $4}')
		tx_end=$( echo "${gencodeLine}" | awk '{print $5}')
		
		if $cached ; then 
			#error " Found (CACHE) at [${tx_gene}] - $tx_contig : $tx_start -> $tx_end"
			error " Found (CACHE)"
		else
			error " Found" 
		fi

		insideTranscript="False"
		if [[ $tx_start -le $variantStart ]] && [[ $variantStart -le $tx_end ]] && [[ $tx_start -le $variantEnd ]] && [[ $variantEnd -le $tx_end ]] ; then
			insideTranscript="True"
		fi

		echo -e "${line}\t${insideTranscript}"
	done < $fname
}

################################################################################

${SCRIPTDIR}/createFuncotationCsvFromFuncotatorVcf.sh FUNCOTATOR_OUT.vcf > FUNCOTATOR_OUT.csv
cut -d, -f1,3,4,5,6,8,9,10,11,12,13,14,17,18,19,21,22 FUNCOTATOR_OUT.csv | sed -e 's#,chr#,#g' | tr ',' '\t' > FUNCOTATOR_OUT.maf.tsv

# Find all the lines that are different between the two files:
# To make our lives easier.
grep -F -x -v -f $funk $onk | grep -v '^#' > onk.tsv
grep -F -x -v -f $onk $funk | grep -v '^#' > funk.tsv

# Create versions without the `other transcripts` field:
rev onk.tsv  | cut -d '	' -f 2- | rev > onk.no_others.tsv
rev funk.tsv | cut -d '	' -f 2- | rev > funk.no_others.tsv

# Create versions that are cross-referenced with their transcripts in Gencode.
# This will show whether the variants actually fall INSIDE the transcripts
# with which they are associated:

crossReferenceWithGencode onk.no_others.tsv  > onk.xrefgencode.tsv
crossReferenceWithGencode funk.no_others.tsv > funk.xrefgencode.tsv

