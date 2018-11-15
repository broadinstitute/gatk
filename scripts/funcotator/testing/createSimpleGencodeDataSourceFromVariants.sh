#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script will create a simple Gencode data source for all the variants 
# in a given VCF file.
# You will have to index the fasta and feature files created by this as well as
# create sequence dictionaries for them.
# It must be internally configured to point at a valid funcotatior data sources
# directory, and will not work for you out-of-the-box unless you are Jonn Smith.
#
# EXAMPLE:
#     ./createSimpleGencodeDataSourceFromVariants.sh hg19 TEST.VCF
#
# AUTHOR: Jonn Smith
#
###############################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( readlink "${BASH_SOURCE[0]}" || echo "${BASH_SOURCE[0]}" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=2
MAXARGS=99

################################################################################

# Change this path to your data sources directory:
DATA_SOURCES_PATH=/Users/jonn/Development/funcotator_dataSources_latest

GENCODE_PATH=${DATA_SOURCES_PATH}/gencode
GENCODE_HG19=${GENCODE_PATH}/hg19
GENCODE_HG38=${GENCODE_PATH}/hg38
GENCODE_HG19_GTF=${GENCODE_HG19}/gencode.v19.annotation.REORDERED.gtf
GENCODE_HG38_GTF=${GENCODE_HG38}/gencode.v28.annotation.REORDERED.gtf
GENCODE_HG19_TX=${GENCODE_HG19}/gencode.v19.pc_transcripts.fa
GENCODE_HG38_TX=${GENCODE_HG38}/gencode.v28.pc_transcripts.fa
GENCODE_HG19_VERSION=19
GENCODE_HG38_VERSION=28

OUT_FOLDER=tmp_gencode

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME REFVERSION VCFFILE [VCFFILE2 VCFFILE3 ...]"
  echo -e "Create a simple Gencode data source for all variants in "
  echo -e "each given VCFFILE."
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo -e ""
  echo -e "VCFFILE       A VCF File containing variants for which to look up the genes in Gencode."
  echo -e "REFVERSION    The reference version in which to look up variants." 
  echo -e "              Must be one of:"
  echo -e "                hg19"
  echo -e "                hg38"
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  VCF FILE DOES NOT EXIST" 
  echo -e "  4  INVALID REFERENCE VERSION" 
  echo -e ""
}


#Display a message to std error:
function error() 
{
  echo "$1" 1>&2 
}

TMPFILELIST=''
function makeTemp()
{
  local f
  f=$( mktemp )
  TMPFILELIST="${TMPFILELIST} $f"
  echo $f
}

function cleanTempVars()
{
  rm -f ${TMPFILELIST}
}

function at_exit()
{
  cleanTempVars
}

################################################################################

trap at_exit EXIT 

################################################################################

function getGeneList() {
	grep '\tgene\t' $GENCODE 
}

################################################################################

#Check given arguments:
if [[ $# -gt $MAXARGS ]] ; then
  usage
  exit 1
elif [[ $# -lt $MINARGS ]] ; then
  usage
  exit 2
fi

################################################################################

[[ -e ${OUT_FOLDER} ]] && error "Output folder already exists: ${OUT_FOLDER}" && exit 8

REF=$( echo $1 | tr '[:upper:]' '[:lower:]' )
[[ "${REF}" != "hg19" ]] && [[ "${REF}" != "hg38" ]] && error "Invalid reference version: ${REF}" && exit 4
shift

VCFFILENAMESLIST=$(makeTemp)
while [ $# -gt 0 ] ; do
	VCFFILE=$1
	[[ ! -f ${VCFFILE} ]] && error "VCF File does not exist: ${VCFFILE}" && exit 3
	echo "${VCFFILE}" >> ${VCFFILENAMESLIST}
	shift
done

################################################################################
# Make sure we point at the right gencode:
GENCODE="$GENCODE_HG19_GTF"
GENCODE_TX="$GENCODE_HG19_TX"
GENCODE_VERSION="$GENCODE_HG19_VERSION"
if [[ "${REF}" == "hg38" ]] ;then 
	GENCODE="$GENCODE_HG38_GTF"
	GENCODE_TX="$GENCODE_HG38_TX"
	GENCODE_VERSION="$GENCODE_HG38_VERSION"
fi

################################################################################
# Set up our output:
mkdir -p ${OUT_FOLDER}/${REF}
OUT_GTF_FILE=${OUT_FOLDER}/${REF}/gencode.v${GENCODE_VERSION}.testVariantSubset.gtf
OUT_TX_FILE=${OUT_FOLDER}/${REF}/gencode.v${GENCODE_VERSION}.testVariantSubset.pc_transcripts.fa
sed -e "s#src_file =.*#src_file = gencode.v${GENCODE_VERSION}.testVariantSubset.gtf#g" \
	  -e "s#gencode_fasta_path =.*#gencode_fasta_path = gencode.v${GENCODE_VERSION}.testVariantSubset.pc_transcripts.fa#g" \
		-e "s#preprocessing_script =\\(.*\\)#preprocessing_script =\\$1 , ${SCRIPTNAME}#g" \
		${GENCODE_HG19}/gencode.config > ${OUT_FOLDER}/${REF}/gencode.config

tmpFile=$( makeTemp )	
tmpGeneList=$( makeTemp )
tmpGeneIdList=$( makeTemp )

# Get our list of genes:
while read VCFFILE ; do 

	error "Processing ${VCFFILE} ..."

	# Filter out the commented out lines:
	rm -f ${tmpFile}
	grep -v '^#' ${VCFFILE} > $tmpFile

	# Go through our VCF file:
	error "  Getting gene info for variants ..."
	while read line ; do 
		CONTIG=$( echo "${line}" | awk '{print $1}' )
		if [[ $CONTIG != chr*  ]] && [[ "${REF}" == "hg19" ]] ; then 
			CONTIG="chr${CONTIG}"
		fi
		START=$( echo "${line}" | awk '{print $2}' )
		${SCRIPTDIR}/getGeneFromGenomicCoordinates.sh ${REF} ${CONTIG} ${START} ${START} 
	done < ${tmpFile} | uniq 

done < ${VCFFILENAMESLIST} | sort -k1.4n -k1,1 -k4,4n -k5,5n | uniq > ${tmpGeneList}
error "Done"

# Create a clean place to put our gene list:
rm -f ${tmpGeneIdList}
awk '{print $10}' ${tmpGeneList} | tr -d ';"' > ${tmpGeneIdList}

totalGenes=$( cat ${tmpGeneIdList} | sort | uniq | wc -l | awk '{print $1}' )
uniqueGeneIdList=$( makeTemp )

# Process the gencode fasta file:
transcriptFile=$( makeTemp )
${SCRIPTDIR}/reformatFastaSequencesToOneLine.sh ${GENCODE_TX} > ${transcriptFile}

# Get our gencode header:
head -n5 ${GENCODE} > ${OUT_GTF_FILE}

# Grep for each gene ID in our gencode file:
i=0
error "Retrieving genes and sequences from GENCODE ..."
while read geneId ; do 

	# make sure we only add each gene once:
	grep -m1 "${geneId}" ${uniqueGeneIdList} &> /dev/null
	r=$?
	[[ $r -eq 0 ]] && continue

	# Arcane awk statement adapted from https://backreference.org/2010/09/11/smart-ranges-in-awk/
	# and https://backreference.org/2010/09/11/smart-ranges-in-awk/
	# Essentially uses a flag to determine whether to print matches and when the flag is no longer
	# true it will exit.
	BANG=!
	awk "${BANG}/${geneId}/{if (p == 1){exit} else{p=0}} /${geneId}/{p=1} p" ${GENCODE} >> ${OUT_GTF_FILE}

	# Get the sequence here:
	grep -A1 "${geneId}" ${transcriptFile} | grep -v '^--' >> ${OUT_TX_FILE}
	
	let i=$i+1
	error "  Processed gene $i of $totalGenes ($geneId)"
	echo "${geneId}" >> ${uniqueGeneIdList}

done < ${tmpGeneIdList}
error "Done"

error "You must index the created feature and fasta files, as well as create a sequence dictionary for the output fasta."

