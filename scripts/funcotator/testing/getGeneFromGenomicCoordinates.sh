#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script will give you the overlapping encode gene entries for a given 
# reference and genomic position.
# It must be internally configured to point at a valid funcotatior data sources
# directory, and will not work for you out-of-the-box unless you are Jonn Smith.
#
# EXAMPLE:
#     ./getGeneFromGenomicCoordinates.sh hg19 1 38 1789
#
# AUTHOR: Jonn Smith
#
################################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( readlink "${BASH_SOURCE[0]}" || echo "${BASH_SOURCE[0]}" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=4
MAXARGS=4

################################################################################

# Change this path to your data sources directory:
DATA_SOURCES_PATH=/Users/jonn/Development/funcotator_dataSources_latest

GENCODE_HG19=${DATA_SOURCES_PATH}/gencode/hg19/gencode.v19.annotation.REORDERED.gtf
GENCODE_HG38=${DATA_SOURCES_PATH}/gencode/hg38/gencode.v28.annotation.REORDERED.gtf

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME REFERENCE_NAME CONTIG START END"
  echo -e "Get the Gencode Gene info for the given location."
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e ""
}


#Display a message to std error:
function error() 
{
  echo "$1" 2>&1 
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

REF=$( echo $1 | tr '[:upper:]' '[:lower:]' )
CONTIG=$2
START=$3
END=$4

################################################################################

GENCODE="$GENCODE_HG19"
if [[ "${REF}" == "hg38" ]] ;then 
	GENCODE="$GENCODE_HG38"
fi

geneList=${REF}_geneInfo.tsv
[[ ! -e ${geneList} ]] && getGeneList > ${geneList}

awk "{ if (\$1 == \"${CONTIG}\" && \$4 <= ${START} && \$5 >= ${END}) {print} }" ${geneList}

