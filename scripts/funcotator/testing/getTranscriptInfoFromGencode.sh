#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script displays Gencode Transript information given a reference version
# and a list of transcript IDs.
# It must be internally configured to point at a valid funcotatior data sources
# directory, and will not work for you out-of-the-box unless you are Jonn Smith.
#
# EXAMPLE:
#     ./getTranscriptInfoFromGencode.sh hg19 ENST00000602725.1 ENST00000424989.1 ENST00000373087.6
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

# Change this path to your data sources directory:
DATA_SOURCES_PATH=/Users/jonn/Development/funcotator_dataSources_latest

GENCODE_HG19=${DATA_SOURCES_PATH}/gencode/hg19/gencode.v19.annotation.REORDERED.gtf
GENCODE_HG38=${DATA_SOURCES_PATH}/gencode/hg38/gencode.v28.annotation.REORDERED.gtf

REF_VER=HG19

################################################################################

txIds=''

################################################################################

#Read args:
while [ $# -gt 0 ] ; do

	case "$1" in
		-38)
			REF_VER=hg38
			;;
		*)
			txIds="$txIds $1"
			;; 
	esac

	#Get next argument in $1:
	shift
done

################################################################################

GENCODE="$GENCODE_HG19"
if [[ "${REF_VER}" == "HG38" ]] ;then 
	GENCODE="$GENCODE_HG38"
fi

for txId in $txIds ; do 
	echo
	TXINFO=$( grep -m1 "$txId" $GENCODE )
	echo $TXINFO | grep "$txId"

	txStart=$(echo $TXINFO | awk '{print $4}')
	txEnd=$(echo $TXINFO | awk '{print $5}')
	txLen=$( echo "$txEnd - $txStart" | bc )

	echo "TX Length: $txLen"

	echo 
	echo "-------------------"
done

