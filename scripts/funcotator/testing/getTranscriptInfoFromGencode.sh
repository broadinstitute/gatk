#!/usr/bin/env bash

################################################################################

#Setup variables for the script:
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=1

################################################################################


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

