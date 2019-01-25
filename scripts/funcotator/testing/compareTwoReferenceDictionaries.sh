#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script will create a table comparing two internally configured reference
# dictionary files.
# It must be internally configured to point at valid reference fasta dictionaries
# and may need to be modified to correctly grab the md5 and contig columns for
# the given files.
# 
# This will not work for you out-of-the-box unless you are Jonn Smith.
#
# EXAMPLE:
#     ./compareTwoReferenceDictionaries.sh
#
# AUTHOR: Jonn Smith
#
###############################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( readlink "${BASH_SOURCE[0]}" || echo "${BASH_SOURCE[0]}" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=0

################################################################################

# Change these to point to reference dictionaries:
D="/Users/jonn/Development/references"
refGRCh37File="${D}/GRCh37.p13.genome.dict"
refHumanG1Kv37File="${D}/human_g1k_v37.dict"
refHg19File="${D}/ucsc.hg19.dict"
refB37File="${D}/Homo_sapiens_assembly19.dict"

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME" 
	echo -e "Compares two internally-selected sequence dictionaries and displays"
	echo -e "the differences in a table."
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
  echo -e "  3  UNKNOWN ARGUMENT"
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

function preserveTempFile()
{
	cp ${!1} ${1}.txt
}

function appendMissingDataToOutputFile() 
{
	local tmpRefFile=$1
	local outFile=$2

	# Get all contigs from hg19 that are not in the joined file and annotate them:
	while read line ; do
		md5=$( echo ${line} | awk {'print $1'} );
		contig=$(  echo ${line} | awk {'print $2'} );
		grep -q "${md5}" ${outFile}
		r=$?
		if [[ $r -ne 0 ]] ; then
			# Now we need to see what other files it's in:
			b37Contig=$( grep "${md5}" ${tmpRefB37} | awk '{print $2}' )
			hg19Contig=$( grep "${md5}" ${tmpRefHg19} | awk '{print $2}' )
			grch37Contig=$( grep "${md5}" ${tmpRefGRCh37} | awk '{print $2}' )
			humanG1kv37Contig=$( grep "${md5}" ${tmpRefHumanG1Kv37} | awk '{print $2}' )

			[[ ${#b37Contig} -eq 0 ]] && b37Contig="----"
			[[ ${#hg19Contig} -eq 0 ]] && hg19Contig="----"
			[[ ${#grch37Contig} -eq 0 ]] && grch37Contig="----"
			[[ ${#humanG1kv37Contig} -eq 0 ]] && humanG1kv37Contig="----"

			echo -e "${md5} ${b37Contig} ${hg19Contig} ${grch37Contig} ${humanG1kv37Contig}" >> ${outFile}
		fi
	done < ${tmpRefFile}
}

################################################################################

trap at_exit EXIT 

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

#Read args:
while [ $# -gt 0 ] ; do

  case "$1" in
    -h|--h|--help|-help|help)
      usage;
      exit 0;
      ;;
    *)
      error "${SCRIPTNAME}: invalid option: $1" 
      simpleUsage
      echo "Try \`$SCRIPTNAME --help' for more information."
      exit 3;
    ;; 
  esac
  
  #Get next argument in $1:
  shift
done

################################################################################

tmpRefA=$( makeTemp )
tmpRefGRCh37=$( makeTemp )
tmpRefHg19=$( makeTemp )
tmpRefB37=$( makeTemp ) 
tmpRefHumanG1Kv37=$( makeTemp )

# Create a file that has two columns for each reference: MD5Sum Contig

echo "Getting contig information from sequence dictionaries..."

# Assumes MD5sum is in column 6 and Contig is in column 2:
cat ${refB37File} | awk '{print $6,$2}'  | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefB37} 

# hg19 has MD5sum is in column 5 and Contig is in column 2:
cat ${refHg19File} | awk '{print $5,$2}' | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefHg19}

# GRCh37 has MD5sum is in column 4 and Contig is in column 2:
cat ${refGRCh37File} | awk '{print $4,$2}' | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefGRCh37}

# human_g1k_v37 has MD5sum is in column 5 and Contig is in column 2:
cat ${refHumanG1Kv37File} | awk '{print $5,$2}' | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefHumanG1Kv37}

# Create temp out file:
tmpJoined=$( makeTemp )

# Join the common data together in a table
echo "Joining sequence dictionary data..."
join -1 1 -2 1 ${tmpRefB37} ${tmpRefHg19} | join -1 1 -2 1 - ${tmpRefGRCh37} | join -1 1 -2 1 - ${tmpRefHumanG1Kv37} | sort -nk2 > ${tmpJoined}

echo "Getting B37 Missing Contigs..."
appendMissingDataToOutputFile ${tmpRefB37} ${tmpJoined}

echo "Getting Hg19 Missing Contigs..."
appendMissingDataToOutputFile ${tmpRefHg19} ${tmpJoined}

echo "Getting GRCh37 Missing Contigs..."
appendMissingDataToOutputFile ${tmpRefGRCh37} ${tmpJoined}

echo "Getting HumanG1Kv37 Missing Contigs..."
appendMissingDataToOutputFile ${tmpRefHumanG1Kv37} ${tmpJoined}

# Display some extra information about the files:
echo -e "B37:  ${refB37File}"
echo -e "HG19: ${refHg19File}"
echo -e "GRCh37: ${refGRCh37File}"
echo -e "human_g1k_v37: ${refHumanG1Kv37File}"
echo ""

# Add headers:
echo -e "MD5 B37 HG19 GRCh37 HumanG1Kv37\n$(cat ${tmpJoined})" > ${tmpJoined}

# Fix columns to be in a better order for comparison:
finalOut=$( makeTemp )
awk '{print $1,$5,$2,$3,$4}' ${tmpJoined} > ${finalOut}

# Make it pretty:
column -t -s' ' ${finalOut} 

exit 0

