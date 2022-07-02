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
D="."
#refAFile="${D}/GRCh37.p13.genome.dict"
#refAFile="${D}/human_g1k_v37.dict"
#refAFile="${D}/ucsc.hg19.dict"
#refBFile="${D}/Homo_sapiens_assembly19.dict"
refAFile="${D}/genome.hg38rg.fa.dict"
refBFile="${D}/Homo_sapiens_assembly38.dict"

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
tmpRefB=$( makeTemp ) 

# Create a file that has two columns: MD5Sum Contig
# Assumes MD5sum is in column 6 and Contig is in column 2:
cat ${refBFile} | awk '{print $6,$2}'  | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefB} 
#cat ${refBFile} | awk '{print $4,$2}'  | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefB} 

# Create a file that has two columns: MD5Sum Contig
# Assumes MD5sum is in column 4 and Contig is in column 2:
#cat ${refAFile} | awk '{print $4,$2}' | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefA}
# Assumes MD5sum is in column 5 and Contig is in column 2:
cat ${refAFile} | awk '{print $5,$2}' | grep '^M5' | sort -nk2 | sed -e 's#M5:##g' -e 's#SN:##g' > ${tmpRefA}

tmpJoined=$( makeTemp )

echo "Joining Files" 1>&2
join -1 1 -2 1 ${tmpRefB} ${tmpRefA} | sort -nk2 > ${tmpJoined} 

echo "Consolidating Ref B" 1>&2
while read line ; do
	md5=$( echo ${line} | awk {'print $1'} );
	sn=$(  echo ${line} | awk {'print $2'} );
	echo "  SN:${sn}" 1>&2
	grep "${md5}" $tmpJoined &> /dev/null
	r=$?
	if [[ $r -ne 0 ]] ; then
		echo -e "${md5} ${sn} ----" >> ${tmpJoined}
	fi
done < ${tmpRefB}

echo "Consolidating Ref A" 1>&2
while read line ; do
	md5=$( echo ${line} | awk {'print $1'} );
	sn=$(  echo ${line} | awk {'print $2'} );
	echo "  SN:${sn}" 1>&2
	grep "${md5}" $tmpJoined &> /dev/null
	r=$?
	if [[ $r -ne 0 ]] ; then
		echo -e "${md5} ---- ${sn}" >> ${tmpJoined}
	fi
done < ${tmpRefA}

echo "Creating Table..." 1>&2
echo -e "MD5\tB37SN (${refBFile})\tHG19SN (${refAFile})"
cat ${tmpJoined} | tr ' ' '\t' 

exit 0

