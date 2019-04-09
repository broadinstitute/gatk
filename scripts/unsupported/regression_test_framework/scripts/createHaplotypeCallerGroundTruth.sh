#!/usr/bin/env bash

################################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( python -c "import os;print os.path.realpath(\"${BASH_SOURCE[0]}\")" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=0
MAXARGS=0
PREREQUISITES="gsutil"

# Determine if this shell is interactive:
ISCALLEDBYUSER=true
[[ "${BASH_SOURCE[0]}" != "${0}" ]] && ISCALLEDBYUSER=false

# Required for the aliased call to checkPipeStatus:
shopt -s expand_aliases

################################################################################

INPUT_BUCKET_URL="gs://haplotypecallerspark-evaluation/inputData"
OUTPUT_BUCKET_URL="gs://haplotypecallerspark-evaluation/groundTruth"

HAPLOTYPE_CALLER_WDL="${SCRIPTDIR}/../wdl/haplotypeCaller.wdl"

CROMSHELL_SCRIPT="/Users/jonn/Development/cromshell/cromshell"

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME"
  echo -e "Run the baseline data through HaplotypeCaller to get ground truth for comparison."
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo -e ""
  echo -e "Runs HaplotypeCaller on the bam files located in:"
  echo -e "  ${INPUT_BUCKET_URL}"
  echo -e ""
  # List prerequisites, if any:
  if [[ ${#PREREQUISITES} -ne 0 ]] ; then
    echo -e "Requires the following programs to be installed:"
    for prereq in ${PREREQUISITES} ; do 
      echo -e "  ${prereq}"
    done
    echo
  fi
  echo -e "  -t dockertag the full name for where to upload the docker image"
  echo -e "  -b BRANCH the branch of gatk off of which to build an image"
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  UNKNOWN ARGUMENT"
  echo -e "  4  MISSING PREREQUISITE"
  echo -e ""
}

#Display a message to std error:
function error() 
{
  echo "$1" 1>&2 
}

# Make a temporary file that will be cleaned after this script finishes. 
function makeTemp()
{
  local f
  f=$( mktemp )
  TMPFILELIST="${TMPFILELIST} $f"
  echo $f
}
TMPFILELIST=''

# Clean all temporary files made in this script:
function cleanTempVars()
{
  rm -f ${TMPFILELIST}
}

# Function run in the exit trap.
function at_exit()
{
  cleanTempVars
}

# Checks the bash built-in PIPESTATUS variable for any failures
# If given strings will output the string corresponding to the failure
# position in PIPESTATUS of any failures in the chain of commands 
# that occurred.
# This function should be called as `checkPipeStatus` as per the 
# alias below it.
function _checkPipeStatus()
{

  local hadFailure=false

  for (( i = 0 ; i < ${#RC[@]} ; i++ )) ; do 
    st=${RC[i]}
    echo "st = ${st}"
    if [ $st -ne 0 ] ; then
      # If we were passed a description string for this error in the pipe, then we print it:
      let argIndex=$i+1
      description=${!argIndex}
      [[ ${#description} -ne 0 ]] && error "$description"
      hadFailure=true  
    fi
  done

  if $hadFailure ; then
    return 10
  fi
  return 0
}
alias checkPipeStatus='RC=( "${PIPESTATUS[@]}" );_checkPipeStatus'

function checkPrerequisites
{
  local missing=""
  local foundMissing=false
  for c in ${@} ; do
    which $c &> /dev/null
    r=$?
    [[ $r -ne 0 ]] && foundMissing=true && missing="${missing}$c "
  done

  if $foundMissing ; then
    error "Error: the following commands could not be found:"
    error "  $missing"
    error "Please install them and try again."
    if [[ $(uname) == 'Darwin' ]] ; then
      error "NOTE: You likely must use homebrew to install them."
    fi
    return 1
  fi
  return 0
}

################################################################################

function createJsonInputsForHaplotypeCaller() 
{
    local doGvcf=$1
    local inputBam=$2
    local inputBai=$3
    local outputVcfName=$4

    echo "{"
    echo "  \"HaplotypeCaller.gatk_docker\": \"broadinstitute/gatk-nightly:2018-12-07-4.0.11.0-88-g2ae01efda-SNAPSHOT\","
    echo ""
    echo "  \"HaplotypeCaller.input_bam\": \"${inputBam}\","
    echo "  \"HaplotypeCaller.input_bam_index\": \"${inputBai}\","
    echo "  \"HaplotypeCaller.out_vcf_name\": \"${outputVcfName}\","
    echo ""
    echo "  \"HaplotypeCaller.interval_list\": \"gs://haplotypecallerspark-evaluation/inputData/whole_exome_illumina_coding_v1.Homo_sapiens_assembly19.targets.interval_list\","
    echo ""
    echo "  \"HaplotypeCaller.ref_fasta\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta\","
    echo "  \"HaplotypeCaller.ref_fasta_index\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai\","
    echo "  \"HaplotypeCaller.ref_fasta_dict\": \"gs://broad-references/hg19/v0/Homo_sapiens_assembly19.dict\","
    echo ""
    echo ""
    echo "  \"HaplotypeCaller.gvcf_mode\": \"${doGvcf}\","
    echo ""
    echo "  \"HaplotypeCaller.disk_space_gb\": \"512\","
    echo "  \"HaplotypeCaller.boot_disk_size_gb\": \"64\","
    echo "  \"HaplotypeCaller.mem_gb\": \"32\""
    echo "}"
}

################################################################################

trap at_exit EXIT 

################################################################################

# Do all the interactive stuff if the script is called by a user.
# The bash equivalent to `if __name__ == '__main__':
if ${ISCALLEDBYUSER} ; then

	#Check given arguments:
	if [[ $# -gt $MAXARGS ]] ; then
		usage
		exit 1
	elif [[ $# -lt $MINARGS ]] ; then
		usage
		exit 2
	fi

  # Make sure we have all the required commands:
  checkPrerequisites $PREREQUISITES 
  r=$?
  [[ $r -ne 0 ]] && exit 4

	#----------------------------------
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

	#----------------------------------
	# Pre-run sanity checks:
    [ ! -f ${CROMSHELL_SCRIPT} ] && error "Cromshell script does not exist: ${CROMSHELL_SCRIPT}" && exit 5

	#----------------------------------
	# Do real work here.
    for inputBamName in $(gsutil ls ${INPUT_BUCKET_URL}/*.bam | sort | uniq) ; do

#        echo "${inputBamName}" | grep -qE 'G96830|G96831|G96832|G94982|359781|412726|445394|472246|506817|538834|572804|603388|633960|656480|679060'
#        echo "${inputBamName}" | grep -qEv 'G96830|G96831|G96832|G94982'
#        r=$?
#        [ $r -eq 0 ] && continue

        for doGvcf in true false ; do
            b=$( basename ${inputBamName} )
            indexFile=$( echo ${b} | sed 's#.bam$#.bai#g' )

            gvcfModeSuffix=""
            outVariantFile=$( echo ${b} | sed 's#.bam$#.HC.vcf#g' )

            if $doGvcf ; then
                outVariantFile=$( echo ${b} | sed 's#.bam$#.HC.g.vcf#g' )
                gvcfModeSuffix=".gvcf"
            fi

#            echo "${inputBamName}" | grep -qEv 'NexPond-412726.ba|506817|633960'
#            r=$?
#            [ $r -eq 0 ] && continue
#            ! $doGvcf && continue

            tmpJsonFile="haplotypeCaller.${b}${gvcfModeSuffix}.json"
            createJsonInputsForHaplotypeCaller $doGvcf ${inputBamName} ${INPUT_BUCKET_URL}/${indexFile} ${outVariantFile} > ${tmpJsonFile}

            echo "Submitting ${tmpJsonFile}"
            ${CROMSHELL_SCRIPT} submit ${HAPLOTYPE_CALLER_WDL} ${tmpJsonFile}
            sleep 1
        done
    done
fi

