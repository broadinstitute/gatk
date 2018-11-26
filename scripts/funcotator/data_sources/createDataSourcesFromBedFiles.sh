#!/usr/bin/env bash

################################################################################
#
# DESCRIPTION:
#
# This script will take a list of BED files and create data sources for
# Funcotator from them.
# This process involves but is not limited to:
#    
#  - Remove tracks
#  - Making coordinates 1-based
#  - Adding column headers to top of file
#  - Creating a config file and directory structure for data sources.
#
# EXAMPLE:
#     ./createDataSourcesFromBedFiles.sh GATK_LAUNCH_SCRIPT DATA_SOURCES_DIRECTORY BED_FILE [BED_FILE BED_FILE ...]
#
# AUTHOR: Jonn Smith
#
################################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( python -c "import os;print os.path.realpath(\"${BASH_SOURCE[0]}\")" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=2
MAXARGS=99
PREREQUISITES="perl grep awk sed"

# Determine if this shell is interactive:
ISCALLEDBYUSER=true
[[ "${BASH_SOURCE[0]}" != "${0}" ]] && ISCALLEDBYUSER=false

# Required for the aliased call to checkPipeStatus:
shopt -s expand_aliases

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME DATA_SOURCES_DIRECTORY BED_FILE [BED_FILE BED_FILE ...]"
  echo -e "Create a data source folder for each given BED_FILE in the given DATA_SOURCES_DIRECTORY" 
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo
  echo -e 'This script will take a list of BED files and create data sources for'
  echo -e 'Funcotator from them.'
  echo -e 'This process involves but is not limited to:'
  echo -e '   '
  echo -e ' - Making coordinates 1-based'
  echo -e ' - Adding column headers to top of file'
  echo -e ' - Creating a config file and directory structure for data sources.'
  echo -e ''
  echo -e 'EXAMPLE:'
  echo -e '    ./createDataSourcesFromBedFiles.sh ~/gatk-4.0.0/gatk funcotator_dataSources.v1.5.20181119g autogen.bed autogen2.bed' 
  echo -e ""
  if [[ ${#PREREQUISITES} -ne 0 ]] ; then
    echo -e "Requires the following programs to be installed:"
    for prereq in ${PREREQUISITES} ; do 
      echo -e "  ${prereq}"
    done
    echo
  fi
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  MISSING PREREQUESITE(S)"
  echo -e "  4  GATK LAUNCH SCRIPT DOESN'T EXIST" 
  echo -e "  5  DATA SOURCES DIRECTORY DOES NOT EXIST"
  echo -e "  6  BED FILE DOES NOT EXIST"
  echo -e "  7  BED FILE DOES NOT HAVE .bed EXTENSION" 
  echo -e "  8  ERROR INDEXING FEATURE FILE" 
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
trap at_exit EXIT 

# Checks the bash built-in PIPESTATUS variable for any failures
# If given strings will output the string corresponding to the failure
# position in PIPESTATUS of any failures in the chain of commands 
# that occurred.
# This function should be called as `checkPipeStatus` as per the 
# alias below it.
function _checkPipeStatus() {

  local hadFailure=false

  for (( i = 0 ; i < ${#RC[@]} ; i++ )) ; do 
    st=${RC[i]}
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

function checkPrerequisites {
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

function printConfigFile() {

  local dsName=$1
  local srcFileName=$2
  local origin=$3

  echo -e "name = $dsName"
  echo -e "version = $(date +%Y-%m-%dT%H:%M:%S)"
  echo -e "src_file = $srcFileName"
  echo -e "origin_location = $origin"
  echo -e "preprocessing_script = $SCRIPTNAME" 
  echo -e ""
  echo -e "# Supported types:"
  echo -e "# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID"
  echo -e "# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location"
  echo -e "# gencode      -- Custom datasource class for GENCODE"
  echo -e "# cosmic       -- Custom datasource class for COSMIC"
  echo -e "# vcf          -- Custom datasource class for Variant Call Format (VCF) files"
  echo -e "type = locatableXSV"
  echo -e ""
  echo -e "# Required field for GENCODE files."
  echo -e "# Path to the FASTA file from which to load the sequences for GENCODE transcripts:"
  echo -e "gencode_fasta_path ="
  echo -e ""
  echo -e "# Required field for simpleXSV files."
  echo -e "# Valid values:"
  echo -e "#     GENE_NAME"
  echo -e "#     TRANSCRIPT_ID"
  echo -e "xsv_key = GENE_NAME"
  echo -e ""
  echo -e "# Required field for simpleXSV files."
  echo -e "# The 0-based index of the column containing the key on which to match"
  echo -e "xsv_key_column = 1"
  echo -e ""
  echo -e "# Required field for simpleXSV AND locatableXSV files."
  echo -e "# The delimiter by which to split the XSV file into columns."
  echo -e "xsv_delimiter = \\\\t"
  echo -e ""
  echo -e "# Required field for simpleXSV files."
  echo -e "# Whether to permissively match the number of columns in the header and data rows"
  echo -e "# Valid values:"
  echo -e "#     true"
  echo -e "#     false"
  echo -e "xsv_permissive_cols = true"
  echo -e ""
  echo -e "# Required field for locatableXSV files."
  echo -e "# The 0-based index of the column containing the contig for each row "
  echo -e "contig_column = 0"
  echo -e ""
  echo -e "# Required field for locatableXSV files."
  echo -e "# The 0-based index of the column containing the start position for each row "
  echo -e "start_column = 1"
  echo -e ""
  echo -e "# Required field for locatableXSV files."
  echo -e "# The 0-based index of the column containing the end position for each row "
  echo -e "end_column = 2"
  echo -e ""
}

################################################################################

# Do all the interactive stuff if the script is called by a user.
# The bash equivalent to `if __name__ == '__main__':
if ${ISCALLEDBYUSER} ; then

  ################################################################################
  # Do some basic pre-processing before we do the work:

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
  [[ $r -ne 0 ]] && exit 3

  # Check for help:
  for arg in "${@}" ; do
    case "$arg" in
      -h|--h|--help|-help|help)
        usage;
        exit 0;
        ;;
    esac
  done

  #Read args:
  

  ################################################################################
  # Do the work here:

	GATK_LAUNCH_SCRIPT=$1
	shift
  # Make sure folder exists:
  [[ ! -f ${GATK_LAUNCH_SCRIPT} ]] && error "GATK launch script not exist: ${GATK_LAUNCH_SCRIPT}" && exit 4 

  DATA_SOURCES_DIR=$1
  shift 
  # Make sure folder exists:
  [[ ! -d ${DATA_SOURCES_DIR} ]] && error "Data sources directory does not exist: ${DATA_SOURCES_DIR}" && exit 5 

  # Get the BED files and process them:
  for BED_FILE in "${@}" ; do 

    error "Creating data source for ${BED_FILE}"

    # Make sure BED_FILE exists:
    [[ ! -f ${BED_FILE} ]] && error "Bed file does not exist: ${BED_FILE}" && exit 6 

    # Make sure BED_FILE has .bed as an extension:
    echo ${BED_FILE} | sed 's#\(.*\)\.[bBeEdD]*#\1#g' &> /dev/null
    checkPipeStatus "problem with echo" "Bed file does not have .bed extension: ${BED_FILE}" 
    [[ ! $? -eq 0 ]] && exit 6

    # Create a temporary directory for this BED file data source:
    tmpDataSourceDir=$( mktemp -d )
    bedDataSourceName=$( basename ${BED_FILE} | sed 's#\(.*\)\.[bBeEdD]*#\1_bed#g' | tr -d ' ' )
    bedFileLocalName="${bedDataSourceName}.tsv"

    # Make folders for HG19 and HG38:     
    mkdir -p ${tmpDataSourceDir}/hg19 ${tmpDataSourceDir}/hg38 

    # Create config files for both HG19 and HG38:
    # Note: they'll be identical so we just make one and copy it:
    printConfigFile ${bedDataSourceName} ${bedFileLocalName} ${BED_FILE} > ${tmpDataSourceDir}/hg19/${bedDataSourceName}.config
    cp ${tmpDataSourceDir}/hg19/${bedDataSourceName}.config ${tmpDataSourceDir}/hg38/${bedDataSourceName}.config 

    # Create header for file:
    echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts' > ${tmpDataSourceDir}/hg19/${bedFileLocalName}

    # Transform the input BED file:
    # Get rid of tracks, add 1 to positions, 
    grep -v '^track' ${BED_FILE} | perl -pe 's#[\t ]+#\t#g' | awk '{ $2=$2+1; $3=$3+1; printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12 }' >>  ${tmpDataSourceDir}/hg19/${bedFileLocalName}

		# Index the tsv file:
		error "Indexing file: ${tmpDataSourceDir}/hg19/${bedFileLocalName}"
		echo "################################################################################"
		${GATK_LAUNCH_SCRIPT} IndexFeatureFile -F ${tmpDataSourceDir}/hg19/${bedFileLocalName} 
		r=$?
		[[ $r -ne 0 ]] && error "Error indexing file: ${tmpDataSourceDir}/hg19/${bedFileLocalName}" && exit 8

    # Copy the files into the hg38 directory:
    cp ${tmpDataSourceDir}/hg19/${bedFileLocalName}* ${tmpDataSourceDir}/hg38/.
		
		# Put the new data source in the data sources folder:
    mv ${tmpDataSourceDir} ${DATA_SOURCES_DIR}/${bedDataSourceName}
		echo "################################################################################"
  done

  error "All data sources created in folder: ${DATA_SOURCES_DIR}"
  error "DONE"

fi

