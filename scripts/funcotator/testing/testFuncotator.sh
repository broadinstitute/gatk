#!/usr/bin/env bash

################################################################################
#
# WARNING: THIS SCRIPT IS UNSUPPORTED!
# USE AT YOUR OWN RISK
#
# DESCRIPTION:
#
# This script builds and runs Funcotator on different data files, depending on
# the options given.
# Must be internally configured to point to reference FASTA files and 
# Funcotator data sources.
#
# Must be run from GATK development directory.
# 
# Will not work for you out-of-the-box unless you are Jonn Smith.
#
# For more information run with --help
# 
# EXAMPLEs:
#     ./testFuncotator.sh
#     ./testFuncotator.sh -c -u
#     ./testFuncotator.sh -c -u -t
#     ./testFuncotator.sh -c -t -38
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

GATKDIR=${SCRIPTDIR}/../../../

BASE_DIR=/Users/jonn/Development

################################################################################

doForceRun=false

doUnitTests=false
doRunLargeTests=false

doClean=false

REF_VER=hg19
OUT_FORMAT=VCF
useV16DataSources=false
useAOUDataSources=false
useCloudDataSources=false

MANUAL_MODE=false

################################################################################

# Change this to point to your funcotator data sources folder:
DATA_SOURCES_PATH_16=${BASE_DIR}/funcotator_dataSources.v1.6.20190124s
DATA_SOURCES_PATH=${BASE_DIR}/funcotator_dataSources_latest
DATA_SOURCES_PATH_GERMLINE=${BASE_DIR}/funcotator_dataSources_germline_latest
HG19=${BASE_DIR}/references/Homo_sapiens_assembly19.fasta
HG38=${BASE_DIR}/references/Homo_sapiens_assembly38.fasta

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME [-c] [-cloud] [-u] [-t] [-19|-38] [-MAF|-VCF] [-AOU]"
  echo -e "Build and run Funcotator."
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo -e "Can clean, run tests, and run large file tests."
  echo -e ""
  echo -e "MUST be run from the GATK development directory." 
  echo -e ""
  echo -e "Will by default (with no options) build GATK/Funcotator."
  echo -e "For large file tests, defaults to hg19 tests with VCF output."
  echo -e ""
  echo -e "The following options are available:"
  echo -e "  -c                                          clean GATK/Funcotator"
  echo -e "  -u                                          run all tests in the Funcotator Package."
  echo -e "                                              (org.broadinstitute.hellbender.tools.funcotator)"
  echo -e "  -t                                          run Funcotator on a large data file"
  echo -e "                                              (internally configured)"
  echo -e "  -f                                          Force a run, ignoring build and file checks."
  echo -e "  -19                                         run with hg19 data sources/reference/input file"
  echo -e "                                              (default)"
  echo -e "  -38                                         run with hg38 data sources/reference/input file"
  echo -e "  -v16                                        use the Funcotator v1.6 data sources."
  echo -e "  -MAF                                        create MAF output"
  echo -e "  -VCF                                        create VCF output (default)"
  echo -e "  -cloud                                      use cloud data sources"
  echo -e "  -AOU                                        use the All of Us/Clinical Pipeline data sources"
  echo -e "  -M REF_VER REFERENCE INPUT DATA_SOURCES     run in MANUAL mode, providing all necessary input"
  echo -e "                                              REF_VER      - a string for the reference version"
  echo -e "                                              REFERENCE    - reference FASTA file"
  echo -e "                                              INPUT        - input VCF file"
  echo -e "                                              DATA_SOURCES - path to FUNCOTATOR data sources folder"
  echo -e ""
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  REQUIRED FILE DOES NOT EXIST"
  echo -e "  4  REQUIRED DIRECTORY DOES NOT EXIST"
  echo -e "  5  INCORRECT ARGUMENTS"
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

function assertFileExists() {
  [[ ! -f $1 ]] && error "Error: File does not exist: $1" && exit 3
}

function assertDirectoryExists() {
  [[ ! -d $1 ]] && error "Error: Directory does not exist: $1" && exit 4
}

################################################################################

trap at_exit EXIT 

################################################################################

function assertInputFilesExist() {
  assertFileExists ${INPUT}
  assertFileExists ${REF}
 
   [[ ! -d $DATA_SOURCES_PATH ]] && error "Warning: Data sources may not exist ${DATA_SOURCES_PATH}" && error "Ignore this if data sources directory is in the cloud."
}

################################################################################

#Read args:
while [ $# -gt 0 ] ; do

  case "$1" in
    -f)
    doForceRun=true
    ;;
    -c)
    doClean=true
    ;;
    -u)
    doUnitTests=true
    ;;
    -19)
    REF_VER=hg19
    ;;
    -38)
    REF_VER=hg38
    ;;
    -VCF)
    OUT_FORMAT=VCF
    ;;
    -MAF)
    OUT_FORMAT=MAF
    ;;
    -v16)
    useV16DataSources=true
    ;;
    -AOU)
    useAOUDataSources=true
    ;;
    -cloud)
    useCloudDataSources=true
    ;;
    -t)
    doRunLargeTests=true
    ;;
    -M)
      shift
      REF_VER=$1
      shift
      REF=$1
      shift
      INPUT=$1
      shift 
      DATA_SOURCES_PATH=$1
      MANUAL_MODE=true
      # Validate our args:
      if [[ ${#REF} -eq 0 ]] || [[ ${#INPUT} -eq 0 ]] || [[ ${#DATA_SOURCES_PATH} -eq 0 ]] ; then
        error "Error: For manual mode you must specify a reference version, reference fasta, input file, and data sources directory." && exit 5
      fi
    ;;
    --help)
      usage
      exit 0
    ;;
    *)
    ;; 
  esac
  
  #Get next argument in $1:
  shift
done

################################################################################

r=1
if ${doForceRun} ; then
  r=0
else
  if ${doClean} ; then
    ${GATKDIR}/gradlew clean compileJava compileTestJava installDist
    r=$?
  else
    ${GATKDIR}/gradlew compileJava compileTestJava installDist
    r=$?
  fi 
  
  if [[ $r -eq 0 ]] && ${doUnitTests} ; then
    echo "################################################################################"
    echo "## Running Unit Tests... "
    ${GATKDIR}/gradlew test \
      --tests org.broadinstitute.hellbender.tools.funcotator* \
      --tests org.broadinstitute.hellbender.utils.codecs.xsvLocatableTable* \
      --tests org.broadinstitute.hellbender.utils.codecs.gtf* \
      --tests org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.SimpleAnnotatedIntervalWriterUnitTest* \
      --tests org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval.AnnotatedIntervalCollectionUnitTest* \
			--stacktrace > >(tee -a FUNCOTATOR.unitTest.stdout.log) 2> >(tee -a FUNCOTATOR.unitTest.stderr.log >&2)
    r=$?
  fi
  
  ################################################################################
  
  if [[ $r -eq 0 ]] && $MANUAL_MODE ; then 
  
    echo "################################################################################"
    echo "## Running MANUAL Test... "
    echo
    echo "########################################"
    echo "## Using Reference: ${REF_VER}              ##"
    echo "########################################"
    
    OUT_FORMAT_LOWER=$( echo "${OUT_FORMAT}" | tr 'A-Z' 'a-z' )
    OUT_FILE_NAME=FUNCOTATOR_OUT.${OUT_FORMAT_LOWER}
  
    assertInputFilesExist
  
    ${GATKDIR}/gatk Funcotator \
      -V ${INPUT} \
      -O ${OUT_FILE_NAME} \
      -R ${REF} \
      --verbosity DEBUG \
      --data-sources-path ${DATA_SOURCES_PATH} \
      --ref-version ${REF_VER} \
      --output-file-format ${OUT_FORMAT} -- --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
    
    r=$?
    exit $r
  fi
fi

if [[ $r -eq 0 ]] && ${doRunLargeTests} ; then
  
  echo "################################################################################"
  echo "## Running Large Tests... "
  echo
  echo
  echo "########################################"
  echo "## Using Reference: ${REF_VER}              ##"
  echo "########################################"

  if [[ "${REF_VER}" == "hg19" ]] ; then
    INPUT=${BASE_DIR}/NON_PUBLIC/0816201804HC0_R01C01.vcf
    #INPUT=${BASE_DIR}/funcotator_debugging/problem_variant_fixed.vcf
    #INPUT=${BASE_DIR}/gatk/src/test/resources/large/funcotator/regressionTestVariantSet1.vcf
    #INPUT=${BASE_DIR}/gatk/src/test/resources/large/funcotator/regressionTestVariantSet2.vcf
    #INPUT=${BASE_DIR}/gatk/src/test/resources/large/funcotator/regressionTestHg19Large.vcf
    #INPUT=${BASE_DIR}/gatk/hg38_trio_liftoverb37.vcf
    #INPUT=${BASE_DIR}/gatk/tmp.vcf
    #INPUT=${BASE_DIR}/gatk/tmp2.vcf
    #INPUT=${BASE_DIR}/data_to_run/problem_samples/splice_site_should_not_be_splice_site/error_case.vcf
    
    #HG19=${BASE_DIR}/references/ucsc.hg19.fasta
    #HG19=${BASE_DIR}/references/ucsc.hg19.fasta
    #HG19=${BASE_DIR}/references/GRCh37.p13.genome.fasta
    REF=$HG19
  else
    INPUT=${BASE_DIR}/FUNCOTATOR_LARGE_TEST_INPUTS/hg38_trio.vcf
		#INPUT=${BASE_DIR}/gatk/src/test/resources/large/funcotator/regressionTestVariantSetHG38.vcf
    #INPUT=${BASE_DIR}/tmp/cohort24_23_seg.subset.vcf
    #INPUT=${BASE_DIR}/gatk/tmp.38.vcf
    #INPUT=${BASE_DIR}/gatk/src/test/resources/org/broadinstitute/hellbender/tools/funcotator/badDataOneAlleleDepthValue_hg38.vcf
    REF=$HG38
  fi

  # Use the AOU data sources if we need them:
	$useV16DataSources && echo "Using v1.6 data sources." && DATA_SOURCES_PATH=${DATA_SOURCES_PATH_16}
  
	# Use the AOU data sources if we need them:
  $useAOUDataSources && echo "Using AOU data sources." && DATA_SOURCES_PATH=${BASE_DIR}/funcotator_dataSources_germline_latest

  # Use cloud data sources if we need them:
  #$useCloudDataSources && echo "Using cloud data sources." && DATA_SOURCES_PATH=${BASE_DIR}/gatk/src/test/resources/large/funcotator/funcotator_dataSources_cloud_gnomad/
  $useCloudDataSources && echo "Using cloud data sources." && DATA_SOURCES_PATH=gs://broad-dsde-methods-jonn/funcotator_dataSources.v1.5.20181119s

  OUT_FORMAT_LOWER=$( echo "${OUT_FORMAT}" | tr 'A-Z' 'a-z' )
  OUT_FILE_NAME=FUNCOTATOR_OUT.${OUT_FORMAT_LOWER}
  
  if ! ${doForceRun} ; then 
    assertInputFilesExist
  fi

  time ${GATKDIR}/gatk Funcotator \
    -V ${INPUT} \
    -O ${OUT_FILE_NAME} \
    -R ${REF} \
    --verbosity DEBUG \
    --data-sources-path ${DATA_SOURCES_PATH} \
    --ref-version ${REF_VER} \
    --output-file-format ${OUT_FORMAT} \
    --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx16g' \
    --cloud-index-prefetch-buffer 40 \
    --cloud-prefetch-buffer 80 > >(tee -a FUNCOTATOR.stdout.log) 2> >(tee -a FUNCOTATOR.stderr.log >&2)

  r=$?
fi

exit $r

