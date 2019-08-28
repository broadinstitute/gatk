#!/bin/bash

# terminate script on error or if a command fails before piping to another command
set -eu
set -o pipefail

########################################################################################################################
# helper functions
show_help() {
cat << EOF
Manage the SV discovery pipeline on Google Cloud Services (GCS) cluster
Syntax
  manage_sv_pipeline.sh [Options] GATK_Folder Project_Name \\
                        GCS_BAM_File GCS_REFERENCE_File \\
                        [args to SV discovery pipeline]
  Options:
    -h or --help:
      display help and exit
    -q or --quiet:
      run without prompting for user input
    -n or --no-fastq
      don't copy fastq files from cluster to cloud storage
    -u [GCS_USER] or --user [GCS_USER]:
      set GCS username (defaults to ${USER})
    -i [INIT_SCRIPT] or --init [INIT_SCRIPT]
      use specified spark init script
      (defaults to \${GATK_FOLDER}/scripts/sv/default_init.sh)
    -s [GCS_SAVE_PATH] or --save [GCS_SAVE_PATH]
      save results in specified bucket/folder
      (defaults to \$PROJECT_NAME/\$GCS_USER)
    -t [MAX_LIFE] or --life [MAX_LIFE]:
      life time of cluster (e.g. "2h" for two hours or "1d" for 1 day)
      (defaults to \$CLUSTER_MAX_LIFE_HOURS if defined, otherwise 4h)
    -d [MAX_IDLE] or --idle [MAX_IDLE]:
      maximum idling time for cluster (e.g. "2h" for two hours or "10m"
      for 10 minutes).
      (defaults to \$CLUSTER_MAX_IDLE_MINUTES if defined, otherwise 60m)
    -L [TOOL_NAME] or --launch-tool [TOOL_NAME]
      specify GATK_SV_TOOL to launch
      (defaults to "StructuralVariationDiscoveryPipelineSpark")

  Mandatory Positional Arguments:
    GATK_Folder: path to local copy of GATK
    Project_Name: name of GCS project
    GCS_BAM_File: path to .bam file hosted in GCS
                  (e.g. gs://bucket/path/to/file.bam)
    GCS_REFERENCE_File: path to reference .2bit file hosted in GCS

  Optional Positional Arguments:
    [args to SV discovery pipeline]: additional arguments will be passed
                  to StructuralVariationDiscoveryPipelineSpark.
                  NOTE: these args will be eval-ed within
                  runWholePipeline.sh, so it is possible to pass
                  variable names as strings, provided they reference
                  variables defined within runWholePipeline.sh

This script manages the SV discovery pipeline by calling other scripts
a) ensure that the gcloud utility is up to date
b) ensure that local GATK .git hash matches local GATK .jar file.
c) create a Google Dataproc cluster used for running the GATK-SV
   pipeline (call create_cluster.sh)
d) run the analysis (call runWholePipeline.sh)
e) copy results to a datetime and git-version stamped directory
   on GCS, along with logged console output (call copy_sv_results.sh)
f) shut down the cluster (or self-terminate based on the configuration)
EOF
}

throw_error() {
    echo "$1" >&2
    exit 1
}

########################################################################################################################
# parses arguments provided to this script
QUIET=${SV_QUIET:-"N"}
GCS_USER=${GCS_USER:-${USER}}
CLUSTER_MAX_LIFE_HOURS=${CLUSTER_MAX_LIFE_HOURS:-4h}
CLUSTER_MAX_IDLE_MINUTES=${CLUSTER_MAX_IDLE_MINUTES:-60m}
GATK_SV_TOOL=${GATK_SV_TOOL:-"StructuralVariationDiscoveryPipelineSpark"}
COPY_FASTQ=${COPY_FASTQ:-"Y"}
SV_UPDATE_GCLOUD=${SV_UPDATE_GCLOUD:-true}

function update_gcloud() {
    if [ $SV_UPDATE_GCLOUD = true ]; then
        if [[ "${QUIET}" == "Y" ]]; then
            gcloud components update --quiet
        else
            gcloud components update
        fi
    fi
}
# try to update gcloud, on error ignore and continue. This can happen if
# a package manager is maintaining gcloud
update_gcloud || true

while [ $# -ge 1 ]; do
    case $1 in
        -h|-\?|--help)
            show_help
            exit 0
            ;;
        -q|--quiet)
            QUIET="Y"
            shift
            ;;
        -n|--no-fastq)
            COPY_FASTQ="N"
            shift
            ;;
        -u|--user)
            if [ $# -ge 2 ]; then
                GCS_USER="$2"
                shift 2
            else
                throw_error "--user requires a non-empty argument"
            fi
            ;;
        --user=?*)
            GCS_USER=${1#*=} # remove everything up to = and assign rest to user
            shift
            ;;
        -i|--init)
            if [ $# -ge 2 ]; then
                INIT_SCRIPT="$2"
                shift 2
            else
                throw_error "--init requires a non-empty argument"
            fi
            ;;
        --init=?*)
            INIT_SCRIPT=${1#*=} # remove everything up to = and assign rest to init
            shift
            ;;
        -s|--save)
            if [ $# -ge 2 ]; then
                GCS_SAVE_PATH="$2"
                shift 2
            else
                throw_error "--save requires a non-empty argument"
            fi
            ;;
        --save=?*)
            GCS_SAVE_PATH=${1#*=} # remove everything up to = and assign rest to save
            shift
            ;;
        -t|--life)
            if [ $# -ge 2 ]; then
                CLUSTER_MAX_LIFE_HOURS="$2"
                shift 2
            else
                throw_error "--life requires a non-empty argument"
            fi
            ;;
        --life=?*)
            CLUSTER_MAX_LIFE_HOURS=${1#*=} # remove everything up to = and assign rest to life
            shift
            ;;
        -d|--idle)
            if [ $# -ge 2 ]; then
                CLUSTER_MAX_IDLE_MINUTES="$2"
                shift 2
            else
                throw_error "--idle requires a non-empty argument"
            fi
            ;;
        --idle=?*)
            CLUSTER_MAX_IDLE_MINUTES=${1#*=} # remove everything up to = and assign rest to idle
            shift
            ;;
        --L|--launch-tool)
            if [ $# -ge 2 ]; then
                GATK_SV_TOOL="$2"
                shift 2
            else
                throw_error "--launch-tool requires a non-empty argument"
            fi
            ;;
        --launch-tool=?*)
            GATK_SV_TOOL=${1#*=} # remove everything up to = and assign rest to GATK_SV_TOOL
            shift
            ;;
        --)   # explicit call to end of all options
            shift
            break
            ;;
        -?*)  # unsupported option
            throw_error "Unknown option \"$1\". use --help for syntax"
            ;;
        *)  # not an option, a positional argument. break out
            break
            ;;
    esac
done

if [ $# -lt 4 ]; then
  show_help
  exit 1
fi

# init variables that MUST be defined
GATK_DIR="$1"
PROJECT_NAME="$2"
GCS_BAM="$3"
GCS_REFERENCE_2BIT="$4"
shift $(($# < 4 ? $# : 4))
SV_ARGS=${*:-${SV_ARGS:-""}} && SV_ARGS=${SV_ARGS:+" ${SV_ARGS}"}

# add GATK SV scripts to PATH
PATH="${GATK_DIR}/scripts/sv:${PATH}"

########################################################################################################################
# upfront sanity checks (exit whenever these checks fail)

export GATK_DIR
export QUIET
source sanity_checks.sh

########################################################################################################################
# various paths
##############################################
# paths on google bucket
GCS_BAM_DIR="$(dirname ${GCS_BAM})"
GCS_REFERENCE_DIR="$(dirname ${GCS_REFERENCE_2BIT})"

# for now assume the reference image is just the reference fasta + ".img"
GCS_REFERENCE_IMAGE="$(echo "${GCS_REFERENCE_2BIT}" | sed 's/.2bit$/.fasta.img/')"
if [ "$(dirname ${GCS_REFERENCE_2BIT})" != "$(dirname ${GCS_REFERENCE_IMAGE})" ]; then
    echo "Reference and reference image must be in same folder"
    exit -1
fi

GCS_SAVE_PATH=${GCS_SAVE_PATH:-"${PROJECT_NAME}-${GCS_USER}"}

# configure caching .jar files
export GATK_GCS_STAGING=${GATK_GCS_STAGING:-"gs://${PROJECT_NAME}/${GCS_USER}/staging/"}

##############################################
# paths on cluster (HDFS, except for reference image which must be regular FS)
CLUSTER_BAM="/data/$(basename ${GCS_BAM})"
CLUSTER_REFERENCE_2BIT="/reference/$(basename ${GCS_REFERENCE_2BIT})"
CLUSTER_REFERENCE_IMAGE="/reference/$(basename ${GCS_REFERENCE_IMAGE})"

##############################################
# store local run log in this file (NOTE: can override by defining SV_LOCAL_LOG_FILE)
SANITIZED_BAM=$(basename "${CLUSTER_BAM}" | awk '{print tolower($0)}' | sed 's/[^a-z0-9]/-/g')
TMPDIR=${TMPDIR:-/tmp/} # $TMPDIR is defined by default on OSX, not linux
LOCAL_LOG_FILE=${SV_LOCAL_LOG_FILE:-"${TMPDIR}sv-discovery-${SANITIZED_BAM}.log"}

########################################################################################################################
# call create_cluster, using default_init
CLUSTER_MAX_LIFE_HOURS=${CLUSTER_MAX_LIFE_HOURS:-4h}
CLUSTER_MAX_IDLE_MINUTES=${CLUSTER_MAX_IDLE_MINUTES:-60m}
INIT_SCRIPT=${INIT_SCRIPT:-"${GATK_DIR}/scripts/sv/default_init.sh"}
# set cluster name (NOTE: can override by defining SV_CLUSTER_NAME) based on user and target bam file and branch name
SANITIZED_GIT_BRANCH="master"
if [[ "$GIT_BRANCH" != "master" ]]; then
    SANITIZED_GIT_BRANCH="feature"
fi

CLUSTER_NAME=${SV_CLUSTER_NAME:-"${GCS_USER}-${SANITIZED_BAM}-${SANITIZED_GIT_BRANCH}"}
echo "Using cluster name \"${CLUSTER_NAME}\""
while true; do
    echo "#############################################################" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    if [ "${QUIET}" == "Y" ]; then
        if gcloud dataproc clusters list --project=${PROJECT_NAME} | grep -q ${CLUSTER_NAME}; then
            yn="N"
        else
            yn="Y"
        fi
    else
        read -p "Create cluster? (yes/no/cancel)" yn
    fi
    case $yn in
        [Yy]*)  if [[ ${INIT_SCRIPT} == gs://* ]]; then
                    INIT_ARGS=${INIT_SCRIPT}
                else
                    INIT_ARGS="${INIT_SCRIPT} gs://${GCS_SAVE_PATH}/init/"
                fi

                echo "create_cluster.sh ${GATK_DIR} ${PROJECT_NAME} ${CLUSTER_NAME} ${CLUSTER_MAX_LIFE_HOURS} ${CLUSTER_MAX_IDLE_MINUTES} ${GCS_REFERENCE_DIR} ${GCS_BAM} ${INIT_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}" | tee -a ${LOCAL_LOG_FILE}
                create_cluster.sh ${GATK_DIR} ${PROJECT_NAME} ${CLUSTER_NAME} ${CLUSTER_MAX_LIFE_HOURS} ${CLUSTER_MAX_IDLE_MINUTES} ${GCS_REFERENCE_DIR} ${GCS_BAM} ${INIT_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}
                break
                ;;
        [Nn]*)  break
                ;;
        [Cc]*)  exit
                ;;
        *)      echo "Please answer yes, no, or cancel."
                ;;
    esac
done

########################################################################################################################
# call runWholePipeline
# set output directory to datetime-git branch-git hash stamped folder
OUTPUT_DIR="/results/${SV_OUTPUT_DIR:-"$(date "+%Y-%m-%d_%H.%M.%S")-${SANITIZED_BAM}-${GIT_BRANCH}-${GATK_JAR_HASH}${UNTRACKED_COMMIT}"}"

while true; do
    echo "#############################################################" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    if [ "${QUIET}" == "Y" ]; then
      yn="Y"
    else
      read -p "Run whole pipeline? (yes/no/cancel)" yn
    fi
    case $yn in
        [Yy]*)  SECONDS=0
                echo "GATK_SV_TOOL=${GATK_SV_TOOL} run_whole_pipeline.sh ${GATK_DIR} ${PROJECT_NAME} ${CLUSTER_NAME} ${OUTPUT_DIR} ${CLUSTER_BAM} ${CLUSTER_REFERENCE_2BIT} ${CLUSTER_REFERENCE_IMAGE} ${SV_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}" | tee -a ${LOCAL_LOG_FILE}
                GATK_SV_TOOL=${GATK_SV_TOOL} run_whole_pipeline.sh ${GATK_DIR} ${PROJECT_NAME} ${CLUSTER_NAME} ${OUTPUT_DIR} ${CLUSTER_BAM} ${CLUSTER_REFERENCE_2BIT} ${CLUSTER_REFERENCE_IMAGE} ${SV_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}
                printf 'Pipeline completed in %02dh:%02dm:%02ds\n' $((${SECONDS}/3600)) $((${SECONDS}%3600/60)) $((${SECONDS}%60))
                break
                ;;
        [Nn]*)  break
                ;;
        [Cc]*)  exit
                ;;
        *)      echo "Please answer yes, no, or cancel."
                ;;
    esac
done

########################################################################################################################
# copy results into gcloud
while true; do
    echo "#############################################################" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    if [ "${QUIET}" == "Y" ]; then
      yn="Y"
    else
      read -p "Copy results? (yes/no/cancel)" yn
    fi
    case $yn in
        [Yy]*)  echo "copy_sv_results.sh ${PROJECT_NAME} ${CLUSTER_NAME} ${OUTPUT_DIR} ${GCS_USER} ${GCS_SAVE_PATH} ${LOCAL_LOG_FILE} 2>&1" | tee -a ${LOCAL_LOG_FILE}
                COPY_FASTQ=${COPY_FASTQ} copy_sv_results.sh ${PROJECT_NAME} ${CLUSTER_NAME} ${OUTPUT_DIR} ${GCS_USER} ${GCS_SAVE_PATH} ${LOCAL_LOG_FILE}${SV_ARGS}
                break
                ;;
        [Nn]*)  break
                ;;
        [Cc]*)  exit
                ;;
        *)      echo "Please answer yes, no, or cancel."
                ;;
    esac
done

########################################################################################################################
# delete cluster
while true; do
    echo "#############################################################"
    if [ "${QUIET}" == "Y" ]; then
      yn="Y"
    else
      read -p "Delete cluster? (yes/no/cancel)" yn
    fi
    case $yn in
        [Yy]*)  echo "gcloud dataproc clusters delete ${CLUSTER_NAME} --project=${PROJECT_NAME} --async --quiet"
                gcloud dataproc clusters delete ${CLUSTER_NAME} --project=${PROJECT_NAME} --async --quiet
                echo
                break
                ;;
        [Nn]*)  break
                ;;
        [Cc]*)  exit
                ;;
        *)      echo "Please answer yes, no, or cancel."
                ;;
    esac
done
