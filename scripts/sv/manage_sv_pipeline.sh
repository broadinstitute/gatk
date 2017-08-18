#!/bin/bash

# terminate script on error or if a command fails before piping to another command
set -eu
set -o pipefail

show_help() {
cat << EOF
Manage the SV discovery pipeline on Google Cloud Services (GCS) cluster
Syntax
  manage_sv_pipeline.sh [Options] GATK_Folder Project_Name \\
                        GCS_BAM_File GCS_FASTA_File \\
                        [args to SV discovery pipeline]
  Options:
    -h or --help:
      display help and exit
    -q or --quiet:
      run without prompting for user input
    -u [GCS_USER] or --user [GCS_USER]:
      set GCS username (defaults to ${USER})
    -i [INIT_SCRIPT] or --init [INIT_SCRIPT]
      use specified spark init script
      (defaults to \${GATK_FOLDER}/scripts/sv/default_init.sh)
    -s [GCS_SAVE_PATH] or --save [GCS_SAVE_PATH]
      save results in specified bucket/folder
      (defaults to \$PROJECT_NAME/\$GCS_USER)

  Mandatory Positional Arguments:
    GATK_Folder: path to local copy of GATK
    Project_Name: name of GCS project
    GCS_BAM_File: path to .bam file hosted in GCS
                  (e.g. gs://bucket/path/to/file.bam)
    GCS_FASTA_File: path to reference .fasta file hosted in GCS

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
f) shut down the cluster
EOF
}

throw_error() {
    echo "$1" >&2
    exit 1
}

QUIET=${SV_QUIET:-"N"}
GCS_USER=${GCS_USER:-${USER}}

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
GCS_REFERENCE_FASTA="$4"
shift $(($# < 4 ? $# : 4))
SV_ARGS=${*:-${SV_ARGS:-""}} && SV_ARGS=${SV_ARGS:+" ${SV_ARGS}"}

INIT_SCRIPT=${INIT_SCRIPT:-"${GATK_DIR}/scripts/sv/default_init.sh"}
GCS_SAVE_PATH=${GCS_SAVE_PATH:-"${PROJECT_NAME}/${GCS_USER}"}

# add GATK SV scripts to PATH
PATH="${GATK_DIR}/scripts/sv:${PATH}"

# configure caching .jar files
export GATK_GCS_STAGING=${GATK_GCS_STAGING:-"gs://${PROJECT_NAME}/${GCS_USER}/staging/"}

# set cluster name based on user and target bam file
# (NOTE: can override by defining SV_CLUSTER_NAME)
SANITIZED_BAM=$(basename "${GCS_BAM}" | awk '{print tolower($0)}' | sed 's/[^a-z0-9]/-/g')
CLUSTER_NAME=${SV_CLUSTER_NAME:-"${GCS_USER}-${SANITIZED_BAM}"}
echo "Using cluster name \"${CLUSTER_NAME}\""

# update gcloud
if [ "${QUIET}" == "Y" ]; then
    gcloud components update --quiet
else
    gcloud components update
fi

# for now assume the reference image is just the reference fasta + ".img"
GCS_REFERENCE_IMAGE="${GCS_REFERENCE_FASTA}.img"
# monkey around with variables to put them in form used by gatk/sv/scripts
GCS_BAM_DIR="$(dirname ${GCS_BAM})"
GCS_BAM="/data/$(basename ${GCS_BAM})"
GCS_REFERENCE_DIR="$(dirname ${GCS_REFERENCE_FASTA})"
if [ "$(dirname ${GCS_REFERENCE_IMAGE})" != "$(dirname ${GCS_REFERENCE_IMAGE})" ]; then
    echo "Reference fasta and reference image must be in same folder"
    exit -1
fi
GCS_REFERENCE_FASTA="/reference/$(basename ${GCS_REFERENCE_FASTA})"
GCS_REFERENCE_IMAGE="/mnt/1/reference/$(basename ${GCS_REFERENCE_IMAGE})"

# store run log in this file
# (NOTE: can override by defining SV_LOCAL_LOG_FILE)
LOCAL_LOG_FILE=${SV_LOCAL_LOG_FILE:-"${TMPDIR}sv-discovery-${SANITIZED_BAM}.log"}

# check if GATK jar was compiled from the current .git hash
GATK_GIT_HASH=$(readlink ${GATK_DIR}/build/libs/gatk-spark.jar | cut -d- -f5 | cut -c2-)
CURRENT_GIT_HASH=$(git -C ${GATK_DIR} rev-parse --short HEAD | cut -c1-7)
if [ "${QUIET}" != "Y" ] && [ "${GATK_GIT_HASH}" != "${CURRENT_GIT_HASH}" ]; then
    while true; do
        read -p "Current git hash does not match GATK git hash. Run anyway?" yn
        case $yn in
            [Yy]*)  break
                    ;;
            [Nn]*)  exit
                    ;;
            *)      echo "Please answer yes or no"
                    ;;
        esac
    done
fi
GIT_BRANCH=$(git -C ${GATK_DIR} branch --contains ${GATK_GIT_HASH} | rev | cut -d" " -f 1 | rev)
# set output directory to datetime-git branch-git hash stamped folder
OUTPUT_DIR="/results/$(date "+%Y-%m-%d_%H.%M.%S")-${GIT_BRANCH}-${GATK_GIT_HASH}"

# call create_cluster, using default_init
while true; do
    echo "#############################################################" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    if [ "${QUIET}" == "Y" ]; then
        yn="Y"
    else
        read -p "Create cluster? (yes/no/cancel)" yn
    fi
    case $yn in
        [Yy]*)  if [[ ${INIT_SCRIPT} == gs://* ]]; then
                    INIT_ARGS=${INIT_SCRIPT}
                else
                    INIT_ARGS="${INIT_SCRIPT} gs://${GCS_SAVE_PATH}/init/"
                fi

                echo "create_cluster.sh ${GATK_DIR} ${PROJECT_NAME} ${CLUSTER_NAME} ${GCS_REFERENCE_DIR} ${GCS_BAM_DIR} ${INIT_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}" | tee -a ${LOCAL_LOG_FILE}
                create_cluster.sh ${GATK_DIR} ${PROJECT_NAME} ${CLUSTER_NAME} ${GCS_REFERENCE_DIR} ${GCS_BAM_DIR} ${INIT_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}
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

# call runWholePipeline
while true; do
    echo "#############################################################" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    if [ "${QUIET}" == "Y" ]; then
      yn="Y"
    else
      read -p "Run whole pipeline? (yes/no/cancel)" yn
    fi
    case $yn in
        [Yy]*)  SECONDS=0
                echo "runWholePipeline.sh ${GATK_DIR} ${CLUSTER_NAME} ${OUTPUT_DIR} ${GCS_BAM} ${GCS_REFERENCE_FASTA} ${GCS_REFERENCE_IMAGE}${SV_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}" | tee -a ${LOCAL_LOG_FILE}
                runWholePipeline.sh ${GATK_DIR} ${CLUSTER_NAME} ${OUTPUT_DIR} ${GCS_BAM} ${GCS_REFERENCE_FASTA} ${GCS_REFERENCE_IMAGE}${SV_ARGS} 2>&1 | tee -a ${LOCAL_LOG_FILE}
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
                copy_sv_results.sh ${PROJECT_NAME} ${CLUSTER_NAME} ${OUTPUT_DIR} ${GCS_USER} ${GCS_SAVE_PATH} ${LOCAL_LOG_FILE}${SV_ARGS}
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
