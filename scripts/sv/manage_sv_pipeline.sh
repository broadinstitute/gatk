#!/bin/bash

# This script manages the SV discovery pipeline by calling other scripts
# a) it creates a Google Dataproc cluster used for running the GATK-SV
#    pipeline (call create_cluster.sh)
# b) it runs the analysis (call runWholePipeline.sh)
# c) it copies results to a datetime and git-version stamped directory
#    on GCS, along with logged console output (call copy_sv_results.sh)
# d) it shuts down the cluster

# terminate script on error or if a command fails before piping to another command
set -eu
set -o pipefail

if [[ "$#" -lt 4 ]]; then
    echo -e \
"Please provide:
  [1] local directory of GATK build (required)
  [2] project name (required)
  [3] GCS path to indexed BAM (required)
  [4] GCS path to reference fasta (required)
      OPTIONAL arguments:
  [5] path to initialization script (local or GCS, defaults to
      \${GATK_DIR}/scripts/sv/default_init.sh if omitted or empty)
  [6] GCS username (defaults to local username if omitted or empty)
  [*] additional arguments to pass to
      StructuralVariationDiscoveryPipelineSpark
To leave a value as default but specify a later value, use an empty
  string. e.g. to use default initialization script but override
  GCS user name:
\$ manage_sv_pipeline.sh /my/gatk/dir my_project gs://mybam.bam \\
     gs://myfasta.fasta \"\" my_gcs_user"
    exit 1
fi

# init variables that MUST be defined
GATK_DIR="$1"
PROJECT_NAME="$2"
GCS_BAM="$3"
GCS_REFERENCE_FASTA="$4"

# init variables with optional overrides / defaults
# passed value -> overrides system value -> overrides default value
INIT_SCRIPT=${5:-${INIT_SCRIPT:-"${GATK_DIR}/scripts/sv/default_init.sh"}}
GCS_USER=${6:-${GCS_USER:-${USER}}}
shift 6
SV_ARGS=${*:-${SV_ARGS:-""}} && SV_ARGS=${SV_ARGS:+" ${SV_ARGS}"}

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
gcloud components update

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
if [ "${GATK_GIT_HASH}" != "${CURRENT_GIT_HASH}" ]; then
        while true; do
                read -p "Current git hash does not match GATK git hash. Run anyway?" yn
                case $yn in
                        [Yy]*)  break
                                ;;
                        [Nn]*)  exit
                                ;;
                        *)      echo "Please answer yes or no"
                        
                esac
        done
fi
GIT_BRANCH=$(git -C ${GATK_DIR} branch --contains ${GATK_GIT_HASH} | rev | cut -d" " -f 1 | rev)
# set output directory to datetime-git branch-git hash stamped folder
OUTPUT_DIR="/results/$(date "+%Y-%m-%d_%H.%M.%S")-${GIT_BRANCH}-${GATK_GIT_HASH}"

# call create_cluster, using default_init
while true; do
    echo "#############################################################" 2>&1 | tee -a ${LOCAL_LOG_FILE}
    read -p "Create cluster? (yes/no/cancel)" yn
    case $yn in
        [Yy]*)  if [[ ${INIT_SCRIPT} == gs://* ]]; then
                        INIT_ARGS=${INIT_SCRIPT}
                else
                        INIT_ARGS="${INIT_SCRIPT} gs://${PROJECT_NAME}/${GCS_USER}/init/$(basename INIT_SCRIPT)"
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
    read -p "Run whole pipeline? (yes/no/cancel)" yn
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
    read -p "Copy results? (yes/no/cancel)" yn
    case $yn in
        [Yy]*)  echo "copy_sv_results.sh ${PROJECT_NAME} ${CLUSTER_NAME} ${OUTPUT_DIR} ${GCS_USER} ${LOCAL_LOG_FILE} 2>&1" | tee -a ${LOCAL_LOG_FILE}
                copy_sv_results.sh ${PROJECT_NAME} ${CLUSTER_NAME} ${OUTPUT_DIR} ${GCS_USER} ${LOCAL_LOG_FILE}${SV_ARGS}
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
    read -p "Delete cluster? (yes/no/cancel)" yn
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
