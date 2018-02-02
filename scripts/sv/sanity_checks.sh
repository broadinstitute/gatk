#!/usr/bin/env bash

set -eu

# sanity 0: must have jar
echo
if [[ ! -f ${GATK_DIR}/build/libs/gatk-spark.jar ]]; then
    echo "Cannot find GATK spark jar, maybe you forgot to build? Given GATK dir.: ${GATK_DIR}"
    exit 1
fi

# sanity 1: check if there's un-committed changes
# there's the possibility that the jar to be used was built with un-committed changes
# but tagged with a git hash that is present in multiple branches,
# this makes $GIT_BRANCH a space-delimited string hence may mess with arg parsing
GIT_BRANCH=$(git -C ${GATK_DIR} rev-parse --abbrev-ref HEAD)
UNTRACKED_COMMIT=$(git diff-index --quiet HEAD -- || echo "-untracked")
if [[ ${UNTRACKED_COMMIT} ]]; then
    echo "There are uncommitted changes in your current branch ($GIT_BRANCH); this might cause problems:"
    echo "  either later argument parsing may be messed up, or "
    echo "  misleading result-tracking--GATK jar could have been built with uncommitted changes but is stamped with the latest commit hash"
    if [[ "${QUIET}" == "Y" ]]; then
        echo "ERROR: Quitting because running with uncommitted changes leads to inconsistent result-tracking"
        exit 1
    fi
    read -p "Want to proceed anyway? (yes/no/cancel)" yn
    case $yn in
        [Yy]*)  ;;
        [Nn]*)  exit
                ;;
        [Cc]*)  exit
                ;;
            *)  echo "Please answer yes, no, or cancel."
                ;;
    esac
fi

# sanity 2: check if GATK jar was compiled from the current .git hash
GATK_JAR_HASH=$(readlink ${GATK_DIR}/build/libs/gatk-spark.jar | awk 'BEGIN {FS="-g"} {print $2}' | cut -d- -f 1)
CURRENT_GIT_HASH=$(git -C ${GATK_DIR} rev-parse --short HEAD | cut -c1-7)
if [[ "${QUIET}" != "Y" ]] && [[ "${GATK_JAR_HASH}" != "${CURRENT_GIT_HASH}" ]]; then
    while true; do
        echo "gatk-spark.jar version (${GATK_JAR_HASH}) does not match current git commit (${CURRENT_GIT_HASH})."
        read -p  "Run anyway? (yes/no)" yn
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
