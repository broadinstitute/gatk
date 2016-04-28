#!/usr/bin/env bash

#Helper script to be called from upload_artifact_with_both_lib_and_dylib.sh
#
#checkout and clone a new copy of gatk in a tmpdir
#build the .so file and return its location
#this can only be run on broad machines
#
# usage build_pair_hmm_so_in_clean_repo.sh <commit hash>

set -e
set -v

COMMIT=$1

VECTOR_LIB="libVectorLoglessPairHMM.so"
LIB_PATH="build/libs/vectorLoglessPairHMM/shared"

export TMPDIR="/broad/hptmp"

GATK_TMP_DIR=`mktemp --tmpdir -d 2>/dev/null || mktemp -d -t 'mytmpdir'`
function finish {
  rm -rf "$GATK_TMP_DIR/gatk"
}
trap finish EXIT

cd "$GATK_TMP_DIR"
GIT_LFS_SKIP_SMUDGE=1 git clone git@github.com:broadinstitute/gatk.git 1>2
cd "gatk"
GIT_LFS_SKIP_SMUDGE=1 git checkout -f "$COMMIT" 1>2

./gradlew copySharedLib 1>2
cp "$LIB_PATH/$VECTOR_LIB" "$GATK_TMP_DIR/$VECTOR_LIB"

echo "$GATK_TMP_DIR/$VECTOR_LIB"
exit 0
