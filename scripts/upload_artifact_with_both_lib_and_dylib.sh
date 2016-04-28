#!/usr/bin/env bash

#Upload a snapshot build of the current HEAD with both .dylib and .so
#This can only be run from within Broad
#this is a temporary workaround for #1645 until we have a solution for cross compiling a dylib and an so on the same machine

set -v
set -e

SERVER=gsa6.broadinstitute.org

LIB_PATH="build/classes/main/lib/"

./gradlew clean copySharedLib
vectorLib=$( ssh $SERVER 'bash -s' < scripts/build_pair_hmm_so_in_clean_repo.sh $( git rev-parse HEAD) )
scp ${SERVER}:${vectorLib} $LIB_PATH

./gradlew uploadArchives

exit 0
