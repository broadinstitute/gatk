#!/bin/bash
#
# A script to generate the GATK tool WDLs for a specified GATK release, and publish them
# to the GATK tool WDL github repository broadinstitute/gatk-tool-wdls.
#
# Must be run from the root of a GATK clone.
#
# Usage: bash scripts/publish_gatk_tool_wdls.sh <gatk_version_tag>
#
# Eg.,   bash scripts/publish_gatk_tool_wdls.sh 4.1.9.0
#
# The specified GATK version will also be used to tag the WDLs themselves in
# broadinstitute/gatk-tool-wdls
#

if [ $# -ne 1 ]; then
  echo "Usage: $0 <gatk_version_tag>"
  exit 1
fi

GATK_VERSION="$1"
GATK_CLONE_ROOT=$( pwd )
GATK_WDL_REPO="git@github.com:broadinstitute/gatk-tool-wdls.git"
WDL_STAGING_AREA="wdl_staging_area"
WDL_GEN_DIR="build/docs/wdlGen"
DOCKSTORE_YML=".dockstore.yml"

function cleanup() {
  cd "${GATK_CLONE_ROOT}"
  rm -rf "${WDL_STAGING_AREA}"
}

function fatal_error() {
  echo "$1" 1>&2
  cleanup
  exit 1
}

function generate_dockstore_yml() {
  local WDL_DIR="$1"
  local YML_FILE="$2"

  printf "version: 1.2\n" > "${YML_FILE}"
  printf "workflows:\n" >> "${YML_FILE}"

  for wdl in ${WDL_DIR}/*.wdl; do
    local WDL_NAME=$( basename "${wdl}" .wdl )
    local WDL_JSON="${WDL_NAME}Inputs.json"

    printf "   - name: %s\n" "${WDL_NAME}" >> "${YML_FILE}"
    printf "     subclass: WDL\n" >> "${YML_FILE}"
    printf "     primaryDescriptorPath: %s\n" "/${WDL_DIR}/${WDL_NAME}.wdl" >> "${YML_FILE}"
    printf "     testParameterFiles:\n" >> "${YML_FILE}"
    printf "     -  %s\n" "/${WDL_DIR}/${WDL_JSON}" >> "${YML_FILE}"
  done
}

# Checkout the GATK version for which we are publishing WDLs
echo "$0: Checking out GATK version ${GATK_VERSION}"
git checkout -f "${GATK_VERSION}"
if [ $? -ne 0 ]; then
  fatal_error "Failed to checkout GATK version ${GATK_VERSION}"
fi

# Get a fresh clone of the GATK WDL repo, and delete all existing WDLs
echo "$0: Cloning ${GATK_WDL_REPO} into ${WDL_STAGING_AREA}"
cleanup
git clone "${GATK_WDL_REPO}" "${WDL_STAGING_AREA}"
if [ $? -ne 0 ]; then
  fatal_error "Failed to clone GATK WDL repo ${GATK_WDL_REPO}"
fi

rm -rf "${WDL_STAGING_AREA}/wdls"
mkdir "${WDL_STAGING_AREA}/wdls"

# Generate the GATK WDLs, copy into our clone of the WDL repo, and delete the test WDLs
echo "$0: Running GATK WDL generation"
./gradlew clean gatkWDLGen
if [ $? -ne 0 ]; then
  fatal_error "Unable to generate the GATK tool WDLs"
fi

echo "$0: Copying WDLs to staging area"
cp ${WDL_GEN_DIR}/* "${WDL_STAGING_AREA}/wdls"
rm -f ${WDL_STAGING_AREA}/wdls/*Test.wdl
rm -f ${WDL_STAGING_AREA}/wdls/*TestInputs.json
rm -f ${WDL_STAGING_AREA}/wdls/*.html

# Switch to the root of our WDL staging clone, and add the WDLs we just generated to git
cd "${WDL_STAGING_AREA}"
git add wdls/*

# Regenerate the .dockstore.yml file
echo "$0: Generating .dockstore.yml"
rm -f "${DOCKSTORE_YML}"
generate_dockstore_yml wdls "${DOCKSTORE_YML}"
git add "${DOCKSTORE_YML}"

# Commit, tag, and push the WDLs
echo "$0: Pushing to github"
git commit -a -m "Version ${GATK_VERSION} of the GATK tool WDLs"
git tag -a -f -m "Version ${GATK_VERSION} of the GATK tool WDLs" "${GATK_VERSION}"
git push -f --tags origin HEAD:master
if [ $? -ne 0 ]; then
  fatal_error "Unable to push WDLs to github repo ${GATK_WDL_REPO}"
fi

# Cleanup and exit
echo "Done! Successfully published GATK tool wdls for version ${GATK_VERSION}"
cd "${GATK_CLONE_ROOT}"
cleanup
git checkout -f master
exit 0
