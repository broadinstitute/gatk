#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

WORKING_DIR=/home/runner/work/gatk/gatk
WDL_DIR=$WORKING_DIR/scripts/vcf_site_level_filtering_wdl
CROMWELL_TEST_DIR=$WORKING_DIR/scripts/vcf_site_level_filtering_cromwell_tests

set -e
echo "Building docker image for VCF Site Level Filtering WDL tests (skipping unit tests)..."

#assume Dockerfile is in root
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR

# IMPORTANT: This code is duplicated in the cnv and M2 WDL test.
if [ ! -z "$CI_PULL_REQUEST" ]; then
  HASH_TO_USE=FETCH_HEAD
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/ -t ${CI_PULL_REQUEST};
  echo "using fetch head:"$HASH_TO_USE
else
  HASH_TO_USE=${CI_COMMIT}
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/;
  echo "using commit:"$HASH_TO_USE
fi
echo "Docker build done =========="

sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" $CROMWELL_TEST_DIR/vcf_site_level_filtering.json >$WORKING_DIR/vcf_site_level_filtering_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" $CROMWELL_TEST_DIR/vcf_site_level_filtering_bgmm.json >$WORKING_DIR/vcf_site_level_filtering_bgmm_mod.json

echo "Running Filtering WDL through cromwell"

cat $WORKING_DIR/vcf_site_level_filtering_mod.json
java -jar $CROMWELL_JAR run $WDL_DIR/JointVcfFiltering.wdl -i $WORKING_DIR/vcf_site_level_filtering_mod.json

cat $WORKING_DIR/vcf_site_level_filtering_bgmm_mod.json
java -jar $CROMWELL_JAR run $WDL_DIR/JointVcfFiltering.wdl -i $WORKING_DIR/vcf_site_level_filtering_bgmm_mod.json
