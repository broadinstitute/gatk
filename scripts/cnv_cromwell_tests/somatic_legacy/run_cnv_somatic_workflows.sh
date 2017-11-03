#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_common_tasks.wdl
ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_copy_ratio_bam_workflow.wdl
ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_allele_fraction_pair_workflow.wdl
ln -fs /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_oncotate.wdl

WORKING_DIR=/home/travis/build/broadinstitute

pushd .
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR/gatk
# IMPORTANT: This code is duplicated in the M2 WDL test.
if [ ${TRAVIS_PULL_REQUEST} != false ]; then
  HASH_TO_USE=FETCH_HEAD
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/ -t ${TRAVIS_PULL_REQUEST};
else
  HASH_TO_USE=${TRAVIS_COMMIT}
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/;
fi
echo "Docker build done =========="

popd

echo "Inserting docker image into json ========"
CNV_CROMWELL_TEST_DIR="${WORKING_DIR}/gatk/scripts/cnv_cromwell_tests/somatic_legacy/"
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wes_tumor-only_workflow.json >cnv_somatic_pair_wes_tumor-only_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wgs_tumor-only_workflow.json >cnv_somatic_pair_wgs_tumor-only_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wes_workflow.json >cnv_somatic_pair_wes_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wgs_workflow.json >cnv_somatic_pair_wgs_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wes_workflow.json >cnv_somatic_panel_wes_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wgs_workflow.json >cnv_somatic_panel_wgs_workflow_mod.json


echo "Running ========"
CROMWELL_JAR="cromwell-0.28.jar"

# Panel WES
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wes_workflow_mod.json
# Panel WGS
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_panel_workflow.wdl cnv_somatic_panel_wgs_workflow_mod.json

# Pair WES
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wes_workflow_mod.json
# Pair WGS
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wgs_workflow_mod.json
# Pair WES tumor-only
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wes_tumor-only_workflow_mod.json
# Pair WGS tumor-only
java -jar ~/${CROMWELL_JAR} run /home/travis/build/broadinstitute/gatk/scripts/cnv_wdl/somatic_legacy/cnv_somatic_pair_workflow.wdl cnv_somatic_pair_wgs_tumor-only_workflow_mod.json