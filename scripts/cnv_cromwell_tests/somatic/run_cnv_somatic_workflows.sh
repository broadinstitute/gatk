#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

WORKING_DIR=/home/runner/work/gatk

ln -fs $WORKING_DIR/scripts/cnv_wdl/cnv_common_tasks.wdl
ln -fs $WORKING_DIR/scripts/cnv_wdl/somatic/cnv_somatic_oncotator_workflow.wdl

pushd .
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR/gatk
# IMPORTANT: This code is duplicated in the M2 WDL test.
if [ ! -z "$CI_PULL_REQUEST" ]; then
  HASH_TO_USE=FETCH_HEAD
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/ -t ${CI_PULL_REQUEST};
else
  HASH_TO_USE=${CI_COMMIT}
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/;
fi
echo "Docker build done =========="
sleep 10

popd

echo "Inserting docker image into json ========"
CNV_CROMWELL_TEST_DIR="${WORKING_DIR}/gatk/scripts/cnv_cromwell_tests/somatic/"
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wes_no-gc_workflow.json >cnv_somatic_panel_wes_no-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wgs_no-gc_workflow.json >cnv_somatic_panel_wgs_no-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wes_do-gc_workflow.json >cnv_somatic_panel_wes_do-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_panel_wgs_do-gc_workflow.json >cnv_somatic_panel_wgs_do-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wes_no-gc_workflow.json >cnv_somatic_pair_wes_no-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wgs_no-gc_workflow.json >cnv_somatic_pair_wgs_no-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wes_do-gc_workflow.json >cnv_somatic_pair_wes_do-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wgs_do-gc_workflow.json >cnv_somatic_pair_wgs_do-gc_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wgs_do-gc_tumor_only_workflow.json > cnv_somatic_pair_wgs_do-gc_tumor_only_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wes_no-gc_tumor_only_workflow.json > cnv_somatic_pair_wes_no-gc_tumor_only_workflow_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" ${CNV_CROMWELL_TEST_DIR}/cnv_somatic_pair_wes_do-gc_workflow_funcotator.json > cnv_somatic_pair_wes_do-gc_workflow_funcotator_mod.json

echo "Running ========"

# Panel WES
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl -i cnv_somatic_panel_wes_no-gc_workflow_mod.json
# Panel WGS
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl -i cnv_somatic_panel_wgs_no-gc_workflow_mod.json
# Panel WES w/ explicit GC correction
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl -i cnv_somatic_panel_wes_do-gc_workflow_mod.json
# Panel WGS w/ explicit GC correction
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_panel_workflow.wdl -i cnv_somatic_panel_wgs_do-gc_workflow_mod.json

# Pair WES
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl -i cnv_somatic_pair_wes_no-gc_workflow_mod.json
# Pair WGS
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl -i cnv_somatic_pair_wgs_no-gc_workflow_mod.json
# Pair WES w/ explicit GC correction
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl -i cnv_somatic_pair_wes_do-gc_workflow_mod.json
# Pair WGS w/ explicit GC correction
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl -i cnv_somatic_pair_wgs_do-gc_workflow_mod.json
# Tumor only WGS w/ explicit GC correction
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl -i cnv_somatic_pair_wgs_do-gc_tumor_only_workflow_mod.json
# Tumor only WES w/o explicit GC correction
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl -i cnv_somatic_pair_wes_no-gc_tumor_only_workflow_mod.json

#### Test FuncotateSegments on a small run.
#  Must tar up test datasources for the WDL
echo "Preparing small gencode-only datasource dir for funcotator test..."
pushd .
cd $WORKING_DIR/gatk/src/test/resources/org/broadinstitute/hellbender/tools/funcotator/
tar zcvf $WORKING_DIR/gatk/ds.tar.gz small_cntn4_ds/gencode_cntn4/
popd

# Pair WES with funcotator
java -jar ${CROMWELL_JAR} run $WORKING_DIR/gatk/scripts/cnv_wdl/somatic/cnv_somatic_pair_workflow.wdl \
  -i cnv_somatic_pair_wes_do-gc_workflow_funcotator_mod.json
#####
