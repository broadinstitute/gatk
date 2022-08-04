#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

WORKING_DIR=/home/runner/work/gatk

set -e

echo "Creating tar.gz for Funcotator datasources =========="
pushd .
FUNCOTATOR_TEST_DS_DIR=${WORKING_DIR}/gatk/src/test/resources/large/funcotator/
cd ${FUNCOTATOR_TEST_DS_DIR}
# First parameter must match Mutect2_Multi.funco_data_sources_tar_gz test_m2_wdl.json
tar zcvf ${WORKING_DIR}/gatk/small_ds_pik3ca.tar.gz small_ds_pik3ca/*
popd

echo "Building docker image for M2 WDL tests (skipping unit tests)..."

#cd $WORKING_DIR/gatk/scripts/docker/
#assume Dockerfile is in root
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR/gatk
# IMPORTANT: This code is duplicated in the cnv WDL test.
if [ ! -z "$CI_PULL_REQUEST" ]; then
  HASH_TO_USE=FETCH_HEAD
  echo "Building pr build with $HASH_TO_USE... ========="
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/ -t ${CI_PULL_REQUEST};
else
  HASH_TO_USE=${CI_COMMIT}
  echo "Building commit build with $HASH_TO_USE... ========="
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/;
fi
echo "Docker build done =========="
echo "Putting the newly built docker image into the json parameters"
cd $WORKING_DIR/gatk/scripts/
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" m2_cromwell_tests/test_m2_wdl.json >$WORKING_DIR/test_m2_wdl_mod.json
echo "JSON FILE (modified) ======="
cat $WORKING_DIR/test_m2_wdl_multi_mod.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" m2_cromwell_tests/test_mitochondria_m2_wdl.json >$WORKING_DIR/test_mitochondria_m2_wdl_mod.json
echo "JSON FILE (modified) ======="
cat $WORKING_DIR/test_mitochondria_m2_wdl_mod.json
echo "=================="

# Create the tumor-only json by removing normal_reads and normal_reads_index from the input json
grep -v 'Mutect2.normal_reads' $WORKING_DIR/test_m2_wdl_mod.json >$WORKING_DIR/test_m2_wdl_mod_to.json
cd $WORKING_DIR/

echo "Running M2 WDL through cromwell (T/N)"
ln -fs $WORKING_DIR/gatk/scripts/mutect2_wdl/mutect2.wdl
sudo java -jar $CROMWELL_JAR run $WORKING_DIR/gatk/scripts/mutect2_wdl/mutect2.wdl -i $WORKING_DIR/test_m2_wdl_mod.json -m $WORKING_DIR/test_m2_wdl.metadata

echo "Running M2 WDL through cromwell (Tumor-only)"
sudo java -jar $CROMWELL_JAR run $WORKING_DIR/gatk/scripts/mutect2_wdl/mutect2.wdl -i $WORKING_DIR/test_m2_wdl_mod_to.json -m $WORKING_DIR/test_m2_wdl_to.metadata

echo "Running Mitochondria M2 WDL through cromwell"
ln -fs $WORKING_DIR/gatk/scripts/mitochondria_m2_wdl/AlignAndCall.wdl
ln -fs $WORKING_DIR/gatk/scripts/mitochondria_m2_wdl/AlignmentPipeline.wdl
sudo java -jar $CROMWELL_JAR run $WORKING_DIR/gatk/scripts/mitochondria_m2_wdl/MitochondriaPipeline.wdl -i $WORKING_DIR/test_mitochondria_m2_wdl_mod.json -m $WORKING_DIR/test_mitochondria_m2_wdl.metadata
