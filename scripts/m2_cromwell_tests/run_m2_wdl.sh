#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

WORKING_DIR=/home/travis/build/broadinstitute

set -e
echo "Building docker image for M2 WDL tests (skipping unit tests)..."
HASH_TO_USE=`git rev-parse ${TRAVIS_BRANCH}`
#cd $WORKING_DIR/gatk-protected/scripts/docker/
#assume Dockerfile is in root
cp -rfp $WORKING_DIR/gatk-protected/scripts/docker/Dockerfile $WORKING_DIR/gatk-protected/
cp -rfp $WORKING_DIR/gatk-protected/scripts/docker/build_docker.sh $WORKING_DIR/gatk-protected/
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR/gatk-protected
sudo bash build_docker.sh  -e $HASH_TO_USE -s -u
echo "Docker build done =========="
echo "Putting the newly built docker image into the json parameters"
cd $WORKING_DIR/gatk-protected/scripts/docker/
sed -r "s/__M2_DOCKER__/broadinstitute\/gatk-protected\:$HASH_TO_USE/g" ../m2_cromwell_tests/test_m2_wdl_multi.json >$WORKING_DIR/test_m2_wdl_multi_mod.json
echo "JSON FILE (modified) ======="
cat $WORKING_DIR/test_m2_wdl_multi_mod.json
echo "=================="

# Create the tumor-only json by using the pair_list_tumor_only file
sed -r "s/\"pair_list/\"pair_list_tumor_only/g" $WORKING_DIR/test_m2_wdl_multi_mod.json >$WORKING_DIR/test_m2_wdl_multi_mod_to.json
cd $WORKING_DIR/

echo "Running M2 WDL through cromwell (T/N)"
ln -fs $WORKING_DIR/gatk-protected/scripts/mutect2_wdl/mutect2.wdl
sudo java -jar ~/cromwell-0.25.jar run $WORKING_DIR/gatk-protected/scripts/mutect2_wdl/mutect2_multi_sample.wdl $WORKING_DIR/test_m2_wdl_multi_mod.json - $WORKING_DIR/test_m2_wdl.metadata

echo "Running M2 WDL through cromwell (Tumor-only)"
sudo java -jar ~/cromwell-0.25.jar run $WORKING_DIR/gatk-protected/scripts/mutect2_wdl/mutect2_multi_sample.wdl $WORKING_DIR/test_m2_wdl_multi_mod_to.json - $WORKING_DIR/test_m2_wdl_to.metadata
