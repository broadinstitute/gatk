#!/bin/bash -l
set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

WORKING_DIR=/home/travis/build/broadinstitute

set -e
echo "Building docker image for CNN WDL tests (skipping unit tests)..."

#assume Dockerfile is in root
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR/gatk

# IMPORTANT: This code is duplicated in the cnv and M2 WDL test.
if [ ${TRAVIS_PULL_REQUEST} != false ]; then
  HASH_TO_USE=FETCH_HEAD
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/ -t ${TRAVIS_PULL_REQUEST};
  echo "using fetch head:"$HASH_TO_USE
else
  HASH_TO_USE=${TRAVIS_COMMIT}
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/;
  echo "using travis commit:"$HASH_TO_USE
fi
echo "Docker build done =========="
sleep 10

cd $WORKING_DIR/gatk/scripts/
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" cnn_variant_wdl/jsons/cnn_score_variants_travis.json >$WORKING_DIR/cnn_score_variants_travis.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" cnn_variant_wdl/jsons/cnn_score_variants_travis_1d.json >$WORKING_DIR/cnn_score_variants_travis_1d.json
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" cnn_variant_wdl/jsons/cram2filtered_travis.json >$WORKING_DIR/cram2filtered_travis.json
echo "JSON FILES (modified) ======="
cat $WORKING_DIR/cnn_score_variants_travis.json
cat $WORKING_DIR/cnn_score_variants_travis_1d.json
cat $WORKING_DIR/cram2filtered_travis.json
echo "=================="


echo "Running CNN Score Variants WDL through cromwell"
ln -fs $WORKING_DIR/gatk/scripts/cnn_variant_wdl/cnn_score_variants.wdl
cd $WORKING_DIR/gatk/scripts/cnn_variant_wdl/
java -jar $CROMWELL_JAR run cnn_score_variants.wdl -i $WORKING_DIR/cnn_score_variants_travis_1d.json
java -jar $CROMWELL_JAR run cnn_score_variants.wdl -i $WORKING_DIR/cnn_score_variants_travis.json
java -jar $CROMWELL_JAR run cram2filtered.wdl -i $WORKING_DIR/cram2filtered_travis.json
