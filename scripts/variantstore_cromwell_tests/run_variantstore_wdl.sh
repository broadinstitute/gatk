#!/usr/bin/env bash

set -e
#cd in the directory of the script in order to use relative paths
script_path=$( cd "$(dirname "${BASH_SOURCE}")" ; pwd -P )
cd "$script_path"

WORKING_DIR=/home/travis/build/broadinstitute
UUID=$(cat /proc/sys/kernel/random/uuid | sed s/-/_/g)

set -e

echo "Building docker image for VariantStore WDL tests (skipping unit tests)..."

#assume Dockerfile is in root
echo "Building docker without running unit tests... ========="
cd $WORKING_DIR/gatk
# IMPORTANT: This code is duplicated in the cnv WDL test.
if [ ${TRAVIS_PULL_REQUEST} != false ]; then
  HASH_TO_USE=FETCH_HEAD
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/ -t ${TRAVIS_PULL_REQUEST};
else
  HASH_TO_USE=${TRAVIS_COMMIT}
  sudo bash build_docker.sh  -e ${HASH_TO_USE} -s -u -d $PWD/temp_staging/;
fi
echo "Docker build done =========="
echo "Putting the newly built docker image into the json parameters"
cd $WORKING_DIR/gatk/scripts/
sed -r "s/__GATK_DOCKER__/broadinstitute\/gatk\:$HASH_TO_USE/g" variantstore_cromwell_tests/import_array_manifest_test.json >$WORKING_DIR/import_array_manifest_test_tmp.json
sed -r "s/__TABLE_NAME__/$UUID/g" $WORKING_DIR/import_array_manifest_test_tmp.json > $WORKING_DIR/import_array_manifest_test_mod.json
echo "JSON FILE (modified) ======="
cat $WORKING_DIR/import_array_manifest_test_mod.json

sed -r "s|__SERVICE_ACCOUNT__|$GOOGLE_APPLICATION_CREDENTIALS|g" variantstore_cromwell_tests/local-with-gcs.conf >$WORKING_DIR/set_up.conf
echo "Updated local_backend.conf with service account"

echo "Running ImportArrayManifest WDL through cromwell"
ln -fs $WORKING_DIR/gatk/scripts/variantstore_wdl/ImportArrayManifest.wdl
sudo java -Dconfig.file=$WORKING_DIR/set_up.conf -jar $CROMWELL_JAR run $WORKING_DIR/gatk/scripts/variantstore_wdl/ImportArrayManifest.wdl -i $WORKING_DIR/import_array_manifest_test_mod.json -m $WORKING_DIR/test_import_manifest_wdl.metadata
