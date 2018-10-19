#!/usr/bin/env bash
# Install and Gcloud and authenticate on travis

export BOTO_CONFIG=/dev/null;
openssl aes-256-cbc -K $encrypted_703d76169d63_key -iv $encrypted_703d76169d63_iv -in resources_for_CI/servicekey.json.enc -out servicekey.json -d;
scripts/install_gcloud.sh;
printf 'y\n' | $GCLOUD_HOME/gcloud components update;
if [[ $TEST_TYPE == cloud ]]; then
   printf 'y\n' | $GCLOUD_HOME/gcloud components install beta;
fi;
$GCLOUD_HOME/gcloud config set project broad-dsde-dev;
$GCLOUD_HOME/gcloud auth activate-service-account --key-file servicekey.json;
$GCLOUD_HOME/gcloud --version
