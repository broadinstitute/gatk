#!/bin/bash
#Install gcloud
if [ ! -d $HOME/gcloud/google-cloud-sdk ]; then
    mkdir -p $HOME/gcloud &&
    wget https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz --directory-prefix=$HOME/gcloud &&
    cd $HOME/gcloud &&
    tar xzf google-cloud-sdk.tar.gz &&
    printf '\ny\n\ny\ny\n' | ./google-cloud-sdk/install.sh &&
    cd $TRAVIS_BUILD_DIR;
fi
