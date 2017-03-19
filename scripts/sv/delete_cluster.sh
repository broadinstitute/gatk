#!/bin/bash

# This script deletes a Google Dataproc cluster used for running the GATK-SV pipeline.

set -eu
if [ "$#" -ne 1 ]; then
    echo "Please provide name of cluster to be shutdown"
    exit 1
fi

CLUSTER_NAME=$1

gcloud dataproc clusters delete "$CLUSTER_NAME" --async
