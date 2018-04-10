#!/bin/bash

# This script deletes a Google Dataproc cluster used for running the GATK-SV pipeline.

set -eu
if [ "$#" -ne 2 ]; then
    echo "Usage: delete_cluster.sh <CLUSTER_NAME> <PROJECT>"
    exit 1
fi

CLUSTER_NAME=$1
PROJECT=$2

gcloud dataproc clusters delete --project $PROJECT "$CLUSTER_NAME" --async
