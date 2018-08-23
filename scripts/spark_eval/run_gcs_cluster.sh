#!/usr/bin/env bash

# Starts a GCS cluster, runs scripts, then deletes the cluster.

if [ -z "$GCS_CLUSTER" ]; then
  echo "Please set the GCS_CLUSTER environment variable to the name of the cluster you would like to start."
  exit 1
fi

# Create cluster
gcloud beta dataproc clusters create "$GCS_CLUSTER" \
  --zone us-central1-a \
  --master-machine-type n1-standard-4 \
  --master-boot-disk-size 500 \
  --num-workers ${NUM_WORKERS:-10} \
  --worker-machine-type n1-standard-16 \
  --worker-boot-disk-size 2000 \
  --image-version 1.2 \
  --max-age 3h \
  --project broad-gatk-collab

# Run scripts
for script in "$@"
do
  SCRIPT_NAME="$script"
  source "$script"
done

# Delete cluster
gcloud dataproc clusters delete --quiet "$GCS_CLUSTER"