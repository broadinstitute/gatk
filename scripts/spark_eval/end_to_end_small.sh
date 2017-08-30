#!/usr/bin/env bash

# Create cluster
gcloud dataproc clusters create "$GCS_CLUSTER" \
  --zone us-central1-a \
  --master-machine-type n1-standard-4 \
  --master-boot-disk-size 500 \
  --num-workers 5 \
  --worker-machine-type n1-standard-16 \
  --worker-boot-disk-size 2000 \
  --image-version 1.1 \
  --project broad-gatk-collab

# Copy data prep script to cluster, then run it
gcloud compute scp prep_data_small_gcs.sh "$GCS_CLUSTER"-m:prep_data_small_gcs.sh \
  --zone us-central1-a
gcloud compute ssh "$GCS_CLUSTER"-m \
  --command "./prep_data_small_gcs.sh /user/$USER/small_spark_eval" \
  --zone us-central1-a

# Run pipeline
source small_pipeline_single_gcs.sh
#source small_pipeline_single_gcs_hdfs.sh

gsutil ls -l gs://gatk-tom-testdata-small/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.vcf

# Delete cluster
#gcloud dataproc clusters delete "$GCS_CLUSTER"