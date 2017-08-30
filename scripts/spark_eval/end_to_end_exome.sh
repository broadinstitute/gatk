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
gcloud compute scp prep_data_exome_gcs.sh "$GCS_CLUSTER"-m:prep_data_exome_gcs.sh \
  --zone us-central1-a
gcloud compute ssh "$GCS_CLUSTER"-m \
  --command "./prep_data_exome_gcs.sh /user/$USER/exome_spark_eval" \
  --zone us-central1-a

# Run pipeline
source exome_pipeline_single_gcs_hdfs.sh

gcloud compute ssh "$GCS_CLUSTER"-m \
  --command "hadoop fs -ls hdfs://${GCS_CLUSTER}-m:8020/user/$USER/exome_spark_eval/out/NA12878.ga2.exome.maq.raw.vcf" \
  --zone us-central1-a

# (Optional) Retrieve output
#gcloud compute ssh "$GCS_CLUSTER"-m \
#  --command "hadoop fs -get hdfs://${GCS_CLUSTER}-m:8020/user/$USER/exome_spark_eval/out/NA12878.ga2.exome.maq.raw.vcf NA12878.ga2.exome.maq.raw.spark.vcf" \
#  --zone us-central1-a
#gcloud compute scp "$GCS_CLUSTER"-m:NA12878.ga2.exome.maq.raw.spark.vcf NA12878.ga2.exome.maq.raw.spark.vcf \
#  --zone us-central1-a

# Delete cluster
gcloud dataproc clusters delete "$GCS_CLUSTER"