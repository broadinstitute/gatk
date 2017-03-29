#!/bin/bash

set -eu

# This script initializes the master and worker nodes on a Google Dataproc
# Spark cluster to prepare them to run the GATK-SV pipeline.
#
# On the master node we launch a hadoop distcp command to copy 
#    the reference, and
#    data sets
# down from a bucket to the hdfs file system on the cluster.
#
# On the worker nodes we copy the bwa index image file 
# from a bucket to each node's local disk.

NAME=$(/usr/share/google/get_metadata_value hostname)
ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
REFLOC=$(/usr/share/google/get_metadata_value attributes/reference)
SAMPLE=$(/usr/share/google/get_metadata_value attributes/sample)

if [[ "${ROLE}" == 'Master' ]]; then
    hadoop distcp "$REFLOC"/* hdfs://"$NAME":8020/reference/ &
    hadoop distcp "$SAMPLE"/* hdfs://"$NAME":8020/data/
else
    if [[ "$REFLOC" =~ .+/$ ]]; then
        mkdir -p /reference && gsutil -m cp "$REFLOC"*.img /reference/
    else
        mkdir -p /reference && gsutil -m cp "$REFLOC"/*.img /reference/
    fi
fi
