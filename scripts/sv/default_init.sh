#!/bin/bash

set -eu

# This script initializes the master and worker nodes on a Google Dataproc
# Spark cluster to prepare them to run the GATK-SV pipeline.
#
# On the worker nodes we copy the bwa index image file
# from a bucket to each node's local disk.

ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
REFLOC=$(/usr/share/google/get_metadata_value attributes/reference)

if [[ ! "${ROLE}" == 'Master' ]]; then
    if  [[ ! $REFLOC =~ .+/$ ]]; then
      REFLOC+="/"
    fi
    mkdir -p /mnt/1/reference && gsutil -m cp "$REFLOC"*.img /mnt/1/reference/
fi
