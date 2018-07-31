#!/bin/bash

set -eu

# This script initializes the master and worker nodes on a Google Dataproc
# Spark cluster to prepare them to run the GATK-SV pipeline.
#
# On the worker nodes we copy the bwa index image file
# from a bucket to each node's local disk.

REFLOC=$(/usr/share/google/get_metadata_value attributes/reference)
if  [[ ! $REFLOC =~ .+/$ ]]; then
    REFLOC+="/"
fi

ROLE=$(/usr/share/google/get_metadata_value attributes/dataproc-role)
if [[ ! "${ROLE}" == 'Master' ]]; then
    # the /mnt/1/ prefix is the default mounting location of the ssd if the cluster is created with a local ssd for the worker nodes
    if [[ -d "/mnt/1/" ]]; then
        mkdir -p /mnt/1/reference && gsutil -m cp "$REFLOC"*.img /mnt/1/reference/
        ln -s /mnt/1/reference /reference
    else
        mkdir -p /reference && gsutil -m cp "$REFLOC"*.img /reference/
    fi
fi
