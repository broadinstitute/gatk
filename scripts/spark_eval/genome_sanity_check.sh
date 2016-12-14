#!/usr/bin/env bash

. utils.sh

time_gatk "CountReadsSpark -I gs://hellbender/q4_spark_eval/WGS-G94982-NA12878.bam" 4 4 4g 4g
