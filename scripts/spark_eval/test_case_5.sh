#!/usr/bin/env bash

. utils.sh

for exec_mem in 1g 500m 256m 128m 64m; do
  time_gatk "CountReadsSpark -I hdfs:///user/tom/q4_spark_eval/WGS-G94982-NA12878.bam" 4 4 $exec_mem 1g
done
