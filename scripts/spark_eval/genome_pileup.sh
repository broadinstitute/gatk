#!/usr/bin/env bash

. utils.sh

time_gatk "PileupSpark -I hdfs:///user/$USER/q4_spark_eval/WGS-G94982-NA12878.bam -O hdfs:///user/$USER/q4_spark_eval/out/pileup" 4 8 48g 4g
