#!/usr/bin/env bash

# Run the pipeline (Mark Duplicates, BQSR, Haplotype Caller) on small data using walkers.

. utils.sh

time_gatk_walker "MarkDuplicates -I ${GATK_HOME:-../..}/src/test/resources/large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam -O $HOME/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.md.bam --METRICS_FILE $HOME/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.metrics.txt --VALIDATION_STRINGENCY LENIENT"

