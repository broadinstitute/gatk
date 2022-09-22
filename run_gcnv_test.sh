#!/bin/bash

GCNV_SCRIPT="/home/asmirnov/gatk/src/main/resources/org/broadinstitute/hellbender/tools/copynumber/cohort_denoising_calling.py"
READ_COUNTS_DIR="/home/asmirnov/gatk/src/test/resources/org/broadinstitute/hellbender/tools/copynumber/gcnv-sim-data/"

export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1

python $GCNV_SCRIPT  --read_count_tsv_files ${READ_COUNTS_DIR}/SAMPLE_000_counts.tsv  ${READ_COUNTS_DIR}/SAMPLE_001_counts.tsv  ${READ_COUNTS_DIR}/SAMPLE_002_counts.tsv \
	--modeling_interval_list ${READ_COUNTS_DIR}/sim_intervals_shard_0.annotated.tsv \
	--ploidy_calls_path ${READ_COUNTS_DIR}/contig-ploidy-calls/ \
        --output_calls_path ./calls \
	--output_model_path ./model \
        --output_tracking_path ./tracking

