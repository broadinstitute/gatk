# Spark Evaluation

This directory contains scripts for testing Spark command lines. It is based on this [test definition document](https://docs.google.com/document/d/1OEfV2XNXdbGQQdW-gRaQYY2QRgGsWHUNKJUSAj3qFlE/edit), but also has scripts for running full pipelines.

## Obtaining and preparing the data

Most of the data can be obtained from the [GATK resource bundle](https://software.broadinstitute.org/gatk/download/bundle).

There is also some data in a GCS bucket for this evaluation: _gs://hellbender/q4_spark_eval/_.

To copy the exome data into the cluster run:

```bash
./prep_data_exome_gcs.sh
```

For the whole genome data, run:

```bash
./prep_data_genome_gcs.sh
```

By default these copy the data into a directory in the current user's directory. To copy to different directory, add an argument like this:

```bash
./prep_data_genome_gcs.sh /data/shared/spark_eval
```

## Running test cases

If you want to run tests from the [test definition document](https://docs.google.com/document/d/1OEfV2XNXdbGQQdW-gRaQYY2QRgGsWHUNKJUSAj3qFlE/edit), then run a command like the following:

```bash
nohup ./test_case_2.sh &
```

The output is saved to a CSV file (one per test case type), which can be analysed using _spark_eval.R_ to create plots.

## Running pipelines

The following shows how to run pipelines - from aligned reads to variants.

### Running the exome pipeline on GCS (with data in HDFS)

For exome data, try `n1-standard-16` GCS worker instances, which have 60GB of memory, and 16 vCPUs. 2000GB of disk space per worker should be sufficient. Use 1 master + 5 workers. The master has lower resource requirements so `n1-standard-4`, 500GB disk is enough.

```bash
export API_KEY=...
export GCS_CLUSTER=...

nohup ./exome_pipeline_gcs_hdfs.sh &
```

This will take less than an hour.

### Running the whole genome pipeline on GCS (with data in HDFS)

For whole genome data, use the same instance types but try 10 workers.

```bash
export API_KEY=...
export GCS_CLUSTER=...

nohup ./genome_md_gcs_hdfs.sh &
nohup ./genome_bqsr_gcs_hdfs.sh &
nohup ./genome_hc_gcs_hdfs.sh &
```

This will take a few hours.
