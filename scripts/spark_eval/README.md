# Spark Evaluation

This directory contains scripts for testing GATK pipelines on Spark - either on a dedicated cluster or on Google Cloud Dataproc.

## TL;DR

```bash
export GCS_CLUSTER=...

# Sanity check on small data (a few mins)
./run_gcs_cluster.sh small_reads-pipeline_gcs.sh

# Run on exome (<1hr)
nohup ./run_gcs_cluster.sh exome_reads-pipeline_gcs.sh &

# Run on genome (<2hrs)
NUM_WORKERS=20 nohup ./run_gcs_cluster.sh copy_genome_to_hdfs_on_gcs.sh genome_reads-pipeline_hdfs.sh &

# Check results
cat results/*
```

## Obtaining and preparing the data

There are three main datasets of increasing size: _small_ , _exome_, and _genome_ (WGS). The _small_ data is useful for sanity checking command lines before running them on the larger _exome_ and whole _genome_ datasets.

The datasets are stored in GCS buckets, so if you run using GCS input and output then there is no initial data preparation.

To copy the exome data into the cluster run:

```bash
./copy_exome_to_hdfs.sh
```

For the whole genome data, run:

```bash
./copy_genome_to_hdfs.sh
```

By default these copy the data into a directory in the current user's directory. To copy to different directory, add an argument like this:

```bash
./copy_genome_to_hdfs.sh /data/shared/spark_eval
```

Most of the data was obtained from the [GATK resource bundle](https://software.broadinstitute.org/gatk/download/bundle).

There is also some data in a GCS bucket for this evaluation: _gs://hellbender/q4_spark_eval/_.

## Running pipelines

The following shows how to run pipelines - from aligned reads to variants. The scripts follow a naming convention to make it easier to understand what they do:

```
<dataset>_<GATK tools>_<source/sink>.sh
```

So `small_reads-pipeline_gcs.sh` will run `ReadsPipelineSpark` on the `small` dataset in GCS (and writing output to GCS).

To run on Dataproc, make sure you set the `GCS_CLUSTER` environment variable:

```bash
export GCS_CLUSTER=...
```

### Running the exome pipeline on Dataproc (with data in HDFS)

For exome data, try `n1-standard-16` Dataproc worker instances, which have 60GB of memory, and 16 vCPUs. 2000GB of disk space per worker should be sufficient. Use 1 master + 5 workers. The master has lower resource requirements so `n1-standard-4`, 500GB disk is enough.

```bash
nohup ./exome_md-bqsr-hc_hdfs.sh &
```

This will take less than an hour.

### Running the whole genome pipeline on Dataproc (with data in HDFS)

For whole genome data, use the same instance types but try 10 workers.

```bash
nohup ./genome_md-bqsr-hc_hdfs.sh &
```

This will take a few hours.

### Running end-to-end

The following starts a GCS cluster, runs the given pipeline, then deletes the cluster.

```bash
nohup ./run_gcs_cluster.sh small_reads-pipeline_gcs.sh &
```

To copy the dataset to HDFS use a copy script first:

```bash
nohup ./run_gcs_cluster.sh copy_small_to_hdfs_on_gcs.sh small_reads-pipeline_hdfs.sh &
```

### More examples

```bash
# Exome Mark Duplicates, BQSR, Haplotype Caller on HDFS
nohup ./run_gcs_cluster.sh copy_exome_to_hdfs_on_gcs.sh exome_md-bqsr-hc_hdfs.sh &

# Exome ReadsSparkPipeline on HDFS
nohup ./run_gcs_cluster.sh copy_exome_to_hdfs_on_gcs.sh exome_reads-pipeline_hdfs.sh &

# Genome Mark Duplicates, BQSR, Haplotype Caller on HDFS using 20 workers
NUM_WORKERS=20 nohup ./run_gcs_cluster.sh copy_genome_to_hdfs_on_gcs.sh genome_md-bqsr-hc_hdfs.sh &

# Genome ReadsSparkPipeline on HDFS using 20 workers
NUM_WORKERS=20 nohup ./run_gcs_cluster.sh copy_genome_to_hdfs_on_gcs.sh genome_reads-pipeline_hdfs.sh &
```

## Running test cases

If you want to run tests from the [test definition document](https://docs.google.com/document/d/1OEfV2XNXdbGQQdW-gRaQYY2QRgGsWHUNKJUSAj3qFlE/edit), then run a command like the following:

```bash
nohup ./test_case_2.sh &
```

The output is saved to a CSV file (one per test case type), which can be analysed using _spark_eval.R_ to create plots.
