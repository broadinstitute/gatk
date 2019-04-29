## Running the Spark WDL

This directory contains WDL for running GATK Spark tools on a
Dataproc cluster.

### Pre-requisites

From the top-level GATK directory, build a local Spark jar with

```bash
./gradlew sparkJar
```

Install Cromwell, and set an environment variable pointing to the JAR
(the path will need adjusting to where you installed it on your local
machine).

```bash
CROMWELL_JAR=~/sw/cromwell-39/cromwell-39.jar
```

### Configuration

Edit the JSON configuration in _spark_wdl/inputs_ to suit your own set
up. At a minimum you need to set

* `gcloud_service_account_key_file` to the GCP service account key file
* `gcloud_project` to the GCP project name you are using
* `output_vcf` to a GCS bucket that you have write access to

Strictly speaking, `gcloud_service_account_key_file` is optional and you
can authorize access to gcloud before running the WDL, with e.g.
`gcloud auth login`.

The WDL will automatically create a cluster and delete it at the end of
the job. A unique name for the cluster will be generated, unless the
`cluster_name` input is specified. If a cluster with that name already
exists then it will be used to run the job, otherwise a new cluster with
that name will be created.

The cluster will be deleted at the end of the job, unless
`delete_cluster` is set to `false`, in which case it will be left
running. Combined with setting the cluster name with `cluster_name`,
this allows clusters to be reused, which can be useful while debugging.

Clusters are created with a maximum age of three hours, so they will be
deleted after this time in any event.

The optional input `gatk_gcs_staging` can be set to a GCS path to stage
the GATK Spark JAR so it doesn't need to be uploaded for every job
(equivalent to the `GATK_GCS_STAGING` environment variable in GATK).

### Running workflows

In the _spark_wdl_ directory, run the workflow for the small dataset
with the following command:

```bash
java -jar $CROMWELL_JAR run ReadsPipelineSpark.wdl -i inputs/ReadsPipelineSpark_small.json
```

You can run on exome or genome data by using a different JSON file.
These run on clusters with workers, cores, and memory that are scaled
appropriately for the size of the input.

```bash
java -jar $CROMWELL_JAR run ReadsPipelineSpark.wdl -i inputs/ReadsPipelineSpark_exome.json
```

or

```bash
java -jar $CROMWELL_JAR run ReadsPipelineSpark.wdl -i inputs/ReadsPipelineSpark_genome.json
```

#### Using local GATK

WDL uses Docker to provide a controlled runtime environment. Sometimes
it's useful to be able to use a local version of GATK, during
development, for example. You can do this by using a local backend that
disables Docker.

You also need to change the WDL inputs in the JSON file for the `gatk`
binary, and the Spark JAR:

```
  "ReadsPipelineSparkWorkflow.gatk": "../../gatk",
  "ReadsPipelineSparkWorkflow.gatk_spark_jar": "../../build/libs/gatk-spark.jar",
```

Then run Cromwell using a special configuration file that disables
Docker for the local backend:

```bash
java -Dconfig.file=local-no-docker.conf -jar $CROMWELL_JAR run ReadsPipelineSpark.wdl -i inputs/ReadsPipelineSpark_small.json
```