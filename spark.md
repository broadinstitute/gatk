GATK4 on Spark
================

[Apache Spark](https://spark.apache.org/) is a fast and general engine for large-scale data processing. GATK4 can run on any Spark cluster, such as an
on-premise Hadoop cluster with HDFS storage and the Spark runtime. You can also use the local runner.


Running GATK4 Commands Locally on Spark (without a Cluster)
------------

This allows you to run GATK4 command lines using Spark running in the local JVM. It is useful for testing
commands against local files using the Spark runner.

First, [download](https://spark.apache.org/downloads.html) and unpack the Spark tarball on your local machine.
Version 1.5.x for Hadoop 2.6 is recommended.

For convenience, add Spark to your `PATH`:

```bash
export SPARK_HOME=/path/to/spark/installation
export PATH=$PATH:$SPARK_HOME/bin
```

The GATK build produces a JAR file that is suitable for running Spark jobs with (note that it differs from the
"fat JAR" that is used for running regular GATK commands). Build the Spark-compatible jar with:

```bash
gradle clean sparkJar
```

You should see a JAR file in the _build/libs_ directory with a `-spark` extension.

Now you can run a pipeline using `spark-submit` in place of the regular `java -jar` invocation. For example, to run
the count reads program to count the number of reads in a BAM file, type the following:

```bash
spark-submit \
  build/libs/gatk-all-*-spark.jar CountReadsSpark \
    --input ./src/test/resources/org/broadinstitute/hellbender/tools/flag_stat.bam \
    --output countreads
```

You can inspect the output with

```bash
cat countreads*
```

Running GATK4 Commands on a Spark Cluster
------------

In this scenario, your input and output files reside on HDFS, and Spark will run in a distributed fashion on the cluster.
The Spark documentation has a good [overview of the architecture](https://spark.apache.org/docs/latest/cluster-overview.html).

Note that if you don't have a dedicated cluster you can run Spark in
[standalone mode](https://spark.apache.org/docs/latest/spark-standalone.html) on a single machine, which exercises
the distributed code paths, albeit on a single node.

The first step is to copy a BAM file to HDFS, which can be achieved using the `hadoop` command as follows.

```bash
hadoop fs -put ./src/test/resources/org/broadinstitute/hellbender/tools/flag_stat.bam flag_stat.bam
```

Next, run `spark-submit` again, but with a few extra configuration properties set, which are suitable for a cluster
environment:

```bash
NAMENODE=<hostname of the HDFS namenode>
SPARK_MASTER=yarn-client
spark-submit \
  --master $SPARK_MASTER \
  --conf spark.driver.userClassPathFirst=true \
  --conf spark.executor.userClassPathFirst=true \
  --conf spark.io.compression.codec=lzf \
  build/libs/gatk-all-*-spark.jar CountReadsSpark \
    --input hdfs://$NAMENODE/user/$USER/flag_stat.bam \
    --output hdfs://$NAMENODE/user/$USER/out/countreads \
    --sparkMaster $SPARK_MASTER
```

This requires a bit more explanation. The environment variables `NAMENODE` and `SPARK_MASTER` are used to tell the
Spark client where to find the cluster. `NAMENODE` is used to construct the HDFS paths for the `--input` and `--output`
arguments, and should point to the HDFS namenode in the cluster you are using.

`SPARK_MASTER` is used to tell the Spark client where to submit the job to. In this case we are running Spark on YARN,
so `yarn-client` is appropriate. You can find more about Spark master URLs
[here](https://spark.apache.org/docs/latest/submitting-applications.html#master-urls).

The `spark-submit` command takes arguments - and these go before the argument specifying the JAR file. In this case,
we've set some configuration variables: `spark.driver.userClassPathFirst` and `spark.executor.userClassPathFirst` are
set which means that our JAR file is put on the classpath ahead of the Spark and Hadoop JAR files and dependencies. We
need to do this so that we don’t get conflicts with older library versions that Hadoop uses.

Also, the property `spark.io.compression.codec` is changed from its default `snappy`, since that can cause classloading
issues on some clusters. (Note that this setting is only for over-the-wire compression, not for file compression.)

When the command completes successfully you will see files in the `out` directory in HDFS:

```bash
hadoop fs -ls -R out/countreads*
```

You can inspect the output with

```bash
hadoop fs -cat out/countreads*
```

### Adjusting the Number of Executors and Memory Settings

For pipelines running on large datasets (such as typical BAM files), you will need to adjust the number of Spark
executors and the amount of memory that each executor uses. For example, the following runs a pipeline on 100 executors
each with one core and 3GB of memory:

```bash
NAMENODE=<hostname of the HDFS namenode>
SPARK_MASTER=yarn-client
spark-submit \
  --master $SPARK_MASTER \
  --driver-memory 3G \
  --num-executors 100 \
  --executor-cores 1 \
  --executor-memory 3G \
  --conf spark.driver.userClassPathFirst=true \
  --conf spark.executor.userClassPathFirst=true \
  --conf spark.io.compression.codec=lzf \
  --conf spark.shuffle.io.preferDirectBufs=false \
  build/libs/gatk-all-*-spark.jar CountReadsSpark \
    --input hdfs://$NAMENODE/user/$USER/bam/NA12877_S1.bam \
    --output hdfs://$NAMENODE/user/$USER/out/countreads \
    --sparkMaster $SPARK_MASTER
```

The input file, _NA12877_S1.bam_, is 121GB in size, which gets split into ~950 pieces of input (by the Hadoop input
format), so roughly ten rounds of 100 concurrent tasks are needed to process the file. On a larger cluster it would be
advisable to increase the number of executors (or executor cores, or both) so that the pipeline completed faster.

While your Spark job is running, the [Spark UI](http://spark.apache.org/docs/latest/monitoring.html) is an excellent place to monitor the  progress. Additionally, if you're running tests, then by adding `-Dgatk.spark.debug=true` you can run a single Spark test and look at the Spark UI (on [http://localhost:4040/](http://localhost:4040/)) as it runs.

You can find more information about tuning Spark and choosing good values for number of executors and memory settings
at the following:

* [Tuning Spark](https://spark.apache.org/docs/latest/tuning.html)
* [How-to: Tune Your Apache Spark Jobs (Part 1)](http://blog.cloudera.com/blog/2015/03/how-to-tune-your-apache-spark-jobs-part-1/)
* [How-to: Tune Your Apache Spark Jobs (Part 2)](http://blog.cloudera.com/blog/2015/03/how-to-tune-your-apache-spark-jobs-part-2/)
 

