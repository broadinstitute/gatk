[![Build Status](https://travis-ci.org/broadinstitute/gatk-protected.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk-protected)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gatk-protected/badge.svg?branch=master&service=github)](https://coveralls.io/github/broadinstitute/gatk-protected?branch=master)
[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

GATK4-Protected (codename Hellbender-protected)
===============================================

GATK4 development of the previously license-protected part of the toolkit. The contents of this repo will be merged into broadinstitute/gatk in the near future.

This README is aimed at developers.  For user information, please see the [GATK 4 forum](http://gatkforums.broadinstitute.org/gatk/categories/gatk-4-alpha)

Read GATK 4 README
------------------------

**Please refer to the [GATK 4 public repo README](https://github.com/broadinstitute/hellbender/blob/master/README.md) for general guidelines and how to setup your development environment.**

Requirements
------------
* R 3.1.3 see additional requirements below: [R package requirements](#r-required-packages)

* Java 8

* (Developers) Gradle 2.13 is needed for building the GATK. We recommend using the `./gradlew` script which will
download and use an appropriate gradle version automatically.

* (Developers) git lfs 1.1.0 (or greater) is needed for testing GATK-Protected builds.  It is needed to download large files for the complete test suite. Run ``git lfs install`` after downloading, followed by ``git lfs pull`` to download the large files. The download is ~500 MB.

#### R Required Packages
R packages can be installed using the install_R_packages.R script inside the scripts directory.

## Building GATK4

* To do a fast build that lets you run GATK tools from within a git clone locally (but not on a cluster), run:
        
        ./gradlew installDist
        
* To do a slower build that lets you run GATK tools from within a git clone both locally and on a cluster, run:

        ./gradlew installAll
     
* To build a fully-packaged GATK jar that can be distributed and includes all dependencies needed for running tools locally, run:

        ./gradlew localJar
        
    * The resulting jar will be in `build/libs` with a name like `gatk-protected-package-VERSION-local.jar`
    
* To build a fully-packaged GATK jar that can be distributed and includes all dependencies needed for running spark tools on a cluster, run:

        ./gradlew sparkJar
        
    * The resulting jar will be in `build/libs` with a name like `gatk-protected-package-VERSION-spark.jar`
    * This jar will not include Spark and Hadoop libraries, in order to allow the versions of Spark and Hadoop installed on your cluster to be used.

* To remove previous builds, run: 

        ./gradlew clean
        
## Running GATK4

* The standard way to run GATK4 tools is via the **`gatk-launch`** wrapper script located in the root directory of a clone of this repository.
    * Requires Python 2.6 or greater.
    * You need to have built the GATK as described in the "Building GATK4" section above before running this script.
    * There are three ways `gatk-launch` can be run:
        * from the root of your git clone after building
        * or, put the `gatk-launch` script within the same directory as fully-packaged GATK jars produced by `./gradlew localJar` and `./gradlew sparkJar`
        * or, the environment variables `GATK_LOCAL_JAR` and `GATK_SPARK_JAR` can be defined, and contain the paths to the fully-packaged GATK jars produced by `./gradlew localJar` and `./gradlew sparkJar` 
    * Can run non-Spark tools as well as Spark tools, and can run Spark tools locally, on a Spark cluster, or on Google Cloud Dataproc.

* For help on using `gatk-launch` itself, run **`./gatk-launch --help`**

* To print a list of available tools, run **`./gatk-launch --list`**.

* To print help for a particular tool, run **`./gatk-launch ToolName --help`**.

* To run a non-Spark tool, or to run a Spark tool locally, the syntax is:
    ```
    ./gatk-launch ToolName toolArguments
    ```
    *Examples:*
    ```
    ./gatk-launch PrintReads -I input.bam -O output.bam

    ./gatk-launch PrintReadsSpark -I input.bam -O output.bam
    ```

* To run a Spark tool on a Spark cluster, the syntax is:
    ```
    ./gatk-launch ToolName toolArguments -- --sparkRunner SPARK --sparkMaster <master_url> additionalSparkArguments
    ```
    *Example:*
    ```
    ./gatk-launch PrintReadsSpark -I hdfs://path/to/input.bam -O hdfs://path/to/output.bam \
        -- \
        --sparkRunner SPARK --sparkMaster <master_url> \
        --num-executors 5 --executor-cores 2 --executor-memory 4g \
        --conf spark.yarn.executor.memoryOverhead=600
    ```

* To run a Spark tool on Google Cloud Dataproc, the syntax is:
    ```
    ./gatk-launch ToolName toolArguments -- --sparkRunner GCS --cluster myGCSCluster additionalSparkArguments
    ```
    *Example:*
    ```
    ./gatk-launch PrintReadsSpark \
          -I gs://my-gcs-bucket/path/to/input.bam \
          -O gs://my-gcs-bucket/path/to/output.bam \
          -- \
          --sparkRunner GCS --cluster myGCSCluster \
          --num-executors 5 --executor-cores 2 --executor-memory 4g \
          --conf spark.yarn.executor.memoryOverhead=600
    ```
    
* **See the [GATK4 public README](https://github.com/broadinstitute/hellbender/blob/master/README.md) for full instructions on using `gatk-launch` to run tools on a Spark/Dataproc cluster.**

## Testing GATK4

* To run the tests, run **`./gradlew test`**.
    * Test report is in `build/reports/tests/test/index.html`. 
    * Note that `git lfs` must be installed and set up as described in the "Requirements" section above
      in order for all tests to pass.

* To run a subset of tests, use gradle's test filtering (see [gradle doc](https://docs.gradle.org/current/userguide/java_plugin.html)), e.g.,
    * `./gradlew test -Dtest.single=SomeSpecificTestClass`
    * `./gradlew test --tests *SomeSpecificTestClass`
    * `./gradlew test --tests all.in.specific.package*`
    * `./gradlew test --tests *SomeTest.someSpecificTestMethod`

* **See the [GATK4 public README](https://github.com/broadinstitute/hellbender/blob/master/README.md) for further information on running tests.**

The CNV case and PoN workflows (description and examples)
---------------------------------------------------------

This can be found [here](http://gatkforums.broadinstitute.org/gatk/discussion/6791/description-and-examples-of-the-steps-in-the-cnv-case-and-cnv-pon-creation-workflows)


Running the CNV case and PoN creation Workflows with premade Queue scripts
--------------------------------------------------------------------------

For [Broad Internal instructions](http://gatkforums.broadinstitute.org/gatk/discussion/6786/howto-run-gatk-cnv-using-premade-queue-scripts-broad-internal)

