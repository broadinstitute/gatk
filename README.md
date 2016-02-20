[![Build Status](https://travis-ci.org/broadinstitute/gatk.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gatk/badge.svg?branch=master)](https://coveralls.io/r/broadinstitute/gatk?branch=master)

**This project is in an early stage of development.  It is subject to change without warning. Do not use this code for production work.**  

###GATK 4 (codename Hellbender)

This repository contains the next generation GATK/Picard engine and the freely available tools (see [LICENSE](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT)). See also the [gatk-protected](https://github.com/broadinstitute/gatk-protected) repository for additional GATK tools that are distributed under a different license.

GATK4 aims to bring together well-established tools from the [GATK](http://www.broadinstitute.org/gatk) and
[Picard](http://broadinstitute.github.io/picard/) codebases under a single simplified, streamlined framework,
and to enable selected tools to be run in a massively parallel way on local clusters or in the cloud using
[Apache Spark](http://spark.apache.org/).

This project is in an alpha development stage and is not yet ready for general use.

If you are looking for the current version of GATK to use in production work, please see the [GATK website](http://www.broadinstitute.org/gatk), where you can download a precompiled executable, read documentation, ask questions and receive technical support.

If you are looking for the codebase of the current production version of GATK, please see either the [GATK development framework repository](https://github.com/broadgsa/gatk/) or the [full GATK tools repository](https://github.com/broadgsa/gatk-protected/).

##Requirements
* Java 8
* Git 2.5 or greater
* Optional, but recommended:
    * Gradle 2.11 or greater, needed for building the GATK. We recommend using the `./gradlew` script which will
      download and use an appropriate gradle version automatically), see examples below.
    * Python 2.6 or greater (needed for running the `gatk-launch` frontend script)
    * R 3.1.3 (needed for producing plots in certain tools, and for running the test suite)
    * [git-lfs](https://git-lfs.github.com/) 1.1.0 or greater (needed to download large files for the complete test suite).
      Run `git lfs install` after downloading, followed by `git lfs pull` to download the large files. The download is ~500 MB.

##Quick Start Guide

* Build the GATK: `./gradlew installAll`
* Get help on running the GATK: `./gatk-launch --help`
* Get a list of available tools: `./gatk-launch --list`
* Run a tool: `./gatk-launch PrintReads -I src/test/resources/NA12878.chr17_69k_70k.dictFix.bam -O output.bam`
* Get help on a particular tool: `./gatk-launch PrintReads --help`

##Building GATK4

* **Recommended method:** To build everything you need to run all GATK tools, run **`./gradlew installAll`**.
    * This does both the standard build and the Spark build described below.

* To _only_ do a standard build to save time, run **`./gradlew installDist`**.
    * The resulting executable will be `build/install/gatk/bin/gatk`,
      and will be suitable for running non-Spark tools and for running Spark tools locally (not on a cluster).
    * This executable is **not** a fully-packaged executable that can be copied across machines, it's just
      for running in-place. To create a packaged version of this build, run `./gradlew shadowJar`. The
      packaged jar will be in `build/libs/` with a name like `gatk-all-*-SNAPSHOT-shadowJar.jar`.

* To _only_ build a special jar suitable for running spark tools on a cluster, run **`./gradlew installSpark`**.
    * The resulting jar will be in `build/libs/` and will have a name like `gatk-all-*-SNAPSHOT-spark.jar`
    * This jar will not include Spark and Hadoop libraries, in order to allow
      the versions of Spark and Hadoop installed on your cluster to be used.

* To remove previous builds, run **`./gradlew clean`**

* For faster gradle operations, add `org.gradle.daemon=true` to your `~/.gradle/gradle.properties` file.
  This will keep a gradle daemon running in the background and avoid the ~6s gradle start up time on every command.

* Gradle keeps a cache of dependencies used to build GATK.  By default this goes in `~/.gradle`.  If there is insufficient free space in your home directory, you can change the location of the cache by setting the `GRADLE_USER_HOME` environment variable.


##Running GATK4

* The standard way to run GATK4 tools is via the **`gatk-launch`** wrapper script located in the root directory of a clone of this repository.
    * Requires Python 2.6 or greater.
    * You need to have built the GATK as described in the "Building GATK4" section above before running this script.
    * Can run non-Spark tools as well as Spark tools, and can run Spark tools locally, on a Spark cluster, or on Google Cloud Dataproc.

* For help on using `gatk-launch` itself, run **`./gatk-launch --help`**

* To print a list of available tools, run **`./gatk-launch --list`**.
    * Spark-based tools will have a name ending in `Spark` (eg., `BaseRecalibratorSpark`) and will be in one of the
      `Spark` categories. All other tools are non-Spark-based.

* To print help for a particular tool, run **`./gatk-launch ToolName --help`**.

* To run a non-Spark tool, or to run a Spark tool locally, the syntax is:
 **`./gatk-launch ToolName toolArguments`**.
    * Examples:

        ```
        ./gatk-launch PrintReads -I input.bam -O output.bam
        ```

        ```
        ./gatk-launch PrintReadsSpark -I input.bam -O output.bam
        ```

* To run a Spark tool on a Spark cluster, the syntax is:

  **`./gatk-launch ToolName toolArguments -- --sparkRunner SPARK --sparkMaster <master_url> additionalSparkArguments`**
    * Examples:

        ```
        ./gatk-launch PrintReadsSpark -I hdfs://path/to/input.bam -O hdfs://path/to/output.bam \
            -- \
            --sparkRunner SPARK --sparkMaster <master_url>
        ```

        ```
        ./gatk-launch PrintReadsSpark -I hdfs://path/to/input.bam -O hdfs://path/to/output.bam \
            -- \
            --sparkRunner SPARK --sparkMaster <master_url> \
            --num-executors 5 --executor-cores 2 --executor-memory 4g \
            --conf spark.yarn.executor.memoryOverhead=600
        ```

    * Note that the Spark-specific arguments are separated from the tool-specific arguments by a `--`.
    * Running a Spark tool on a cluster requires Spark to have been installed from http://spark.apache.org/, since
      `gatk-launch` invokes the `spark-submit` tool behind-the-scenes.

* To run a Spark tool on Google Cloud Dataproc:
    * You must have a [Google cloud services](https://cloud.google.com/) account, and have spun up a Dataproc cluster
      in the [Google Developer's console](https://console.developers.google.com).
    * You need to have installed the Google Cloud SDK from https://cloud.google.com/sdk/, since
      `gatk-launch` invokes the `gcloud` tool behind-the-scenes. As part of the installation, be sure
      that you follow the `gcloud` setup instructions [here](https://cloud.google.com/sdk/gcloud/).
    * Your inputs to the GATK need to be in Google Cloud Storage buckets, and should be specified on
      your GATK command line using the syntax `gs://my-gcs-bucket/path/to/my-file`

    Once you're set up, you can run a Spark tool on your Dataproc cluster using a command of the form:

    **`./gatk-launch ToolName toolArguments -- --sparkRunner GCS --cluster myGCSCluster additionalSparkArguments`**

    * Examples:

        ```
        ./gatk-launch PrintReadsSpark \
            -I gs://my-gcs-bucket/path/to/input.bam \
            -O gs://my-gcs-bucket/path/to/output.bam \
            -- \
            --sparkRunner GCS --cluster myGCSCluster
        ```

        ```
        ./gatk-launch PrintReadsSpark \
            -I gs://my-gcs-bucket/path/to/input.bam \
            -O gs://my-gcs-bucket/path/to/output.bam \
            -- \
            --sparkRunner GCS --cluster myGCSCluster \
            --num-executors 5 --executor-cores 2 --executor-memory 4g \
            --conf spark.yarn.executor.memoryOverhead=600
        ```

    * Note that the spark-specific arguments are separated from the tool-specific arguments by a `--`.
    * If you want to avoid uploading the GATK jar to GCS on every run, set the `GATK_GCS_STAGING`
      environment variable to a bucket you have write access to (eg., `export GATK_GCS_STAGING=gs://<my_bucket>/`)
    * For further instructions on using GATK with Google Cloud Dataproc, see our
      [Cloud Dataproc wiki page](https://github.com/broadinstitute/gatk/wiki/Running-GATK-on-Cloud-Dataproc)

* If you don't want to run the GATK via the `gatk-launch` script, it's possible to run non-Spark and local Spark
  tools directly using the `build/install/gatk/bin/gatk` executable after a `gradle installDist`, and to run Spark tools
  on a cluster or the cloud by building a Spark jar with `gradle installSpark` and passing the resulting jar in `build/libs/`
  directly to either `spark-submit` or `gcloud`.

##Testing GATK4

* To run all tests, run **`gradle test`**.
    * Test report is in `build/reports/tests/index.html`.
    * What will happen depends on the value of the `CLOUD` environment variable: if it's `false` or
      unset then only local tests are run, if it's `mandatory` then it'll run only the cloud tests,
      and if it's `together` it will run both sets of tests.
    * Note that `git lfs` must be installed and set up as described in the "Requirements" section above
      in order for all tests to pass.
    * Cloud tests require being logged into `gcloud` and authenticated with a project that has access
      to the test data.

* To run just a single test class, run **`gradle test -Dtest.single=MyTestClass`**
    * example: `gradle test -Dtest.single=ReadUtilsUnitTest`

* To run tests and compute coverage reports, run **`gradle jacocoTestReport`**. The report is then in `build/reports/jacoco/test/html/index.html`.
  (IntelliJ 14 has a good coverage tool that is preferable for development).

* We use [Travis-CI](https://travis-ci.org/broadinstitute/gatk) as our continuous integration provider.

    * Before merging any branch make sure that all required tests pass on travis.
    * Every travis build will upload the test results to our gatk google bucket.
      A link to the uploaded report will appear at the very bottom of the travis log.
      Look for the line that says `See the test report at`.
      If TestNG itself crashes there will be no report generated.

##General guidelines for GATK4 developers

* **Do not put private or restricted data into the repo.**

* **Try to keep datafiles under 100kb in size.** Larger test files should go into `src/test/resources/large`, and must be
  managed using `git lfs` by running `git lfs track <file>` on each new large file before commit.

* GATK4 is BSD licensed.  The license is in the top level LICENSE.TXT file.  Do not add any additional license text or accept files with a license included in them.

* Each tool should have at least one good end-to-end integration test with a check for expected output, plus high-quality unit tests for all non-trivial utility methods/classes used by the tool. Although we have no specific coverage target, coverage should be extensive enough that if tests pass, the tool is guaranteed to be in a usable state.

* All newly written code must have good test coverage (>90%).

* All bug fixes must be accompanied by a regression test.

* All pull requests must be reviewed before merging to master (even documentation changes).

* Don't issue or accept pull requests that introduce warnings. Warnings must be addressed or suppressed.

* Don't issue or accept pull requests that significantly decrease coverage (less than 1% decrease is sort of tolerable). 

* Don't override `clone()` unless you really know what you're doing. If you do override it, document thoroughly. Otherwise, prefer other means of making copies of objects.

* Don't use `toString()` for anything other than human consumption (ie. don't base the logic of your code on results of `toString()`.)

* For logging, use [org.apache.logging.log4j.Logger](https://logging.apache.org/log4j/2.0/log4j-api/apidocs/org/apache/logging/log4j/Logger.html)

* We mostly follow the [Google Java Style guide](http://google-styleguide.googlecode.com/svn/trunk/javaguide.html)

* Git: Don't push directly to master - make a pull request instead. 

* Git: Rebase and squash commits when merging.

* If you push to master or mess the commit history, you owe us 1 growler or tasty snacks at happy hour. If you break the master build, you owe 3 growlers (or lots of tasty snacks). Beer may be replaced by wine (in the color and vintage of buyer's choosing) in proportions of 1 growler = 1 bottle. 

##Note on 2bit Reference
* Note: Some GATK Spark tools by default require the reference file to be in 2bit format (notably `BaseRecalibratorSpark`,`BQSRPipelineSpark` and `ReadsPipelineSpark`). You can convert your fasta to 2bit by using the `faToTwoBit` utility from [UCSC](http://hgdownload.soe.ucsc.edu/admin/exe/) - see also the [documentation for `faToTwoBit`](https://genome.ucsc.edu/goldenpath/help/blatSpec.html#faToTwoBitUsage).

##R Dependency
Certain GATK tools may optionally generate plots if R is installed.  We recommend **R v3.1.3** if you want to produce plots.  If you are uninterested in plotting, R is still required by several of the unit tests.  Plotting is currently untested and should be viewed as a convinience rather than a primary output.  

R installation is not part of the gradle build.  See http://cran.r-project.org/ for general information on installing R for your system.
* for ubuntu see these [ubuntu specific instructions](http://cran.r-project.org/bin/linux/ubuntu/README)
* for OSX we recommend installation through [homebrew](http://brew.sh/)
```
brew tap homebrew/science
brew install R
```

The plotting R scripts require certain R packages to be installed. You can install these by running `scripts/install_R_packages.R`.  Either run it as superuser to force installation into the sites library or run interactively and create a local library.
```
sudo Rscript scripts/install_R_packages.R
```
**or**
```
R 
source("scripts/install_R_packages.R")
```

##Creating a GATK project in the IntelliJ IDE:

* Ensure that you have `gradle` and the Java 8 JDK installed

* Install the TestNG plugin (in preferences)

* Clone the GATK repository using git

* In IntelliJ, go to File -> "Import Project"

* Select the root directory of your GATK clone, then "Ok"

* Select "Import project from external model", then "Gradle", then "Next"

* Ensure that "Gradle project" points to the build.gradle file in the root of your GATK clone

* Select "Use auto-import" and "Use default gradle wrapper".

* Click "Finish"

* After downloading project dependencies, IntelliJ should open a new window with your GATK project

* In File -> "Project Structure" -> "Project", set the "Project SDK" to your Java 1.8 JDK, and "Project language level" to 8 (you may need to add your Java 8 JDK under "Platform Settings" -> SDKs if it isn't there already). Then click "Apply"/"Ok".


##Setting up debugging in IntelliJ

* Follow the instructions above for creating an IntelliJ project for GATK

* Go to Run -> "Edit Configurations", then click "+" and add a new "Application" configuration

* Set the name of the new configuration to something like "GATK debug"

* For "Main class", enter `org.broadinstitute.hellbender.Main`

* Ensure that "Use classpath of module:" is set to use the "gatk" module's classpath

* Enter the arguments for the command you want to debug in "Program Arguments"

* Click "Apply"/"Ok"

* Set breakpoints, etc., as desired, then select "Run" -> "Debug" -> "GATK debug" to start your debugging session

* In future debugging sessions, you can simply adjust the "Program Arguments" in the "GATK debug" configuration as needed

##Updating the Intellij project when dependencies change
If there are dependency changes in `build.gradle` it is necessary to refresh the gradle project. This is easily done with the following steps.

* Open the gradle tool window  ( "View" -> "Tool Windows" -> "Gradle" )
* Click the refresh button in the Gradle tool window.  It is in the top left of the gradle view and is represented by two blue arrows.

##Uploading Archives to Sonatype (to make them available via maven central)
To upload snapshots to Sonatype you'll need the following:

* You must have a registered account on the sonatype JIRA (and be approved as a gatk uploader)
* You need to configure several additional properties in your `/~.gradle/gradle.properties` file

* If you want to upload a release instead of a snapshot you will additionally need to have access to the gatk signing key and password

```
#needed for snapshot upload
sonatypeUsername=<your sonatype username>
sonatypePassword=<your sonatype password>

#needed for signing a release
signing.keyId=<gatk key id>
signing.password=<gatk key password>
signing.secretKeyRingFile=/Users/<username>/.gnupg/secring.gpg
```

To perform an upload, use
```
./gradlew uploadArchives
```

Currently all builds are considered snapshots.  The archive name is based off of `git describe`.

###Further Reading on Spark

[Apache Spark](https://spark.apache.org/) is a fast and general engine for large-scale data processing.
GATK4 can run on any Spark cluster, such as an on-premise Hadoop cluster with HDFS storage and the Spark
runtime, as well as on the cloud using Google Dataproc.

In a cluster scenario, your input and output files reside on HDFS, and Spark will run in a distributed fashion on the cluster.
The Spark documentation has a good [overview of the architecture](https://spark.apache.org/docs/latest/cluster-overview.html).

Note that if you don't have a dedicated cluster you can run Spark in
[standalone mode](https://spark.apache.org/docs/latest/spark-standalone.html) on a single machine, which exercises
the distributed code paths, albeit on a single node.

While your Spark job is running, the [Spark UI](http://spark.apache.org/docs/latest/monitoring.html) is an excellent place to monitor the  progress.
Additionally, if you're running tests, then by adding `-Dgatk.spark.debug=true` you can run a single Spark test and
look at the Spark UI (on [http://localhost:4040/](http://localhost:4040/)) as it runs.

You can find more information about tuning Spark and choosing good values for important settings such as the number
of executors and memory settings at the following:

* [Tuning Spark](https://spark.apache.org/docs/latest/tuning.html)
* [How-to: Tune Your Apache Spark Jobs (Part 1)](http://blog.cloudera.com/blog/2015/03/how-to-tune-your-apache-spark-jobs-part-1/)
* [How-to: Tune Your Apache Spark Jobs (Part 2)](http://blog.cloudera.com/blog/2015/03/how-to-tune-your-apache-spark-jobs-part-2/)


###How to contribute
(Note: section inspired by, and some text copied from, [Apache Parquet](https://github.com/apache/parquet-mr))
 
We welcome all contributions to the GATK project. The contribution can be a [issue report]( https://github.com/broadinstitute/gatk/issues) 
or a [pull request](https://github.com/broadinstitute/gatk/pulls). If you're not a committer, you will 
need to [make a fork](https://help.github.com/articles/fork-a-repo/) of the gatk repository 
and [issue a pull request](https://help.github.com/articles/be-social/) from your fork.

To become a committer, you need to make several high-quality code contributions and be approved by the current committers.

For ideas on what to contribute, check issues labeled ["Help wanted (Community)"](https://github.com/broadinstitute/gatk/issues?q=is%3Aopen+is%3Aissue+label%3A%22Help+Wanted+%28Community%29%22). Comment on the issue to indicate you're interested in contibuting code and for sharing your questions and ideas.

To contribute a patch:
* Break your work into small, single-purpose patches if possible. It’s much harder to merge in a large change with a lot of disjoint features.
* Submit the patch as a GitHub pull request against the master branch. For a tutorial, see the GitHub guides on [forking a repo](https://help.github.com/articles/fork-a-repo/) and [sending a pull request](https://help.github.com/articles/be-social/). If applicable, include the issue number in the pull request name.
* Make sure that your code passes all our tests. You can run the tests with `./gradlew test` in the root directory.
* Add tests for all new code you've written. We prefer unit tests but high quality integration tests that use small amounts of data are acceptable.
* Follow the [**General guidelines for GATK4 developers**](https://github.com/broadinstitute/gatk#general-guidelines-for-gatk4-developers).

We tend to do fairly close readings of pull requests, and you may get a lot of comments. Some things to consider:
* Write tests for all new code.
* Document all classes and public methods.
* For all public methods, check validity of the arguments and throw `IllegalArgumentException` if invalid.
* Use braces for control constructs, `if`, `for` etc.
* Make classes, variables, parameters etc `final` unless there is a strong reason not to.
* Give your operators some room. Not `a+b` but `a + b` and not `foo(int a,int b)` but `foo(int a, int b)`.
* Generally speaking, stick to the [Google Java Style guide](http://google-styleguide.googlecode.com/svn/trunk/javaguide.html)

Thank you for getting involved!

##Discussions
* [GATK forum](http://gatkforums.broadinstitute.org/) for general discussions on how to use the GATK.
* [Issue tracker](https://github.com/broadinstitute/gatk/issues) to report errors and enhancement ideas. 
* Discussions also take place in github pull requests
* For committers, we have a publicly-visible google group [gatk-dev-public](https://groups.google.com/a/broadinstitute.org/forum/?hl=en#!forum/gatk-dev-public)
* For committers, we have a hipchat room at the Broad called 'Hellbender (aka GATK4)'.

##Authors
The authors list is maintained in the [AUTHORS](https://github.com/broadinstitute/gatk/edit/master/AUTHORS) file. 
See also the [Contributors](https://github.com/broadinstitute/gatk/graphs/contributors) list at github. 

##License
Licensed under the BSD License. See the [LICENSE.txt](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT) file.
