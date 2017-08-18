[![Build Status](https://travis-ci.org/broadinstitute/gatk.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk)
[![codecov](https://codecov.io/gh/broadinstitute/gatk/branch/master/graph/badge.svg)](https://codecov.io/gh/broadinstitute/gatk)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/org.broadinstitute/gatk/badge.svg)](https://maven-badges.herokuapp.com/maven-central/org.broadinstitute/gatk)
[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)


***This project is in a pre-release stage of development.  It is subject to change without warning. Do not use this code for production work.***

***If you are looking for the current version of GATK to use in production work (ie., GATK3), please see the [GATK website](http://www.broadinstitute.org/gatk), where you can download a precompiled executable, read documentation, ask questions and receive technical support.***

### GATK 4

This repository contains the next generation of the Genome Analysis Toolkit (GATK). The contents
of this repository are 100% open source and released under the BSD 3-Clause license (see [LICENSE.TXT](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT)).

GATK4 aims to bring together well-established tools from the [GATK](http://www.broadinstitute.org/gatk) and
[Picard](http://broadinstitute.github.io/picard/) codebases under a streamlined framework,
and to enable selected tools to be run in a massively parallel way on local clusters or in the cloud using
[Apache Spark](http://spark.apache.org/). It also contains many newly developed tools not present in earlier
releases of the toolkit.

## Table of Contents
* [Requirements](#requirements)
* [Quick Start Guide](#quickstart)
* [Downloading GATK4](#downloading)
* [Building GATK4](#building)
* [Running GATK4](#running)
    * [Passing JVM options to gatk-launch](#jvmoptions)
    * [Running GATK4 with inputs on Google Cloud Storage](#gcs)
    * [Running GATK4 Spark tools on a Spark cluster](#sparkcluster)
    * [Running GATK4 Spark tools on Google Cloud Dataproc](#dataproc)
    * [Note on 2bit Reference](#2bit)
    * [Using R to generate plots](#R)
    * [Running the CNV workflows](#cnv_workflows)
* [For GATK Developers](#developers)
    * [General guidelines for GATK4 developers](#dev_guidelines)
    * [Testing GATK4](#testing)
    * [Using Git LFS to download and track large test data](#lfs)
    * [Creating a GATK project in the IntelliJ IDE](#intellij)
    * [Setting up debugging in IntelliJ](#debugging)
    * [Updating the Intellij project when dependencies change](#intellij_gradle_refresh)
    * [Setting up profiling using JProfiler](#jprofiler)
    * [Uploading Archives to Sonatype](#sonatype)
    * [Building GATK4 Docker images](#docker_building)
    * [Releasing GATK4](#releasing_gatk)
    * [Generating GATK4 documentation](#gatkdocs)
    * [Using Zenhub to track github issues](#zenhub)
* [Further Reading on Spark](#spark_further_reading)
* [How to contribute to GATK](#contribute)
* [Discussions](#discussions)
* [Authors](#authors)
* [License](#license)

## <a name="requirements">Requirements</a>
* Java 8
* Git 2.5 or greater
* Optional, but recommended:
    * Gradle 3.1 or greater, needed for building the GATK. We recommend using the `./gradlew` script which will
      download and use an appropriate gradle version automatically (see examples below).
    * Python 2.6 or greater (needed for running the `gatk-launch` frontend script)
    * R 3.1.3 (needed for producing plots in certain tools, and for running the test suite)
    * [git-lfs](https://git-lfs.github.com/) 1.1.0 or greater (needed to download large files for the complete test suite).
      Run `git lfs install` after downloading, followed by `git lfs pull` from the root of your git clone to download the large files. The download is several hundred megabytes.

## <a name="quickstart">Quick Start Guide</a>

* Build the GATK: `./gradlew bundle` (creates `gatk-VERSION.zip` in `build/`)
* Get help on running the GATK: `./gatk-launch --help`
* Get a list of available tools: `./gatk-launch --list`
* Run a tool: `./gatk-launch PrintReads -I src/test/resources/NA12878.chr17_69k_70k.dictFix.bam -O output.bam`
* Get help on a particular tool: `./gatk-launch PrintReads --help`

## <a name="downloading">Downloading GATK4</a>

You can download and run pre-built versions of GATK4 from the following places:

* Starting with the beta release, a zip archive with everything you need to run GATK4 can be downloaded for each release from the [github releases page](https://github.com/broadinstitute/gatk/releases).

* Starting with the beta release, you can download a GATK4 docker image from [our dockerhub repository](https://hub.docker.com/r/broadinstitute/gatk/). We also host unstable nightly development builds on [this dockerhub repository](https://hub.docker.com/r/broadinstitute/gatk-nightly/).
    * Within the docker image, run gatk-launch commands as usual from the default startup directory (/gatk).

## <a name="building">Building GATK4</a>

* **To do a full build of GATK4, run:**

        ./gradlew bundle
        
  Equivalently, you can just type:
  
        ./gradlew
        
    * This creates a zip archive in the `build/` directory with a name like `gatk-VERSION.zip` containing a complete standalone GATK distribution, including our launcher `gatk-launch`, both the local and spark jars, and this README.    
    * You can also run GATK commands directly from the root of your git clone after running this command.

* **Other ways to build:**
    * `./gradlew installDist`  
        * Does a *fast* build that only lets you run GATK tools from inside your git clone, and locally only (not on a cluster). Good for developers! 
    * `./gradlew installAll`
        * Does a *semi-fast* build that only lets you run GATK tools from inside your git clone, but works both locally and on a cluster. Good for developers!
    * `./gradlew localJar`
        * Builds *only* the GATK jar used for running tools locally (not on a Spark cluster). The resulting jar will be in `build/libs` with a name like `gatk-package-VERSION-local.jar`, and can be used outside of your git clone.
    * `./gradlew sparkJar`
        * Builds *only* the GATK jar used for running tools on a Spark cluster (rather than locally). The resulting jar will be in `build/libs` with a name like `gatk-package-VERSION-spark.jar`, and can be used outside of your git clone. 
        * This jar will not include Spark and Hadoop libraries, in order to allow the versions of Spark and Hadoop installed on your cluster to be used.

* **To remove previous builds, run:** 

        ./gradlew clean

* For faster gradle operations, add `org.gradle.daemon=true` to your `~/.gradle/gradle.properties` file.
  This will keep a gradle daemon running in the background and avoid the ~6s gradle start up time on every command.

* Gradle keeps a cache of dependencies used to build GATK.  By default this goes in `~/.gradle`.  If there is insufficient free space in your home directory, you can change the location of the cache by setting the `GRADLE_USER_HOME` environment variable.

## <a name="running">Running GATK4</a>

* The standard way to run GATK4 tools is via the **`gatk-launch`** wrapper script located in the root directory of a clone of this repository.
    * Requires Python 2.6 or greater (this includes Python 3.x)
    * You need to have built the GATK as described in the [Building GATK4](#building) section above before running this script.
    * There are several ways `gatk-launch` can be run:
        * Directly from the root of your git clone after building
        * By extracting the zip archive produced by `./gradlew bundle` to a directory, and running `gatk-launch` from there
        * Manually putting the `gatk-launch` script within the same directory as fully-packaged GATK jars produced by `./gradlew localJar` and/or `./gradlew sparkJar`
        * Defining the environment variables `GATK_LOCAL_JAR` and `GATK_SPARK_JAR`, and setting them to the paths to the GATK jars produced by `./gradlew localJar` and/or `./gradlew sparkJar` 
    * `gatk-launch` can run non-Spark tools as well as Spark tools, and can run Spark tools locally, on a Spark cluster, or on Google Cloud Dataproc.
    * ***Note:*** running with `java -jar` directly and bypassing `gatk-launch` causes several important system properties to not get set, including htsjdk compression level!
    
* For help on using `gatk-launch` itself, run **`./gatk-launch --help`**

* To print a list of available tools, run **`./gatk-launch --list`**.
    * Spark-based tools will have a name ending in `Spark` (eg., `BaseRecalibratorSpark`). Most other tools are non-Spark-based.

* To print help for a particular tool, run **`./gatk-launch ToolName --help`**.

* To run a non-Spark tool, or to run a Spark tool locally, the syntax is: **`./gatk-launch ToolName toolArguments`**.

* Examples:

  ```
  ./gatk-launch PrintReads -I input.bam -O output.bam
  ```

  ```
  ./gatk-launch PrintReadsSpark -I input.bam -O output.bam
  ```

#### <a name="jvmoptions">Passing JVM options to gatk-launch</a>

* To pass JVM arguments to GATK, run `gatk-launch` with the `--javaOptions` argument: 

    ```
    ./gatk-launch --javaOptions "-Xmx4G" <rest of command>
     
    ./gatk-launch --javaOptions "-Xmx4G -XX:+PrintGCDetails" <rest of command>
    ```

#### <a name="gcs">Running GATK4 with inputs on Google Cloud Storage:</a>

* Many GATK4 tools can read BAM or VCF inputs from a Google Cloud Storage bucket. Just use the "gs://" prefix:
  ```
  ./gatk-launch PrintReads -I gs://mybucket/path/to/my.bam -L 1:10000-20000 -O output.bam
  ```
* ***Important:*** You must set up your credentials first for this to work! There are three options:
    * Option (a): run in a Google Cloud Engine VM
        * If you are running in a Google VM then your credentials are already in the VM and will be picked up by GATK, you don't need to do anything special.
    * Option (b): use your own account
        * Install [Google Cloud SDK](https://cloud.google.com/sdk/)
        * Log into your account:
        ```
        gcloud auth application-default login
        ```
        * Done! GATK will use the application-default credentials you set up there.
    * Option (c): use a service account
        * Create a new service account on the Google Cloud web page and download the JSON key file
        * Install [Google Cloud SDK](https://cloud.google.com/sdk/)
        * Tell gcloud about the key file:
        ```
        gcloud auth activate-service-account --key-file "$PATH_TO_THE_KEY_FILE"
        ```
        * Set the `GOOGLE_APPLICATION_CREDENTIALS` environment variable to point to the file
        ```
        export GOOGLE_APPLICATION_CREDENTIALS="$PATH_TO_THE_KEY_FILE"
        ```
        * Done! GATK will pick up the service account. You can also do this in a VM if you'd like to override the default credentials.

#### <a name="sparkcluster">Running GATK4 Spark tools on a Spark cluster:</a>

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

* You can also omit the "--num-executors" argument to enable [dynamic allocation](https://spark.apache.org/docs/latest/job-scheduling.html#dynamic-resource-allocation) if you configure the cluster properly (see the Spark website for instructions).
* Note that the Spark-specific arguments are separated from the tool-specific arguments by a `--`.
* Running a Spark tool on a cluster requires Spark to have been installed from http://spark.apache.org/, since
   `gatk-launch` invokes the `spark-submit` tool behind-the-scenes.
* Note that the examples above use YARN but we have successfully run GATK4 on Mesos as well.

#### <a name="dataproc">Running GATK4 Spark tools on Google Cloud Dataproc:</a>
  * You must have a [Google cloud services](https://cloud.google.com/) account, and have spun up a Dataproc cluster
    in the [Google Developer's console](https://console.developers.google.com). You may need to have the "Allow API access to all Google Cloud services in the same project" option enabled (settable when you create a cluster).
  * You need to have installed the Google Cloud SDK from [here](https://cloud.google.com/sdk/), since
    `gatk-launch` invokes the `gcloud` tool behind-the-scenes. As part of the installation, be sure
      that you follow the `gcloud` setup instructions [here](https://cloud.google.com/sdk/gcloud/).
  * Your inputs to the GATK when running on dataproc are typically in Google Cloud Storage buckets, and should be specified on
    your GATK command line using the syntax `gs://my-gcs-bucket/path/to/my-file`
  * You can run GATK4 jobs on Dataproc from your local computer or from the VM (master node) on the cloud.

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
  * When using Dataproc you can access the web interfaces for YARN, Hadoop and HDFS. Follow [these instructions](https://cloud.google.com/dataproc/cluster-web-interfaces) to create an SSH tunnel and connect with your browser.
  * Note that the spark-specific arguments are separated from the tool-specific arguments by a `--`.
  * If you want to avoid uploading the GATK jar to GCS on every run, set the `GATK_GCS_STAGING`
    environment variable to a bucket you have write access to (eg., `export GATK_GCS_STAGING=gs://<my_bucket>/`)
  * Dataproc Spark clusters are configured with [dynamic allocation](https://spark.apache.org/docs/latest/job-scheduling.html#dynamic-resource-allocation) so you can omit the "--num-executors" argument and let YARN handle it automatically.

#### <a name="2bit">Note on 2bit Reference</a>
* Note: Some GATK Spark tools by default require the reference file to be in 2bit format (notably `BaseRecalibratorSpark`,`BQSRPipelineSpark` and `ReadsPipelineSpark`). You can convert your fasta to 2bit by using the `faToTwoBit` utility from [UCSC](http://hgdownload.soe.ucsc.edu/admin/exe/) - see also the [documentation for `faToTwoBit`](https://genome.ucsc.edu/goldenpath/help/blatSpec.html#faToTwoBitUsage).

#### <a name="R">Using R to generate plots</a>
Certain GATK tools may optionally generate plots if R is installed.  We recommend **R v3.1.3** if you want to produce plots.  If you are uninterested in plotting, R is still required by several of the unit tests.  Plotting is currently untested and should be viewed as a convenience rather than a primary output.

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

#### <a name="cnv_workflows">Running the CNV workflows</a>

* A walkthrough and examples for the CNV workflows can be found [here](http://gatkforums.broadinstitute.org/gatk/discussion/9143)

## <a name="developers">For GATK Developers</a>

#### <a name="dev_guidelines">General guidelines for GATK4 developers</a>

* **Do not put private or restricted data into the repo.**

* **Try to keep datafiles under 100kb in size.** Larger test files should go into `src/test/resources/large` (and subdirectories) so that they'll be stored and tracked by git-lfs as described [above](#lfs).

* GATK4 is BSD licensed.  The license is in the top level LICENSE.TXT file.  Do not add any additional license text or accept files with a license included in them.

* Each tool should have at least one good end-to-end integration test with a check for expected output, plus high-quality unit tests for all non-trivial utility methods/classes used by the tool. Although we have no specific coverage target, coverage should be extensive enough that if tests pass, the tool is guaranteed to be in a usable state.

* All newly written code must have good test coverage (>90%).

* All bug fixes must be accompanied by a regression test.

* All pull requests must be reviewed before merging to master (even documentation changes).

* Don't issue or accept pull requests that introduce warnings. Warnings must be addressed or suppressed.

* Don't issue or accept pull requests that significantly decrease coverage (less than 1% decrease is sort of tolerable). 

* Don't use `toString()` for anything other than human consumption (ie. don't base the logic of your code on results of `toString()`.)

* Don't override `clone()` unless you really know what you're doing. If you do override it, document thoroughly. Otherwise, prefer other means of making copies of objects.

* For logging, use [org.apache.logging.log4j.Logger](https://logging.apache.org/log4j/2.0/log4j-api/apidocs/org/apache/logging/log4j/Logger.html)

* We mostly follow the [Google Java Style guide](http://google-styleguide.googlecode.com/svn/trunk/javaguide.html)

* Git: Don't push directly to master - make a pull request instead. 

* Git: Rebase and squash commits when merging.

* If you push to master or mess up the commit history, you owe us 1 growler or tasty snacks at happy hour. If you break the master build, you owe 3 growlers (or lots of tasty snacks). Beer may be replaced by wine (in the color and vintage of buyer's choosing) in proportions of 1 growler = 1 bottle. 

#### <a name="testing">Testing GATK</a>

* Before running the test suite, be sure that you've installed `git lfs` and downloaded the large test data, following the [git lfs setup instructions](#lfs)

* To run the test suite, run **`./gradlew test`**.
    * Test report is in `build/reports/tests/test/index.html`.
    * What will happen depends on the value of the `TEST_TYPE` environment variable: 
       * unset or any other value         : run non-cloud unit and integration tests, this is the default
       * `cloud`, `unit`, `integration`, `spark`   : run only the cloud, unit, integration, or Spark tests
       * `all`                            : run the entire test suite
    * Cloud tests require being logged into `gcloud` and authenticated with a project that has access
      to the cloud test data.  They also require setting several certain environment variables.
      * `HELLBENDER_JSON_SERVICE_ACCOUNT_KEY` : path to a local JSON file with [service account credentials](https://cloud.google.com/storage/docs/authentication#service_accounts) 
      * `HELLBENDER_TEST_PROJECT` : your google cloud project 
      * `HELLBENDER_TEST_APIKEY` : your google cloud API key
      * `HELLBENDER_TEST_STAGING` : a gs:// path to a writable location
      * `HELLBENDER_TEST_INPUTS` : path to cloud test data, ex: gs://hellbender/test/resources/ 
    * Setting the environment variable `TEST_VERBOSITY=minimal` will produce much less output from the test suite 

* To run a subset of tests, use gradle's test filtering (see [gradle doc](https://docs.gradle.org/current/userguide/java_plugin.html)):
    * You can use `test.single` when you just want to run a specific test class:
        * `./gradlew test -Dtest.single=SomeSpecificTestClass`
    * You can also use `--tests` with a wildcard to run a specific test class, method, or to select multiple test classes:
        * `./gradlew test --tests *SomeSpecificTestClass`
        * `./gradlew test --tests *SomeTest.someSpecificTestMethod`
        * `./gradlew test --tests all.in.specific.package*`

* To run tests and compute coverage reports, run **`./gradlew jacocoTestReport`**. The report is then in `build/reports/jacoco/test/html/index.html`.
  (IntelliJ has a good coverage tool that is preferable for development).

* We use [Travis-CI](https://travis-ci.org/broadinstitute/gatk) as our continuous integration provider.

    * Before merging any branch make sure that all required tests pass on travis.
    * Every travis build will upload the test results to our GATK Google Cloud Storage bucket.
      A link to the uploaded report will appear at the very bottom of the travis log.
      Look for the line that says `See the test report at`.
      If TestNG itself crashes there will be no report generated.

* We use [Broad Jenkins](https://gatk-jenkins.broadinstitute.org/view/Performance/) for our long-running tests and performance tests.
    * To add a performance test (requires Broad-ID), you need to make a "new item" in Jenkins and make it a "copy" instead of a blank project. You need to base it on either the "-spark-" jobs or the other kind of jobs and alter the commandline. 

* To output stack traces for `UserException` set the environment variable `GATK_STACKTRACE_ON_USER_EXCEPTION=true`

#### <a name="lfs">Using Git LFS to download and track large test data</a>

We use [git-lfs](https://git-lfs.github.com/) to version and distribute test data that is too large to check into our repository directly. You must install and configure it in order to be able to run our test suite.

* After installing [git-lfs](https://git-lfs.github.com/), run `git lfs install`
    * This adds hooks to your git configuration that will cause git-lfs files to be checked out for you automatically in the future.
    
* To manually retrieve the large test data, run `git lfs pull` from the root of your GATK git clone.
    * The download is several hundred megabytes.
    
* To add a new large file to be tracked by git-lfs, simply:
    * Put the new file(s) in `src/test/resources/large` (or a subdirectory)
    * `git add` the file(s), then `git commit -a`
    * That's it! Do ***not*** run `git lfs track` on the files manually: all files in `src/test/resources/large` are tracked by git-lfs automatically. 

#### <a name="intellij">Creating a GATK project in the IntelliJ IDE (last tested with version 2016.2.4):</a>

* Ensure that you have `gradle` and the Java 8 JDK installed

* You may need to install the TestNG and Gradle plugins (in preferences)

* Clone the GATK repository using git

* In IntelliJ, click on "Import Project" in the home screen or go to File -> New... -> Project From Existing Sources...

* Select the root directory of your GATK clone, then click on "OK"

* Select "Import project from external model", then "Gradle", then click on "Next"

* Ensure that "Gradle project" points to the build.gradle file in the root of your GATK clone

* Select "Use auto-import" and "Use default gradle wrapper".

* Make sure the Gradle JVM points to Java 1.8 

* Click "Finish"

* After downloading project dependencies, IntelliJ should open a new window with your GATK project

* Make sure that the Java version is set correctly by going to File -> "Project Structure" -> "Project". Check that the "Project SDK" is set to your Java 1.8 JDK, and "Project language level" to 8 (you may need to add your Java 8 JDK under "Platform Settings" -> SDKs if it isn't there already). Then click "Apply"/"Ok".

#### <a name="debugging">Setting up debugging in IntelliJ</a>

* Follow the instructions above for creating an IntelliJ project for GATK

* Go to Run -> "Edit Configurations", then click "+" and add a new "Application" configuration

* Set the name of the new configuration to something like "GATK debug"

* For "Main class", enter `org.broadinstitute.hellbender.Main`

* Ensure that "Use classpath of module:" is set to use the "gatk" module's classpath

* Enter the arguments for the command you want to debug in "Program Arguments"

* Click "Apply"/"Ok"

* Set breakpoints, etc., as desired, then select "Run" -> "Debug" -> "GATK debug" to start your debugging session

* In future debugging sessions, you can simply adjust the "Program Arguments" in the "GATK debug" configuration as needed

#### <a name="intellij_gradle_refresh">Updating the Intellij project when dependencies change</a>
If there are dependency changes in `build.gradle` it is necessary to refresh the gradle project. This is easily done with the following steps.

* Open the gradle tool window  ( "View" -> "Tool Windows" -> "Gradle" )
* Click the refresh button in the Gradle tool window.  It is in the top left of the gradle view and is represented by two blue arrows.

#### <a name="jprofiler">Setting up profiling using JProfiler</a>

   * Running JProfiler standalone:
       * Build a full GATK4 jar using `./gradlew localJar`
       * In the "Session Settings" window, select the GATK4 jar, eg. `~/gatk/build/libs/gatk-package-4.alpha-196-gb542813-SNAPSHOT-local.jar` for "Main class or executable JAR" and enter the right "Arguments"
       * Under "Profiling Settings", select "sampling" as the "Method call recording" method.

   * Running JProfiler from within IntelliJ:
       * JProfiler has great integration with IntelliJ (we're using IntelliJ Ultimate edition) so the setup is trivial.   
       * Follow the instructions [above](#intellij) for creating an IntelliJ project for GATK  
       * Right click on a test method/class/package and select "Profile" 

#### <a name="sonatype">Uploading Archives to Sonatype (to make them available via maven central)</a>
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

#### <a name="docker_building">Building GATK4 Docker images</a>

Please see the [the Docker README](scripts/docker/README.md) in ``scripts/docker``.  This has instructions for the Dockerfile in the root directory.

#### <a name="releasing_gatk">Releasing GATK4</a>

Please see the [How to release GATK4](https://github.com/broadinstitute/gatk/wiki/How-to-release-GATK4) wiki article for instructions on releasing GATK4.

#### <a name="gatkdocs">Generating GATK4 documentation</a>

To generate GATK documentation, run `./gradlew gatkDoc`

* Generated docs will be in the `build/docs/gatkdoc` directory.

#### <a name="zenhub">Using Zenhub to track github issues</a>

We use [Zenhub](https://www.zenhub.com/) to organize and track github issues.

* To add Zenhub to github, go to the [Zenhub home page](https://www.zenhub.com/) while logged in to github, and click "Add Zenhub to Github"

* Zenhub allows the GATK development team to assign time estimates to issues, and to mark issues as Triaged/In Progress/In Review/Blocked/etc.

## <a name="spark_further_reading">Further Reading on Spark</a>

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

## <a name="contribute">How to contribute to GATK</a>
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

## <a name="discussions">Discussions</a>
* [GATK forum](http://gatkforums.broadinstitute.org/) for general discussions on how to use the GATK and support questions.
* [Issue tracker](https://github.com/broadinstitute/gatk/issues) to report errors and enhancement ideas. 
* Discussions also take place in [GATK pull requests](https://github.com/broadinstitute/gatk/pulls)

## <a name="authors">Authors</a>
The authors list is maintained in the [AUTHORS](https://github.com/broadinstitute/gatk/edit/master/AUTHORS) file. 
See also the [Contributors](https://github.com/broadinstitute/gatk/graphs/contributors) list at github. 

## <a name="license">License</a>
Licensed under the BSD License. See the [LICENSE.txt](https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT) file.
