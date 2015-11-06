[![Build Status](https://travis-ci.org/broadinstitute/gatk.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gatk/badge.svg?branch=master)](https://coveralls.io/r/broadinstitute/gatk?branch=master)

**This project is in an early stage of development.  It is subject to change without warning. Do not use this code for production work.**  

###GATK 4 (codename Hellbender)

The public parts of the next generation of GATK/Picard methods engine and tools.

This project is in a pre-alpha development stage and is not yet ready for general use. 

If you are looking for the current version of GATK to use in production work, please see the [GATK website](http://www.broadinstitute.org/gatk), where you can download a precompiled executable, read documentation, ask questions and receive technical support.

If you are looking for the codebase of the current production version of GATK, please see either the [GATK development framework repository](https://github.com/broadgsa/gatk/) or the [full GATK tools repository](https://github.com/broadgsa/gatk-protected/).


##Requirements
* Java 8

* Gradle 2.7

* R 3.1.3

* Git 2.5

* The [git-lfs plugin](https://git-lfs.github.com/). Run `git lfs init` after installing, followed by `git lfs pull` to download the large files.  We recommend v0.6.0 but newer versions may also work.


##Installation
To build and run all tests, run `gradle check`. Test report is in `build/reports/tests/index.html`.

To only build, run `gradle installDist`.

To run all tests, run `gradle test`. 
What will happen depends on the value of the `CLOUD` environment variable: if it's `false` or unset then only local tests are run. If it's `mandatory` then it'll run only the cloud tests. 

To run a single test class, run something like this, `gradle test -Dtest.single=ReadUtilsUnitTest`.

To run tests and compute coverage reports, run `gradle jacocoTestReport`. The report is then in `build/reports/jacoco/test/html/index.html`. (IntelliJ 14 has a good coverage tool that is preferable for development).


Note: for faster gradle operations, add `org.gradle.daemon=true` to your `~/.gradle/gradle.properties` file.  This will keep a gradle daemon running in the background and avoid the ~6s gradle start up time on every command.  


##Running GATK4
To run the main program locally, run `build/install/gatk/bin/gatk`.

For spark tools, first build the Spark jar using `gradle sparkJar`. The jar can be found in `build/libs/` and will have
a name that matches the pattern `gatk-all-4.pre-alpha-7-*-SNAPSHOT-spark.jar`, where `*` is the
hash of the current commit.

There are three main ways to run GATK4 with Spark,
* Locally, using the regular gatk jar and `--sparkMaster 'local[*]'`

* On an on-premises cluster by copying the spark jar to the master and running `spark-submit`

* On Google Cloud Dataproc. For detailed instructions, see our [Clould Dataproc wiki page](https://github.com/broadinstitute/gatk/wiki/Running-GATK-on-Cloud-Dataproc)

##General guidelines for GATK4 developers

* **Do not put private or restricted data into the repo.**

* **Try to keep datafiles under 100kb in size.** Larger test files should go into `src/test/resources/large`, and must be
  managed using `git lfs` by running `git lfs track <file>` on each new large file before commit.

* GATK4 is  BSD licensed.  The license is in the top level LICENSE.TXT file.  Do not add any additional license text or accept files with a license included in them.

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

##Tests
We use [Travis-CI](https://travis-ci.org/broadinstitute/gatk) as our continuous integration provider.

* Before merging any branch make sure that all required tests pass on travis.
* Every travis build will upload the test results to our gatk google bucket.  A link to the uploaded report will appear at the very bottom of the travis log.  Look for the line that says `See the test report at`.
If TestNG itself crashes there will be no report generated.

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
gradle uploadArchives
```

Currently all builds are considered snapshots.  The archive name is based off of `git describe`.

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
* Make sure that your code passes all our tests. You can run the tests with `gradle test` in the root directory.
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
