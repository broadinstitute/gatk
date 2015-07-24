[![Build Status](https://travis-ci.org/broadinstitute/hellbender.svg?branch=master)](https://travis-ci.org/broadinstitute/hellbender)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/hellbender/badge.svg?branch=master)](https://coveralls.io/r/broadinstitute/hellbender?branch=master)

**This project is in an early stage of development.  It is subject to change without warning.**  

Hellbender
================

The public parts of the next generation of GATK/Picard methods engine and tools.


Requirements
------------
* Java 8

* Gradle 2.2.1

* R 3.1.3 


Installation
------------
To build and run all tests, run `gradle check`. Test report is in `build/reports/tests/index.html`.

To only build, run `gradle installApp`.

To run all tests, run `gradle test`. 

To run a single test class, run something like this, `gradle test -Dtest.single=ReadUtilsUnitTest`.

To run tests and compute coverage reports, run `gradle jacocoTestReport`. The report is then in `build/reports/jacoco/test/html/index.html`. (IntelliJ 14 has a good coverage tool that is preferable for development).

To run the main program, run `build/install/hellbender/bin/hellbender`.

Note: for faster gradle operations, add `org.gradle.daemon=true` to your `~/.gradle/gradle.properties` file.  This will keep a gradle daemon running in the background and avoid the ~6s gradle start up time on every command.  

General guidelines for Hellbender developers
----------------

* **Do not put private or restricted data into the repo.**

* **Try to keep datafiles under 100kb in size.**

* Hellbender is  BSD licensed.  The license is in the top level LICENSE.TXT file.  Do not add any additional license text or accept files with a license included in them.

* Hellbender is written in Java 8. Enjoy the new features in particular streams, lambdas, methods in interfaces.

* Untested code is considered **non-existent** and thus is subject to removal at any time (exception is main methods, or super corner conditions, or `toString()` code).

* Each tool should have at least one good end-to-end integration test with a check for expected output, plus high-quality unit tests for all non-trivial utility methods/classes used by the tool. Although we have no specific coverage target, coverage should be extensive enough that if tests pass, the tool is guaranteed to be in a usable state.

* All newly written code must have good test coverage (>90%).

* Don't issue or accept pull requests that introduce warnings. Warnings need to be addressed or suppressed.

* Don't issue or accept pull requests that significantly decrease coverage (less than 1% decrease is sort of tolerable). 

* Don't override `clone()` unless you really know what you're doing. If you do override it, document thoroughly. Otherwise, prefer other means of making copies of objects.

* Don't use `toString()` for anything other than human consumption (ie. don't base the logic of your code on results of `toString()`.)

* We mostly follow the Google Java Style guide: http://google-styleguide.googlecode.com/svn/trunk/javaguide.html

* Git: Don't push directly to master - make a pull request instead. 

* Git: Rebase and squash commits when merging.

* If you push to master or mess the commit history, you owe us 1 growler or tasty snacks at happy hour. If you break the master build, you owe 3 growlers (or lots of tasty snacks). Beer may be replaced by wine (in the color and vintage of buyer's choosing) in proportions of 1 growler = 1 bottle. 

Tests
----------------
We use [Travis-CI](https://travis-ci.org/broadinstitute/hellbender) as our continuous integration provider.

* Before merging any branch make sure that all required tests pass on travis.  Currently tests in the `cloud` and `bucket` groups as well as the tests running on spark are not required to pass before merging.
* Every travis build will upload the test results to our hellbender google bucket.  A link to the uploaded report will appear at the very bottom of the travis log.  Look for the line that says `See the test report at`.
If TestNG itself crashes there will be no report generated.

R Dependency
----------------
Certain Hellbender tools may optionally generate plots if R is installed.  We recommend **R v3.1.3** if you want to produce plots.  If you are uninterested in plotting, R is still required by several of the unit tests.  Plotting is currently untested and should be viewed as a convinience rather than a primary output.  

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

Creating a Hellbender project in the IntelliJ IDE:
----------------

* Ensure that you have `gradle` and the Java 8 JDK installed

* Install the TestNG plugin (in preferences)

* Clone the Hellbender repository using git

* In IntelliJ, go to File -> "Import Project"

* Select the root directory of your Hellbender clone, then "Ok"

* Select "Import project from external model", then "Gradle", then "Next"

* Ensure that "Gradle project" points to the build.gradle file in the root of your Hellbender clone

* Select "Use local gradle distribution", and enter your Gradle home directory in the "Gradle home" box. This will be the directory that contains the `bin` directory where the **actual** gradle executable (**not** merely a symlink to it!) lives. For example, if the actual gradle executable is `/usr/local/Cellar/gradle/2.2.1/libexec/bin/gradle`, you would enter `/usr/local/Cellar/gradle/2.2.1/libexec/` as your gradle home directory.

* Click "Finish"

* After downloading project dependencies, IntelliJ should open a new window with your Hellbender project

* In File -> "Project Structure" -> "Project", set the "Project SDK" to your Java 1.8 JDK, and "Project language level" to 8 (you may need to add your Java 8 JDK under "Platform Settings" -> SDKs if it isn't there already). Then click "Apply"/"Ok".


Setting up debugging in IntelliJ
----------------

* Follow the instructions above for creating an IntelliJ project for Hellbender

* Go to Run -> "Edit Configurations", then click "+" and add a new "Application" configuration

* Set the name of the new configuration to something like "Hellbender debug"

* For "Main class", enter `org.broadinstitute.hellbender.Main`

* Ensure that "Use classpath of module:" is set to use the "hellbender" module's classpath

* Enter the Hellbender arguments for the command you want to debug in "Program Arguments"

* Click "Apply"/"Ok"

* Set breakpoints, etc., as desired, then select "Run" -> "Debug" -> "Hellbender debug" to start your debugging session

* In future debugging sessions, you can simply adjust the "Program Arguments" in the "Hellbender debug" configuration as needed

Updating the Intellij project when dependencies change
-------------------
If there are dependency changes in `build.gradle` it is necessary to refresh the gradle project. This is easily done with the following steps.

* Open the gradle tool window  ( "View" -> "Tool Windows" -> "Gradle" )
* Click the refresh button in the Gradle tool window.  It is in the top left of the gradle view and is represented by two blue arrows.

Uploading Archives to Sonatype
-------------------
To upload snapshots to sonatype you'll need the following:

*You must have a registered account on the sonatype JIRA (and be approved as a hellbender uploader)
*You need to configure several additional properties in your `/~.gradle/gradle.properties` file

*If you want to upload a release instead of a snapshot you will additionally need to have access to the hellbender signing key and password

```
#needed for snapshot upload
sonatypeUsername=<your sonatype username>
sonatypePassword=<your sonatype password>

#needed for signing a release
signing.keyId=<hellbender key id>
signing.password=<hellbender key password>
signing.secretKeyRingFile=/Users/<username>/.gnupg/secring.gpg
```

To perform an upload, use
```
gradle uploadArchives
```

Currently all builds are considered snapshots.  The archive name is based off of `git describe`.

