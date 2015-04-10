[![Build Status](https://travis-ci.org/broadinstitute/hellbender.svg?branch=master)](https://travis-ci.org/broadinstitute/hellbender)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/hellbender/badge.svg?branch=master)](https://coveralls.io/r/broadinstitute/hellbender?branch=master)

**This project is in an early stage of development.  It is subject to change without warning.**  

Hellbender
================

The public parts of the next generation of GATK/Picard/Prometheus methods engine and tools.


To build, run `gradle install` (gradle 2.2.1 recommended).

To run the main program, run `build/install/hellbender/bin/hellbender`.

To run all tests, run `gradle test`. Report is in `build/reports/tests/index.html`.

To run a single test class, run something like this, `gradle test -Dtest.single=ReadUtilsUnitTest`.

To run tests and compute coverage reports, run `gradle jacocoTestReport`. The report is then in `build/reports/jacoco/test/html/index.html`. (IntelliJ 14 has a good coverage tool that is preferable for development).

For faster gradle operations, add `org.gradle.daemon=true` to your `~/.gradle/gradle.properties` file.  This will keep a gradle daemon running in the background and avoid the ~6s gradle start up time on every command.  

General guidelines for Hellbender developers
----------------

* **Do not put private or restricted data into the repo.**

* **Try to keep datafiles under 100kb in size.**

* Hellbender is  BSD licensed.  The license is in the top level LICENSE.TXT file.  Do not add any additional license text or accept files with a license included in them.

* Untested code is considered **non-existent** and thus is subject to removal at any time (exception is main methods, or super corner conditions, or `toString()` code).

* Don't accept pull requests that introduce warnings. Warnings need to be addressed or suppressed.

* Don't accept pull requests that significantly decrease coverage (less than 1% decrease is sort of tolerable). 

* All newly written code must have good test coverage (>90%).

* Hellbender is written in Java 8. Enjoy the new features in particular streams, lambdas, methods in interfaces.

* Don't override `clone()` unless you really know what you're doing. If you do override it, document thoroughly. Otherwise, prefer other means of making copies of objects.

* Don't override or call `finalize()` even when you think you know what you're doing.

* Don't use `toString()` for anything other than human consumption (ie. don't base your code on it.)

* We mostly follow the Google Java Style guide: http://google-styleguide.googlecode.com/svn/trunk/javaguide.html

R Dependency
----------------
Certain Hellbender tools rely on R for generating plots. **R v3.1.3** must be installed for all unit tests to work (other versions of R may work, but there is no guarantee).  

R installation is not part of the gradle build.  See http://cran.r-project.org/ for general information on installing R for your system.
* for [ubuntu specific instructions](http://cran.r-project.org/bin/linux/ubuntu/README))
* for OSX we recommend installation through [homebrew](http://brew.sh/)

* For more information on R see http://www.ihater.org/

Gradle installApp will attempt to install R packages, this may potentially fail if your R library is not writable.  If this happens, the package installation script is the package installation script is `scripts/install_R_packages.R`.  Either run it as superuser to force installation into the sites library or run interactively and create a local library.
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
