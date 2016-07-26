[![Build Status](https://travis-ci.org/broadinstitute/gatk-protected.svg?branch=master)](https://travis-ci.org/broadinstitute/gatk-protected)
[![Coverage Status](https://coveralls.io/repos/broadinstitute/gatk-protected/badge.svg?branch=master&service=github)](https://coveralls.io/github/broadinstitute/gatk-protected?branch=master)

GATK4-Protected (codename Hellbender-protected)
===============================================

GATK4 development of the license-protected part of the toolkit

This README is aimed at developers.  For user information, please see the [GATK 4 forum](http://http://gatkforums.broadinstitute.org/gatk/categories/gatk-4-alpha)

Requirements
------------
* R 3.1.3 see additional requirements below: [R package requirements](#r-required-packages)

* Java 8

* (Developers) Gradle 2.13 is needed for building the GATK. We recommend using the `./gradlew` script which will
download and use an appropriate gradle version automatically.

* (Developers) git lfs 1.1.0 (or greater) is needed for testing GATK-Protected builds.  It is needed to download large files for the complete test suite. Run ``git lfs install`` after downloading, followed by ``git lfs pull`` to download the large files. The download is ~500 MB.

Read GATK 4 README
------------------------

Please refer to the GATK 4 public repo readme file for general guidelines and how to set-up your development environment:

https://github.com/broadinstitute/hellbender/blob/master/README.md


#### R Required Packages
R packages can be installed using the install_R_packages.R script inside the scripts directory.

### Create the jar file

`` ./gradlew shadowJar ``

A jar file will appear in ``build/libs``.

The CNV case and PoN workflows (description and examples)
---------------------------------------------------------

This can be found [here](http://gatkforums.broadinstitute.org/gatk/discussion/6791/description-and-examples-of-the-steps-in-the-cnv-case-and-cnv-pon-creation-workflows)


Running the CNV case and PoN creation Workflows with premade Queue scripts
--------------------------------------------------------------------------

For [Broad Internal instructions](http://gatkforums.broadinstitute.org/gatk/discussion/6786/howto-run-gatk-cnv-using-premade-queue-scripts-broad-internal)

