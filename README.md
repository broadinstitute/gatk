[![Build Status](https://magnum.travis-ci.com/broadinstitute/hellbender.svg?token=WFzCX7pDpMhnHx5RX8kq&branch=master)](https://magnum.travis-ci.com/broadinstitute/hellbender)

Hellbender
================

The (eventually public) parts of the next generation of GATK/Picard/Prometheus methods engine and tools.


To build, run `gradle install`.

To run the main program, run `build/install/hellbender/bin/hellbender`.

To run tests and compute coverage reports, run `gradle jacocoTestReport`. The report is then in `build/reports/jacoco/test/html/index.html`.

----------------
Notes to devs

* Untested code is considered **non existent** and thus is subject to removal at any time (exception is main methods, or super corner conditions, or `toString()` code). New code without tests should not be accepted at pull requests. 

* Don't use toString for anything other than human consumption (ie dont base your code on it.)

* Do use Java 8 features liberally, in particular streams, functional programming etc.

* Don't override `clone()` unless you really know what you're doing. If you do override it, document thoroughly to show that you do know what you're doing. 


