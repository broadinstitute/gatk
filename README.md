[![Build Status](https://magnum.travis-ci.com/broadinstitute/hellbender.svg?token=WFzCX7pDpMhnHx5RX8kq&branch=master)](https://magnum.travis-ci.com/broadinstitute/hellbender)

Hellbender
================

The (eventually public) parts of the next generation of GATK/Picard/Prometheus methods engine and tools.


To build, run `gradle build`.

To run the main program, run `build/install/hellbender/bin/hellbender`.

To run tests and compute coverage reports, run `gradle jacocoTestReport`. The report is then in `build/reports/jacoco/test/html/index.html`.

Note to devs: untested code is consider **non existent** and thus is subject to removal at any time (exception is main methods). New code without tests should not be accepted at pull requests.
