# M2 Automated Tests for WDL

**This directory is for GATK devs only**

This directory contains scripts for running M2 WDL tests in the automated travis build environment.

Please note that this tests whether the WDL will complete successfully.

# FAQ
## Does this test Oncotator?

Yes, BUT it tests it without any datasources.

## Does this test FilterByOrientationBias?

Yes.

For C/T (FFPE deamination) and G/T (OxoG)

## Does this test the GATK4 jar override parameter (gatk4_jar_override)?

No

## Does this test filtering of common germline variants?

Only dbSNP

## Does this test the contamination portion of the WDL?

No

## Does this test that results make sense if the workflow completes?

No, this is the wrong place for such a test.

## Does this test the sensitivity and precision?

No, this is the wrong place for such a test.

## Does this test tumor only?

Yes, with the same caveats as above.