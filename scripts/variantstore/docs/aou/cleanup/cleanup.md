# AoU Callset Cleanup

## Overview

The current Variants policy for AoU callsets is effectively to retain all versions of all artifacts forever. As the  storage costs for these artifacts can be significant (particularly from Delta onward), the Variants team would like to make the cost of retaining artifacts more clear so conscious choices can be made about what to keep and what to delete.

As a general rule, any artifacts that have clearly become obsolete (e.g. VDSes with known issues that have been  superseded by corrected versions, obsolete sets of prepare tables, etc.) should be deleted ASAP. If it's not clear to the Variants team whether an artifact should be cleaned up or not, [we should calculate the monthly cost to preserve the artifact (e.g. the sum of all relevant GCS or BigQuery storage costs) as well as the cost to regenerate the artifact](Cost.md).

Reach out to leadership with these numbers for his verdict on whether to keep or delete.

## Specific AoU GVS Artifacts

During the course of creating AoU callsets several large and expensive artifacts are created:

* Pilot workspace / dataset
    * for the Delta callset the Variants team created an AoU 10K workspace and dataset to pilot the Hail/VDS creation
    * these processes will mature to the point where some or all of the contents of this workspace and dataset should be deleted
* Production BigQuery dataset
    * for each previous callset, there was (at least) one new dataset created
    * the dream is to keep the same dataset for multiple callsets and just add new samples, regenerate the filter and create new deliverables, but that has yet to happen because of new features requested for each callset (e.g. update to Dragen version, addition of ploidy data, different requirements to use Hail...etc.)
* Prepare tables
    * needed for VCF or PGEN extract
    * only variant tables are used by `GvsCallsetStatistics.wdl` for callset statistics deliverable
* Sites-only VCFs
    * the VAT is created from a sites-only VCF and its creation is the most resource-intensive part of the VAT pipeline
    * clean up any failed runs once the VAT has been delivered and accepted
* Avro files (Delta onward)
    * huge, several times larger than the corresponding Hail VDS
    * as long as the VDS has been delivered and accepted, they can be deleted
* Hail VariantDataset (VDS)
    * we used to create it in the Terra workspace and then copy it
    * WDL was updated to create in the "deliverables" GCS bucket so there is only one copy of each one 
    * clean up any failed runs once the VDS has been delivered and accepted
* PGEN/VCF Intermediate Files
    * PGEN: multiple versions of the PGEN files are created by GvsExtractCallsetPgenMerged.wdl because it delivers files split by chromosome
    * VCF: only one version of the VCF files and indices are created, but check for failed runs

## Internal sign-off protocol

The Variants team currently has the following VDS internal sign-off protocol:

1. Generate a VDS for the candidate callset into the "delivery" bucket.
1. Open up the VDS in a beefy notebook and confirm the "shape" looks right.
1. Run `GvsPrepareRangesCallset.wdl` to generate a prepare table of VET data
1. Run `GvsCallsetStatistics.wdl` to generate callset statistics for the candidate callset using the prepare VET created in the preceding step
1. Copy the output of `GvsCallsetStatistics.wdl` into the "delivery" bucket.
1. Email the paths to the VDS and callset statistics to Lee/Wail for QA / approval

