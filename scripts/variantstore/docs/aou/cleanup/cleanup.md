# AoU Callset Cleanup

## Overview

As a general rule, any artifacts that have clearly become obsolete (e.g. VDSes with known issues that have been
superseded by corrected versions, obsolete sets of prepare tables, etc.) should be deleted ASAP. If it's not clear
whether an artifact should be cleaned up or not, [calculate the monthly cost to preserve the artifact](Cost.md)
(i.e. the sum of all relevant GCS or BigQuery storage costs) as well as the cost to regenerate the artifact.
Reach out to leadership with these numbers for their verdict on whether to keep or delete.

## Specific AoU GVS Artifacts

During the course of creating AoU callsets several large and expensive artifacts are created:

* Pilot workspace / dataset
    * for the Delta callset the Variants team created an AoU 10K workspace and dataset to pilot the Hail/VDS creation
    * the dataset, workspace data tables, and submission files have been deleted to save money, but the workspace has been kept around for future testing
* Production BigQuery dataset
    * for each previous callset, there was (at least) one new dataset created
    * the dream is to keep the same dataset for multiple callsets and just add new samples, regenerate the filter and create new deliverables, but that has yet to happen because of new features requested for each callset (e.g. update to Dragen version, addition of ploidy data, different requirements to use Hail...etc.)
    * if there are datasets from previous callsets that aren't needed anymore (check with AoU and Lee/Wail), they should be deleted
* Prepare tables
    * needed for VCF or PGEN extract
    * only variant tables are used by `GvsCallsetStatistics.wdl` for callset statistics deliverable
    * These tables are created with a TTL by default, but if they are not needed anymore, definitely delete ASAP as these are huge and expensive.
      Please be particularly mindful to delete these tables in a timely manner if the TTL has been removed for whatever reason.

## Submissions directory cleanup

