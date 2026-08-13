# AoU Callset Cleanup

## Overview

As a general rule, any AoU callset artifacts that have clearly become obsolete (e.g. VDSes with known issues that have
been superseded by corrected versions, obsolete sets of prepare tables, etc.) should be deleted ASAP. If it's not clear
whether an artifact should be cleaned up or not, [calculate the monthly cost to preserve the artifact](Cost.md)
(i.e. the sum of all relevant GCS or BigQuery storage costs) as well as the cost to regenerate the artifact.
Reach out to leadership with these numbers for their verdict on whether to keep or delete.

## Specific AoU GVS Artifacts

During the course of creating AoU callsets several large and expensive artifacts are created:

* Pilot workspace / dataset
    * For the Delta callset the Variants team created an AoU 10K workspace and dataset to pilot Hail/VDS creation.
      Similarly for Foxtrot there was a "Foxtrot Batch Testing" workspace created to pilot our usage of the then-new
      Google Batch API at AoU scale.
    * These datasets, workspace data tables, and submission files have been cleaned up to save money, but the workspaces
      have been kept around for future reference.
* Production BigQuery dataset
    * For callsets up to and including Echo, at least one completely new dataset was created. For Foxtrot, we copied
      the Echo dataset as the starting point for the Foxtrot dataset.
    * If there are datasets from previous callsets that aren't needed anymore (check with Lee/Liz), they should be
      deleted. Currently we still have the Delta, Echo, and Foxtrot datasets.
* Prepare tables
    * Needed for VCF or PGEN extract.
    * Only the variant data prepare tables are used by `GvsCallsetStatistics.wdl` for the callset statistics deliverable.
    * These tables are created with a TTL by default, but if they are not needed anymore please delete ASAP as these are
      usually huge and expensive. Please be particularly mindful to delete prepare tables in a timely manner if the TTL
      has been removed for whatever reason.

## Submissions directory cleanup

When cleaning up some or all submissions in a workspace, the following commands can be used to remove all files except
for the relatively small but important "detritus" files that can be useful for future inspection. All commands below
can be run in a terminal running in the appropriate AoU workspace. First enumerate all files under a submissions path:

```
gcloud storage ls -r gs://<workspace bucket>/submissions[/<optional submission_id>] > submission_files.txt
```

Without the `[/<optional submission_id>]` bit, all submissions in the workspace bucket would be enumerated for cleanup.

Clean up selectively:

```
grep -E -v '/rc$|/stdout$|/stderr$|\.log$|^$|:$' submission_files.txt | (gcloud storage rm -I 2>/dev/null) &
```

The `grep` above exempts the `rc`, `stdout`, `stderr` and `*.log` files from cleanup, and ignores blank and "directory"
line artifacts from the preceding `ls -r` step.
