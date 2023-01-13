# AoU callset cleanup

## Overview

The current Variants policy for AoU callsets is effectively to retain all versions of all artifacts forever. As the
storage costs for these artifacts can be significant (particularly from Delta onward), the Variants team would like to
make the cost of retaining artifacts more clear so conscious choices can be made about what to keep and what to delete.

As a general rule, any artifacts that have clearly become obsolete (e.g. VDSes with known issues that have been
superseded by corrected versions, obsolete sets of "prepare" tables, etc.) should be deleted ASAP. If it's not clear to
the Variants team whether an artifact should be cleaned up or not, we should calculate the monthly cost to preserve the
artifact (e.g. the sum of all relevant GCS or BigQuery storage costs) as well as the cost to regenerate the artifact.
Reach out to Lee with these numbers for his verdict on whether to keep or delete.

## Specific AoU GVS Artifacts

During the course of creating AoU callsets several large and expensive artifacts are created:

* Pilot workspace / dataset
    * For the AoU Delta callset the Variants team created an AoU 10K workspace and dataset to pilot the Hail-related
      processes we were using for the first time. At some point these processes will mature to the point where some or
      all of the contents of this workspace and dataset should be deleted. This is likely not an issue for discussion
      with Lee as it is internal to the Variants team's development process, but we should be mindful to clean up the
      parts of this that we are done using promptly.
* Production BigQuery dataset
    * The plan from Delta forward is to use the same production BigQuery dataset for all future callsets. This was also
      the plan historically as well, but in practice that didn't work out for various reasons. For Delta in particular
      the Variants team was forced to create a new dataset due to the use of drop state NONE for Hail compatibility.
      What should the Variants team do with all the various datasets for previous callsets?
      See [Topics for discussion with Lee](#topics-for-discussion-with-lee) below.
* "Prepare" tables
    * These tables are used for callset statistics only starting with Delta but were previously used for extract as well
      in pre-Delta callsets. This more limited beginning with Delta might allow the Variants team to more aggressively
      clean up this large set of tables.
* Terra workspace
    * It seems that VAT workflows may generate large amounts of data under the submissions "folder". e.g. ~10 TiB of
      data under this folder in the AoU 10K workspace (!). At the time of this writing the VAT process is not fully
      defined so this item may especially benefit from updates.
* Avro files (Delta onward)
    * These are huge, several times larger than the corresponding Hail VDS. It's not clear that there's any point to
      keeping these files around unless there was a bug in the Hail GVS import code that would require a patch and
      re-import.
* Hail VariantDataset (VDS) (Delta onward)
    * The Variants team creates a copy of the VDS and then delivers a copy to the AoU preprod datasets bucket. That copy
      of the VDS seems to stay in the delivery location for at least a few days, but it's not clear if that copy gets
      cleaned up after AoU later copies the VDS to a production bucket. The Variants team should not rely on this copy
      of the VDS being available long-term.

## Internal sign-off protocol

The Variants team currently has the following VDS internal sign-off protocol:

* Generate a VDS for the candidate callset
* Run `GvsPrepareRangesCallset` to generate a "prepare" table of VET data
* Generate callset statistics for the candidate callset using the "prepare" VET created in the preceding step
* Forward VDS and callset statistics to Lee for QA / approval

### Successful internal sign-off

Successful completion of this internal sign-off protocol currently
gates [delivery of the VDS to AoU](../vds/delivery/Delivering%20a%20VDS.md). Internal sign-off might also be used to
trigger cleanup of the "prepare" tables used in callset statistics generation (As of the Delta callset AoU uses a Hail
VDS callset workflow and does not use these "prepare" tables for extract). AoU Delta regenerate versus store costs for "
prepare" tables:

### Should successful delivery of the VDS trigger a timer for deletion of the Variants team's copy of the VDS?

Anecdotally the Variants team has been able to access the copy of the VDS delivered to AoU for at least a couple of days
after we "throw it over the wall", but we don't know if that copy of the VDS will remain accessible in the long term.

* Current Delta "prepare" tables are named like `delta_v2_no_aian_%` and `delta_v2_no_ext_aian_%`. There are multiple
  tables named with each prefix. While all the tables of a given prefix should be deleted, only the `%__VET_DATA` tables
  are large.

### Failed internal signoff

If the VDS that is delivered to Lee does not get his approval, the Variants team may need to do one or more of the
following:

* Generate new callset statistics => delete old "prepare" tables
* Create a new VDS based on the rejected VDS => delete copy of rejected VDS once new VDS is generated
* Re-import from Avro files => only if there's a bug in Hail import code?

### Prepare table storage versus regeneration

* Running
  GvsPrepareRangesCallset: [$1,803.14](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0)
* Storing "prepare" data: $878.39 / month (currently two "prepare" sets in Delta)

QUESTION: Is there a point to keeping this "prepare" data around after callset statistics have been generated? If we
needed to remove samples wouldn't we need to generate a fresh set of "prepare" tables?

### Avro file storage versus regeneration

The Avro files generated from the Delta callset onward are very large, several times the size of the final Hail VDS.
For the ~250K sample Delta callset the Avro files consume nearly 80 TiB of GCS storage while the delivered VDS is
"only" about 26 TiB.

Approximate figures for the ~250K sample Delta callset:

* Avro storage cost: $1568 / month (might be lowered if could be copied to a colder bucket).
    * `76.61 TiB  gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/submissions/c86a6e8f-71a1-4c38-9e6a-f5229520641e/GvsExtractAvroFilesForHail/efb3dbe8-13e9-4542-8b69-02237ec77ca5/call-OutputPath/avro`
* [Avro generation cost](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0):
  $3000, 12 hours runtime.

### Hail VariantDataset (VDS) storage versus regeneration

The Hail VDS generated for the Delta callset consumes about 26 TiB of space in GCS at a cost of approximately $500 /
month. Recreating the VDS from Avro files would take around 10 hours at about $100 / hour in cluster time for a total of
about $1000. Note that re-creating the VDS requires Avro files; if we have not retained the Avro files per the step
above, we would need to regenerate those as well which would add significantly to the cost.

Approximate figures for the ~250K samples Delta callset:

* VDS storage cost: ~$500 / month. Note AoU should have exact copies of the VDSes we have delivered for Delta, though
  it's not certain that these copies will remain accessible to the Variants team in the long term.
    * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/2022-10-19/dead_alleles_removed_vs_667_249047_samples/gvs_export.vds`
        * First version of the callset, includes many samples that were later removed (multiple rounds of AI/AN,
          controls)
    * ~~`gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_control_filtered_2022_12_13.vds/`~~
    * ~~`gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_filtered_2022_11_15.vds/`~~
    * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/01-04-2023-correct-GT` (LATEST)
* VDS regeneration cost: $1000 (~10 hours @ ~$100 / hour cluster cost) + $3000 to regenerate Avro files if necessary.

# Topics for discussion with Variants team

* Is there a point to keeping "prepare" tables once callset statistics have been generated?
* Are there things we can clean up in the AoU 10K workspace?
* Are there things we can clean up in the AoU 10K dataset? There are many TiBs of "prepare" tables here.
* Is there a point to keeping Avro files? We can additionally discuss this with Lee with costs 

# Topics for discussion with Lee

* Disposition toward expensive pre-Delta callsets [VS-747](https://broadworkbench.atlassian.net/browse/VS-747)
  * If we removed these we might not be able to generate exact copies 
* Should we keep the first VDS generated above? `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/2022-10-19/dead_alleles_removed_vs_667_249047_samples/gvs_export.vds`
  * Storage $500 / month, regeneration $3000 Avro + $1000 VDS
* Should we keep a copy of the most recent VDS? `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/01-04-2023-correct-GT`
  * Storage $500 / month, regeneration $3000 Avro + $1000 VDS
