# AoU callset cleanup

During the normal course of creating an AoU callset several large and expensive artifacts are created:

* BigQuery dataset
    * The plan from Delta forward is to use the same BigQuery dataset for all future callsets. This was also the
      plan historically as well, but in practice that didn't work out for various reasons. For Delta in particular we
      were compelled to create a new dataset to support the use of drop state NONE for Hail compatibility.
    * "Prepare" tables for callset statistics starting with Delta (also used for extract in pre-Delta callsets).
* Terra workspace
    * It seems that VAT workflows may generate large amounts of data under the submissions "folder". e.g. ~10 TiB of
      data under this folder in the AoU 10K workspace (!)
* Avro files (Delta onward)
    * These are huge, several times larger than the corresponding Hail VDS which itself is tens of TiB.
* Hail VariantDataset (VDS) (Delta onward)
    * ~26 TiB for Delta
* Pilot workspaces
    * For Delta we created an AoU 10K workspace to pilot the various Hail-related processes we were using for the first
      time. We should delete this at some point but it's not currently clear when.

The current Variants policy for AoU callsets is effectively to retain all versions of all artifacts forever. As the
storage costs for these artifacts can be significant (particularly from Delta onward), the Variants team would like to
make the cost of retaining artifacts more clear so conscious choices can be made about what to keep and what to delete.

As a general rule, any of these large artifacts that have clearly become obsolete (e.g. VDS versions with known issues
that have been superseded by corrected versions, "prepare" tables for old callset statistics runs etc., entire
lower-scale workspaces like AoU 10K) should be deleted ASAP. If it's not clear whether an artifact should be cleaned up
or not, figure out the monthly cost to preserve the artifact (e.g. GCS or BigQuery storage costs) as well as the cost to
regenerate the artifact. Reach out to Lee with these numbers for his verdict on whether to keep or delete.

## Signoff protocol

The Variants team currently has the following VDS signoff protocol:

* Generate a VDS for the candidate callset
* Run `GvsPrepareRangesCallset` to generate a "prepare" table for VET data
* Generate callset statistics for the candidate callset using the "prepare" VET created in the preceding step
* Forward VDS and callset statistics to Lee for QA / approval

### Successful internal signoff

Successful completion of this internal signoff protocol currently
gates [delivery of the VDS to AoU](../vds/delivery/Delivering%20a%20VDS.md). Internal signoff might also be used to
trigger cleanup of the "prepare" tables used in callset statistics generation (As of the Delta callset AoU uses a Hail
VDS callset workflow and does not use these "prepare" tables for extract). AoU Delta regenerate versus store costs for
"prepare" tables:

* Running
  GvsPrepareRangesCallset: [$1,803.14](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0)
* Storing "prepare" data: $878.39 / month

PROPOSED: Successful completion of the internal signoff protocol should trigger a timer for deletion of the "prepare"
tables used to generate callset statistics.

QUESTION: Should successful delivery of the VDS trigger a timer for deletion of the Variants team's copy of the VDS?
Anecdotally the Variants team has been able to access the copy of the VDS delivered to AoU for at least a couple of days
after we "throw it over the wall", but we don't know if that copy of the VDS would remain accessible long term.

PROPOSED: Regeneration of callset statistics for any reason (e.g. removal of samples) should trigger the deletion of the
previous versions of "prepare" tables.

Current Delta "prepare" tables are named like `delta_v2_no_aian_%` and `delta_v2_no_ext_aian_%`. There are multiple
tables named with each prefix. While all the tables of a given prefix should be deleted, only the `%__VET_DATA` tables
are large.

* `delta_v2_no_aian__VET_DATA` $878.39 / month
* `delta_v2_no_ext_aian__VET_DATA` $878.28 / month

### Failed internal signoff

If the VDS that is delivered to Lee does not get his approval, the Variants team may need to do one or more of the
following:

* Generate new callset statistics
* Create a new VDS based on the rejected VDS
* Re-import from Avro files

### Avro file storage versus regeneration

The Avro files generated from the Delta callset onward are very large, several times the size of the final Hail VDS.
For the ~250K sample Delta callset the Avro files consume nearly 80 TiB of GCS storage while the delivered VDS is
"only" about 26 TiB.

Approximate figures for the ~250K sample Delta callset:

* Avro storage cost: $1568 / month.
    * `76.61 TiB  gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/submissions/c86a6e8f-71a1-4c38-9e6a-f5229520641e/GvsExtractAvroFilesForHail/efb3dbe8-13e9-4542-8b69-02237ec77ca5/call-OutputPath/avro`
* [Avro generation cost](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0):
  $3000, 12 hours runtime.

### Hail VariantDataset (VDS) storage versus regeneration

The Hail VDS generated for the Delta callset consumes about 26 TiB of space in GCS at a cost of approximately $500 /
month. Recreating the VDS from Avro files would take around 10 hours at about $100 / hour in cluster time for a total of
about $1000. Note that re-creating the VDS requires Avro files; if we have not retained the Avro files per the step
above, we would need to regenerate those as well which would add significantly to the cost.

Approximate figures for the ~250K samples Delta callset:

* VDS storage cost: ~$500 / month. Note AoU should have exact copies of the VDSes we have delivered for Delta.
    * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/submissions/c86a6e8f-71a1-4c38-9e6a-f5229520641e/GvsExtractAvroFilesForHail/efb3dbe8-13e9-4542-8b69-02237ec77ca5/call-OutputPath/2022-10-19-6497f023/dead_alleles_removed_vs_667_249047_samples`
      ==> `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/2022-10-19/dead_alleles_removed_vs_667_249047_samples/gvs_export.vds`
        * First version of the callset, includes many samples that were later removed (multiple rounds of AI/AN,
          controls)
    * ~~`gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_control_filtered_2022_12_13.vds/`~~
    * ~~`gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_filtered_2022_11_15.vds/`~~
    * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/01-04-2023-correct-GT` (LATEST)
* VDS regeneration cost: $1000 (~10 hours @ ~$100 / hour cluster cost) + $3000 to regenerate Avro files if necessary.

## New BigQuery Dataset required

If it turns out that we need to create a new BigQuery dataset for the callset we're working on, the Variants team should
ask Lee what we should do with the old dataset(s). Present Lee with the costs of retaining the old datasets to help him
make his decision.

## New set of "prepare" tables required

If a new set of prepare tables is required it should be safe to delete the previous set.

## New Avro files generated

If new Avro files are generated (fully or in part) the superseded Avro files should be deleted.

## New VDS generated

If a new version of a VDS is being created for a callset, consider deleting the previous version.

# Topics for discussion with Lee

* Disposition toward pre-Delta callsets [VS-747](https://broadworkbench.atlassian.net/browse/VS-747)
    * These are sizeable and expensive.
