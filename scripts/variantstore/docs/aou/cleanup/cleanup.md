# Callset cleanup

During the normal course of creating an AoU callset several large and expensive artifacts are created:

* BigQuery dataset (new dataset or adding to an existing dataset)
* Terra workspace
* Avro files (Delta onward)
* Hail VariantDataset (VDS) (Delta onward)

The current Variants policy for AoU callsets is effectively to retain all artifacts for all callset versions forever. As
the storage costs for these artifacts can be very significant (particularly from Delta onward), the Variants team would
like to make the cost of retaining artifacts more clear so AoU can make conscious choices about what they want to keep.


## Internal signoff protocol

The Variants team currently has the following internal VDS signoff protocol:

* Generate a VDS for the candidate callset
* Generate callset statistics for the candidate callset
* Forward VDS and callset statistics to Lee for QA / approval

Successful completion of this internal signoff protocol currently
gates [delivery of the VDS to AoU](../vds/delivery/Delivering%20a%20VDS.md). Internal signoff might also be used to
trigger cleanup of the "prepare" tables used in callset statistics generation. Delta regenerate versus store costs for
"prepare" tables:

* Running GvsPrepareRangesCallset: [$1,803.14](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0)
* Storing "prepare" data: $878.39 / month

PROPOSED: Successful completion of the internal signoff protocol should trigger a timer for deletion of the "prepare"
tables used to generate callset statistics.

PROPOSED: Regeneration of callset statistics for any reason (e.g. removal of samples) should trigger the deletion of the
previous versions of "prepare" tables.

Current Delta "prepare" tables are named like `delta_v2_no_aian_%` and `delta_v2_no_ext_aian_%`. There are multiple
tables named with each prefix. While all the tables of a given prefix should be deleted, only the `%__VET_DATA` tables
are large.

* `delta_v2_no_aian__VET_DATA` $878.39 / month
* `delta_v2_no_ext_aian__VET_DATA` $878.28 / month

Generating callset statistics

## AoU signoff protocol

Analogous to the existing internal signoff protocol, the Variants team would like to have an external AoU callset
signoff protocol. The successful completion of this protocol could gate removal of intermediate data such as Avro files
and pre-delivery copies of the VDS.

There are tradeoffs in cost and time between retaining intermediate callset data versus regenerating it if required,
details of which are discussed below. Ultimately AoU is the party paying for these cloud resources so AoU should
decide what they would like the cleanup policies to be.

### Avro file storage versus regeneration

The Avro files generated from the Delta callset onward are very large, significantly larger than the final Hail VDS.
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

* VDS storage cost: $500 / month. Note AoU should have exact copies of the two VDSes generated for Delta.
  * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/submissions/c86a6e8f-71a1-4c38-9e6a-f5229520641e/GvsExtractAvroFilesForHail/efb3dbe8-13e9-4542-8b69-02237ec77ca5/call-OutputPath/2022-10-19-6497f023/dead_alleles_removed_vs_667_249047_samples`
  * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_control_filtered_2022_12_13.vds/`
  * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/delta_ai_an_filtered_2022_11_15.vds/`
  * `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/01-04-2023-correct-GT` (LATEST)
* VDS regeneration cost: $1000 (~10 hours @ ~$100 / hour cluster cost) + $3000 to regenerate Avro files if necessary.
