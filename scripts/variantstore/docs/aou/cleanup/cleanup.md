# Callset cleanup

During the normal course of creating an AoU callset several large and expensive artifacts are created:

* BigQuery dataset (new dataset or adding to an existing dataset)
* Terra workspace
* Avro files (Delta onward)
* Hail VariantDataset (VDS) (Delta onward)

The current Variants policy for AoU callsets is effectively to retain all artifacts for all callset versions forever. As
the storage costs for these artifacts can be very significant (particularly from Delta onward), the Variants team would
like to make the cost of retaining artifacts more clear so AoU can make conscious choices about what they want to keep.

## Establish an AoU callset signoff protocol

### Existing internal signoff protocol

The Variants team currently has the following internal VDS signoff protocol:

* Generate a VDS for the candidate callset
* Generate callset statistics for the candidate callset
* Forward VDS and callset statistics to Lee for QA / approval

Successful completion of this internal signoff protocol gates
the [delivery of the VDS to AoU](../vds/delivery/Delivering%20a%20VDS.md).

POLICY FOR PREPARE TABLES

### Future AoU signoff protocol

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

* Avro storage cost: $1500 / month.
* [Avro generation cost](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0):
  $3000, 12 hours runtime.

### Hail VariantDataset (VDS) storage versus regeneration

The Hail VDS generated for the Delta callset consumes about 26 TiB of space in GCS at a cost of approximately $500 /
month. Recreating the VDS from Avro files would take around 10 hours at about $100 / hour in cluster time for a total of
about $1000. Note that re-creating the VDS requires Avro files; if we have not retained the Avro files per the step
above, we would need to regenerate those as well which would add significantly to the cost.

Approximate figures for the ~250K samples Delta callset:

* VDS storage cost: $500 / month. Note AoU should have exact copies of the two VDSes generated for Delta.
* VDS regeneration cost: $1000 (~10 hours @ ~$100 / hour cluster cost) + $3000 to regenerate Avro files if necessary.
