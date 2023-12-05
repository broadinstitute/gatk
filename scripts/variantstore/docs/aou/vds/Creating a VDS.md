# Creating a VDS for AoU

**NOTE:** This doc is very much a work in process and is serving as a way to keep track of Hail commands that were used for Delta after the initial VDS delivery.
This example walks through the procedure we followed for creating a VDS.

## Create avro files for VDS extract
Run GvsExtractAvroFilesForHail to create a directory of avro files that will then be used to create a VDS

## Configuring the Terra notebook cluster

Please follow the instructions
in [AoU Delta VDS Cluster Configuration](cluster/AoU%20Delta%20VDS%20Cluster%20Configuration.md) to set up a Delta-scale
cluster.

## Using the created notebook, run the python script to create a VDS from the avro files created by the GvsExtractAvroFilesForHail WDL
These avro files will be passed in using the --avro-path parameter
The vds-path parameter tells the WDL where to write the final VDS to
There is quite a bit of temporary data created and while it will be cleaned up later, we require that you give the python script a location.
NOTE: in the WDL created to run this script, the temp directory can be easily defaulted to a date/random dir in the workspace directory

```python ./hail_gvs_import.py \
--avro-path gs://fc-<bucket id>/submissions/<workflow id>/GvsExtractAvroFilesForHail/<submission id>/call-OutputPath/avro/ \
--vds-path gs://fc-<bucket id>/<name of the new vds>/ \
--temp-path gs://fc-<bucket id>/<name of the temp path>/
```

## Validate the VDS to ensure that it is ready to be shared

Copy the [VDS Validation python script](../../../wdl/extract/vds_validation.py) to the notebook environment.
Run it with the following arguments:

`--vds-path`: the GCS path to the newly-created VDS

`--temp-path`: a GCS path for temp files to live


An example run (on the Anvil 3k) is illustrated below

```
python3 vds_validation.py --vds-path 'gs://fc-3ecd704c-4d2c-4360-bae1-093e214abce2/vds/vds_GQ0_max_ref_1k' --temp-path 'gs://fc-62891290-7fa3-434d-a39e-c64eeee4db8d/temp'
```