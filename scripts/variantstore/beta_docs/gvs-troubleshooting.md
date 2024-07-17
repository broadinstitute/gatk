# Troubleshooting the Genomic Variant Store workflow

Generally, if you have started the GVS workflow and it failed after ingestion, or due to any reason that wasn't cloud failure, we recommend that you delete your BigQuery dataset and recreate it with a new name to re-run. If it fails during ingestion due to cloud issues, it can be restarted and will pick up and ingest only the remaining samples.

## Access and auth errors
1. `Failed to create dataset: Access Denied: User does not have bigquery.datasets.create permission.`
   1. See steps 5 and 6 in the [quickstart](./gvs-quickstart.md) to assign the right permissions on your google project for your Terra proxy user.
2. `BadRequestException: 400 Bucket is a requester pays bucket but no user project provided.`
   1. GVS can ingest data from a requester pays bucket by setting the optional `billing_project_id` input variable. This variable takes a string of a Google project ID to charge for the egress of the GVCFs and index files.

## Runtime errors
1. My workflow failed during ingestion, can I restart it?
   1. If it fails during ingestion, yes, the GvsBeta workflow is restartable and will pick up where it left off.
2. Duplicate sample names error: `ERROR: The input file ~{sample_names_file} contains the following duplicate entries:`
   1. The GVS requires that sample names are unique because the sample names are used to name the samples in the VCF, and VCF format requires unique sample names. 
   2. After deleting or renaming the duplicate sample, you can restart the workflow without any clean up.
3. `BulkIngestGenomes/GvsBulkIngestGenomes/hash/call-ImportGenomes/GvsImportGenomes/hash/call-GetUningestedSampleIds/gvs_ids.csv Required file output '/cromwell_root/gvs_ids.csv' does not exist.`
   1. If you've attempted to run GVS more than once in the same BigQuery dataset, you may see this error. Please delete the dataset and create a new one. We recommend naming the new dataset something different than the one you deleted.
4. AssignIds failure with error message: `BigQuery error in mk operation: Not found: Dataset`
   1. This is saying that GVS was unable to find the BigQuery dataset specified in the inputs. If you haven't created a BigQuery dataset prior to running the workflow, you can follow the steps in [the quickstart](./gvs-quickstart.md). If you created it and still see this error, check the naming of the dataset matches your input specified and that the google project in the inputs is correct. Lastly, confirm you have set up the correct permissions for your Terra proxy account following the instructions in the quickstart. 

## Reblocking
1. `htsjdk.tribble.TribbleException$MalformedFeatureFile: Unable to parse header with error: Your input file has a malformed header: We never saw the required CHROM header line (starting with one #) for the input VCF file, for input source: file:///cromwell_root/v1_[uuid]`
   1. If you are running ReblockGVCF from a TDR snapshot, you will see this error if you did not check the “Convert DRS URLs to Google Cloud Storage Paths (gs://)" box before exporting the snapshot.
2. GVS is running very slowly!
   1. If your GVS workflow is running very slowly compared to the example runtimes in the workspace, you may have run GVS on GVCFs that have not been reblocked. Confirm your GVCFs are reblocked.