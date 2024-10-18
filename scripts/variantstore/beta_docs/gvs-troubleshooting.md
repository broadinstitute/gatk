# Troubleshooting the Genomic Variant Store workflow

Generally, if you have started the GVS workflow and it failed after ingestion, or due to any reason that wasn't cloud failure, we recommend that you delete your BigQuery dataset and recreate it with a new name to re-run. If it fails during ingestion due to cloud issues, it can be restarted and will pick up and ingest only the remaining samples.

## Access and auth errors
1. `Failed to create dataset: Access Denied: User does not have bigquery.datasets.create permission.`
   1. See steps 5 and 6 in the [quickstart](./gvs-quickstart.md) to assign the right permissions on your google project for your Terra proxy user.
2. `BadRequestException: 400 Bucket is a requester pays bucket but no user project provided.`
   1. GVS can ingest data from a requester pays bucket by setting the optional `billing_project_id` input variable. This variable takes a string of a Google project ID to charge for the egress of the GVCFs and index files.

## Ingestion-Specific Issues"
1. GVS is running very slowly!
   1. If your GVS workflow is running very slowly compared to the example runtimes in the workspace, you may have run GVS on GVCFs that have not been reblocked. Confirm your GVCFs are reblocked.
1. My workflow failed during ingestion, can I restart it?
   1. If it fails during ingestion, yes, the GvsBeta workflow is restartable and will pick up where it left off.


## Runtime errors
1.  `Duplicate sample names error: ERROR: The input file ~{sample_names_file} contains the following duplicate entries:`
   1. The GVS requires that sample names are unique because the sample names are used to name the samples in the VCF, and VCF format requires unique sample names.
   1. After deleting or renaming the duplicate sample, you can restart the workflow without any clean up.
1. During Ingest: `Required file output '/cromwell_root/gvs_ids.csv' does not exist.`
   1. If you've attempted to run GVS more than once in the same BigQuery dataset, you may see this error. Please delete the dataset and create a new one. We recommend naming the new dataset something different than the one you deleted.
1. AssignIds failure with error message: `BigQuery error in mk operation: Not found: Dataset`
   1. This is saying that GVS was unable to find the BigQuery dataset specified in the inputs. If you haven't created a BigQuery dataset prior to running the workflow, you can follow the steps in the quickstart. If you created it and still see this error, check the naming of the dataset matches your input specified and that the google project in the inputs is correct. Lastly, confirm you have set up the correct permissions for your Terra proxy account following the instructions in the quickstart.
1. Ingest failure with error message: `raise ValueError("vcf column not in table")`
   1. if you have given an incorrect name for the vcf column or the vcf index column
   1. You can simply restart the workflow with the correct names
1. Ingest failure with error message: `Invalid resource name projects/gvs_internal; Project id: gvs_internal.`
   1. This occurs if you have given the incorrect name of the project.
   1. Restart the workflow with the correct name
1. Ingest failure with `Max id is 0. Exiting.`
   1. You will want to completely start over and delete your BQ dataset--and then re-create it. It can have the exact same name.
1. Ingest failure: There is already a list of sample names. This may need manual cleanup. Exiting.
   1. Clean up the BQ dataset manually by deleting it and recreating it fresh
   1. Make sure to keep the call caching on and run it again
1. Ingest failure with error message: `A USER ERROR has occurred: Cannot be missing required value for `___
   1. (e.g. alternate_bases.AS_RAW_MQ, RAW_MQandDP or RAW_MQ)
   1. This means that there is at least one incorrectly formatted sample in your data model. Confirm your GVCFs are reblocked. If the incorrectly formatted samples are a small portion of your callset and you wish to just ignore them, simply delete the from the data model and restart the workflow without them. There should be no issue with starting from here as none of these samples were loaded.
1. Extract failure with OSError: Is a directory. If you point your extract to a directory that doesn’t already exist, it will not be happy about this. Simply make the directory and run the workflow again.
1. Ingest failure with: `Lock table error`
   1. This means that the lock table has been created, but that the ingest has failed soon after or that perhaps during manual cleanup from another failure, some underlying data was deleted
   1. The lock table can simply be deleted -- `sample_id_assignment_lock` -- and the ingest can be kicked off again
