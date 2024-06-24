# Troubleshooting the Genomic Variant Store workflow

## Access and auth errors
1. `Failed to create dataset: Access Denied: User does not have bigquery.datasets.create permission.`
   1. See steps 5 and 6 in the [quickstart](./gvs-quickstart.md) to assign the right permissions on your google project for your Terra proxy user.
2. `BadRequestException: 400 Bucket is a requester pays bucket but no user project provided.`
   1. GVS can ingest data from a requester pays bucket by setting the optional `billing_project_id` input variable. This variable takes a string of a Google project ID to charge for the egress of the GVCFs and index files.                                                                                                               | String  |


## Runtime errors
1. My workflow failed during ingestion, can I restart it?
   1. If it fails during ingestion, yes, the GvsBeta workflow is restartable and will pick up where it left off.
2. Duplicate sample names error: `ERROR: The input file ~{sample_names_file} contains the following duplicate entries:`
   1. The GVS requires that sample names are unique because the sample names are used to name the samples in the VCF, and VCF format requires unique sample names. 
   2. After deleting or renaming the duplicate sample, you can restart the workflow without any clean up.
