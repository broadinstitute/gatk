# Frequently Asked Questions about the Genomic Variant Store

1. Is the resulting callset different from the WARP Joint Calling Pipeline?
   1. Yes, to enable scale, cost, and runtime efficiencies in GVS, we have optimized what data is output in the final VCFs.
2. How do I confirm the samples included in the callset?
   1. Alongside the output vcfs there is a file called `sample-name-list.txt`
2. How can I know which sharded vcf files should be merged together? Is there any file indicating the relationship between sharded vcf and chromosome?
   3. Yes, please see details on the [outputs of gvs](./gvs-outputs.md).

## Runtime
1. For running GVS, do we need to make a separate BigQuery dataset for every callset we want to create, or just once and then that can hold multiple cohorts across GVS runs?
   1. New data = new dataset in BQ. If adding data to an existing callset that hasn't been cleaned up, you could run GVS again in the same workspace and it would add the new data only and recreate the callset.


