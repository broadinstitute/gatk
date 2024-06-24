# Frequently Asked Questions about the Genomic Variant Store

1. Is the resulting callset different from the WARP Joint Calling Pipeline?
   1. Yes, to enable scale, cost, and runtime efficiencies in GVS, we have optimized what data is output in the final VCFs. Please see details on the [outputs of gvs](./gvs-outputs.md).
   2. Additionally, GVS defaults to use GATK VETS over the prior tool, GATK VQSR. See the [release notes](https://github.com/broadinstitute/gatk/blob/ah_var_store/scripts/variantstore/docs/release_notes/VETS_Release.pdf).
2. How do I confirm the samples included in the callset?
   1. Alongside the output vcfs there is a file called `sample-name-list.txt`
3. How can I know which sharded vcf files should be merged together? Is there any file indicating the relationship between sharded vcf and chromosome?
   1. Yes, please see details on the [outputs of gvs](./gvs-outputs.md).
4. There are so many variants! I thought this was filtered.
   1. The GVS 'soft filters' the data in that we tag site-level filters and allele specific filters onto the data for the user to be able to filter based on their research needs. Please see details on the [outputs of gvs](./gvs-outputs.md).

## Runtime
1. For running GVS, do we need to make a separate BigQuery dataset for every callset we want to create, or just once and then that can hold multiple cohorts across GVS runs?
   1. New data = new dataset in BQ. If adding data to an existing callset that hasn't been cleaned up, you could run GVS again in the same workspace and it would add the new data only and recreate the callset.

## VETS
1. What training resources and truth data are used for VETS model training?

The following resources are used to train the models for both SNPs and INDELS:

| Resource     | Training | Calibration | Details                                                                                                | Data Location                                                                                                 |
|--------------|----------|-------------|--------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------|
| HapMap       | TRUE     | TRUE        | This resource is a SNP callset that has been validated to a very high degree of confidence.            | `gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz`                                       |
| Omni         | TRUE     | TRUE        | This resource is a set of polymorphic SNP sites produced by the Omni genotyping array.                 | `gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz`                                    |
| 1000 Genomes | TRUE     | FALSE       | This resource is a set of high-confidence SNP sites produced by the 1000 Genomes Project.              | `gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz`                |
| Mills        | TRUE     | TRUE        | This resource is an Indel callset that has been validated to a high degree of confidence.              | `gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`             |
| Axiom        | TRUE     | FALSE       | This resource is an Indel callset based on the Affymetrix Axiom array on 1000 Genomes Project samples. | `gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz`  |
