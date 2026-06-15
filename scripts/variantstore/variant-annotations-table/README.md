# Creating the Variant Annotations Table (VAT)

The pipeline takes in a Hail Variant Dataset (VDS) or a sites only VCF, creates a queryable table in BigQuery, and outputs a bgzipped TSV file with the contents of that table.


### VAT WDLs

- [GvsCreateVATfromVDS.wdl](/scripts/variantstore/variant-annotations-table/GvsCreateVATfromVDS.wdl) creates a sites only VCF from a VDS if no input sites only VCF is specified, and then uses that and an ancestry file TSV to build the VAT.
- [GvsValidateVAT.wdl](/scripts/variantstore/variant-annotations-table/GvsValidateVAT.wdl) checks and validates the created VAT and prints a report of any failing validation.

### Using the Nirvana reference image in Terra (AoU Echo and later callsets)

The Variants team created a Cromwell reference image containing all reference files used by Nirvana 3.18.1. This is
useful to avoid having to download tens of GiBs of Nirvana references in each shard of the scattered `AnnotateVCF` task.
In order to use this reference disk, the 'Use reference disks' option in Terra must be selected as shown below:

![Terra Use reference disks](Reference%20Disk%20Terra%20Opt%20In.png)

### Run GvsCreateVATfromVDS

- **Note:** in order for this workflow to run successfully the 'Use reference disks' option must be selected in Terra workflow
configuration. If this option is not selected the `AnnotateVCF` tasks will fail because the Nirvana reference files will not be found. If you forget and it fails, re-run it with call-caching on and all the same inputs, and it will resume at the right point.
- **Note:** due to an [open issue with GCP Batch](https://partnerissuetracker.corp.google.com/issues/449751210) it is expected that some shards of `GenerateVepAndLofteeAnnotations` will fail
with 50002 errors that are caused by this task consuming all the memory on the VM, leaving the Batch agent unable to checkin with the Batch service.
If and when the workflow fails in this manner, adjust the value in the `memory` runtime attribute of the `GenerateVepAndLofteeAnnotations` task and rerun `GvsCreateVATfromVDS` with call caching enabled.
Do not change any task inputs or call caching will break, just edit the value of the runtime attribute directly.
For the Foxtrot dataset, all but around 40 shards were able to complete with the baseline 8 GiB of memory; modifying this task to use 32 GiB of memory and rerunning the workflow enabled those remaining shards to succeed.
- The `ancestry_file` input is the GCS path of the TSV file that maps samples (by `sample_name`) to subpopulations.
- You will want to run this workflow with the same `dataset_name`, `project_id`, and `filter_set_name` as `GvsCreateVds.wdl`.
- For `output_path` use a unique GCS path with a trailing slash (probably in the workspace bucket). This will be used to store the intermediate files for the pipeline.
- The `vds_path` input is the same value that was set for `vds_destination_path` in `GvsCreateVds.wdl`.
- This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.

Optional inputs of note:

- `split_intervals_scatter_count`: If you want to override the step function that decides this based on the number of samples (for Delta we used a scatter of 500, and for Echo a scatter of 1000).
- `vat_version`: if you are creating multiple VATs for one callset, you can distinguish between them (and not overwrite others) by passing in increasing numbers
- If you are debugging a Hail-related issue, you may want to set `leave_hail_cluster_running_at_end` to `true` and refer to [the suggestions for debugging issues with Hail](../docs/aou/HAIL_DEBUGGING.md). 

There are several temporary tables that are created in addition to the main VAT table. The Genes, VT (variant transcripts), and intermediate `_w_dups` VAT tables all have a time to live of 24 hours. The VEP/LOFTEE raw and cooked tables have a time to live of 3 days. The final VAT table is (re)created fresh each time so that there is no risk of duplicates.

Variants may be filtered out of the VAT (that were in the VDS) for the following reasons:

- they are hard-filtered out based on the initial soft filtering from the GVS extract (site- and GT-level filtering)
- they have excess alternate alleles, currently this filters out sites with >= 100 alternate alleles
- they are spanning deletions
- they are duplicate variants; they are tracked via the `GvsCreateVATfromVDS` workflow's scattered `RemoveDuplicatesFromSitesOnlyVCF` task and then merged into one file by the `MergeTsvs` task


### Run GvsValidateVAT

This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option. The `project_id` and `dataset_name` are the same as those used for `GvsCreateVATfromVDS`, and `vat_table_name` is `filter_set_name` + "_vat" (+ "_v" + `vat_version`, if used).

### VID to Participant ID Mapping Table

Once the VAT has been created, you will need to create a database table mapping the VIDs (Variant IDs) from that table to all the participants in the dataset that share that VID. This table is used by the AoU Researcher Workbench, and will need to be copied over to a location specified by them.

1. First run `GvsCreateParticipantMappingTable.wdl` to create this participant ID mapping table and most of its data.
   Specify the `project_id`, `dataset` and `vat_table_name` for the VAT created above. Also specify the `participant_mapping_table_name` that
   will be created to hold the VID to participant mapping information.
1. Next run `GvsMapUnmappedVIDs.wdl` to recover participant ID mappings for "unmapped VIDs" (a.k.a "pseudo VIDs"). These unmapped VIDs have
   VAT table entries but no corresponding `alt_allele` entries and correspond to data that appeared in
   input VCFs only in non-left aligned representations. During the process of creating the VAT these variants were left-aligned and thus
   disconnected from their `alt_allele` representations. Specify the following parameters:
   1. `project_id`, `dataset`, `vat_table_name`, `participant_mapping_table_name`: use the same values as in the step above for `GvsCreateParticipantMappingTable.wdl`.
   1. `sites_only_vcf`: the GCS path to the sites-only VCF file that was generated in the process of creating the VAT. This corresponds to the `output_file_path` output of the `CopySitesOnlyVcf` task in `GvsCreateVATfromVDS.wdl`.
   1. `unmapped_vid_mapping_table_name`: the name to use for a table that will hold the mapping information from input
       position/ref/alt to left-aligned position/ref/alt. This should be a new table.
1. Finally run `GvsMapDroppedDuplicateVIDs.wdl` to recover all participant ID mappings for VIDs which had multiple variant synonyms with AC != 0. The logic in `GvsCreateVATfromVDS` currently preserves a left-aligned version of the
   synonym with the highest AC, but before running `GvsMapDroppedDuplicateVIDs.wdl` the actual participant mapping will only contain samples whose input synonym was left-aligned, which do not necessarily correspond to the synonym with highest AC.
   Specify the following parameters:
   1. `project_id`, `dataset`, `vat_table_name`, `participant_mapping_table_name`: use the same values as in the step above for `GvsCreateParticipantMappingTable.wdl`.
   1. `sites_only_vcf`: the GCS path to the sites-only VCF file that was generated in the process of creating the VAT. This corresponds to the `output_file_path` output of the `CopySitesOnlyVcf` task in `GvsCreateVATfromVDS.wdl`.
       This should be the same value as specified for the `GvsMapUnmappedVIDs.wdl` step above.
   1. `filtered_synonyms`: the GCS path to the file of variant synonyms that were filtered as duplicates. This corresponds to the value of the `output_file` output of the `MergeDroppedSynonyms` task in `GvsCreateVATfromVDS.wdl`.
   1. `duplicate_mapping_table_name`: the name to use for a table that will hold the mapping information from input
      position/ref/alt to left-aligned position/ref/alt. This should be a new table.

### Delivery Steps
1. Once the VAT table is created and a TSV is exported, the AoU Researcher Workbench team should be notified of its creation and permission should be granted so that several members of the team have view permission.
    - Grant `BigQuery Data Viewer` permission to specific people's PMI-OPS accounts. This will include members of the AoU Researcher Workbench team.
    - Copy the tarred and bgzipped export of the VAT into the pre-delivery bucket.
    - Send an email out notifying the AoU Researcher Workbench team of the readiness of the VAT. Additionally, a RW Jira ticket will be made by project management to request copying the VAT to pre-prod.
    - A document describing how this information was shared (for previous callsets) is located [here](https://docs.google.com/document/d/1caqgCS1b_dDJXQT4L-tRxjOkLGDgRNkO9eac1xd9ib0/edit)
1. Copy the created mapping table to the dataset specified by the All of Us DRC. Further details of this process are included in the Google Doc linked in the previous step.
