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
- `use_manual_clinvar_update`: set to `true` to annotate with a manually rebuilt ClinVar database (currently the July 2025 release) instead of the October 2023 ClinVar that ships with Nirvana 3.18.1. This works with or without the 'Use reference disks' option. The database location is set by `manual_clinvar_path_prefix`. See [the ClinVar rebuild notes](../docs/vat_manual_ClinVar_updating/summary.md) for how the database was built and validated.
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

Run `GvsCreateParticipantMappingFromVDS.wdl` to build it. The workflow reads carriers straight from the VDS the VAT was built from, on a Dataproc cluster, and files each one under the VAT VID of its left-aligned representation, so carriers whose input gVCFs held a non-left-aligned form are mapped along with everyone else. Specify:

1. `vds_path`: the VDS the VAT was built from.
1. `project_id`, `dataset_name`: where to create the mapping tables.
1. `fq_vat_table`: the VAT created above, as `project.dataset.table`.
1. `fq_sample_table`: the callset's `sample_info`, as `project.dataset.table`. The workflow fails early unless the VDS holds exactly the active samples in it: none withdrawn, no controls, none it does not list, and none missing.
1. `participant_mapping_table_name`: the name of the delivered view. The workflow creates `<name>_base` (`het_ids`, `hom_ids` and `hemi_ids` per VID), `<name>_provenance` (the VIDs some of whose carriers came from a non-left-aligned representation) and the view `<name>`, with the `vid`, `person_ids` shape the Researcher Workbench expects. None of these may already exist. Its staging tables (`<name>_pairs`, `<name>_normalized_pairs` and `<name>_key_to_vid`) expire after seven days; genome-wide the pairs table is around 77 TB.

The Hail step's outputs default to a path in the workspace bucket named for the dataset and the view; set `mapping_inputs_path` to put them elsewhere. Relaunching with the same VDS and path resumes from the Hail step's checkpoint after a cluster failure. `interval` (e.g. `chr21`) restricts a test build to one contig.

The workflow's last task checks the tables against the VAT and fails if any of `mismatch_outside_provenance`, `below_inside_provenance`, `hom_mismatch_outside_provenance`, `mapped_not_in_vat`, `vat_carriers_without_row`, `provenance_not_mapped` or `vids_with_duplicate_person` is nonzero. The counts are in the `acceptance_checks` output. `above_inside_provenance` is expected to be nonzero: those are the VIDs whose non-left-aligned carriers were added to the VAT's left-aligned count.

### Delivery Steps
1. Once the VAT table is created and a TSV is exported, the AoU Researcher Workbench team should be notified of its creation and permission should be granted so that several members of the team have view permission.
    - Grant `BigQuery Data Viewer` permission to specific people's PMI-OPS accounts. This will include members of the AoU Researcher Workbench team.
    - Copy the tarred and bgzipped export of the VAT into the pre-delivery bucket.
    - Send an email out notifying the AoU Researcher Workbench team of the readiness of the VAT. Additionally, a RW Jira ticket will be made by project management to request copying the VAT to pre-prod.
    - A document describing how this information was shared (for previous callsets) is located [here](https://docs.google.com/document/d/1caqgCS1b_dDJXQT4L-tRxjOkLGDgRNkO9eac1xd9ib0/edit)
1. Copy the created mapping table to the dataset specified by the All of Us DRC. Further details of this process are included in the Google Doc linked in the previous step.
