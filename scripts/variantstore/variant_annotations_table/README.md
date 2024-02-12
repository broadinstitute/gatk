# Creating the Variant Annotations Table

The pipeline takes in a Hail Variant Dataset (VDS), creates a queryable table in BigQuery, and outputs a bgzipped TSV file containing the contents of that table. That TSV file is the sole AoU deliverable.


### VAT WDLs

- [GvsCreateVATfromVDS.wdl](/scripts/variantstore/wdl/GvsCreateVATfromVDS.wdl) creates a sites only VCF  from a VDS and then uses that and an ancestry file TSV to build the variant annotations table.
- [GvsValidateVAT.wdl](/scripts/variantstore/variant_annotations_table/GvsValidateVAT.wdl) checks and validates the created VAT and prints a report of any failing validation.

### Run GvsCreateVATfromVDS

- **Note:** in order for this workflow to run successfully the 'Use reference disks' option must be selected in Terra workflow
configuration. If this option is not selected the `AnnotateVCF` tasks will refuse to run. If you forget and it fails, re-run it with call-caching on and all the same inputs, and it will resume at the right point.
- The `ancestry_file` input is the GCS path of the TSV file that maps samples (by `sample_name`) to subpopulations.
- You will want to run this workflow with the same `dataset_name`, `project_id`, and `filter_set_name` as `GvsCreateVds.wdl`.
- For `hail_generate_sites_only_script_path` refer back to the run of `GvsExtractAvroFilesForHail`. The `hail_create_vat_inputs_script` output of the `GenerateHailScripts` task is GCS the path to use.
- For `output_path` use a unique GCS path with a trailing slash (probably in the workspace bucket). This will be used to store the intermediate files for the pipeline.
- The `vds_path` input is the same value that was set for `vds_destination_path` in `GvsCreateVds.wdl`.
- This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option.

Optional input of note:

- `vat_version`: if you are creating multiple VATs for one callset, you can distinguish between them (and not overwrite others) by passing in increasing numbers
- If you are debugging a hail-related issue, you may want to set `leave_hail_cluster_running_at_end` to `true` and refer to [the guidelines for debugging issues with Hail](../docs/aou/HAIL_DEBUGGING.md). 

There are two temporary tables that are created in addition to the main VAT table: the Genes and VT tables. They have a time to live of 24 hours.  The VAT table is created by that query fresh each time so that there is no risk of duplicates.

Variants may be filtered out of the VAT (that were in the VDS) for the following reasons:

- they are hard-filtered out based on the initial soft filtering from the GVS extract (site- and GT-level filtering)
- they have excess alternate alleles, currently that cut off is 50 alternate alleles
- they are spanning deletions
- they are duplicate variants; they are tracked via the `GvsCreateVATfromVDS` workflow's scattered `RemoveDuplicatesFromSites` task and then merged into one file by the `MergeTsvs` task

### Remove duplicates from the VAT table.

There is a bug in the way we build the VAT table from the shards. If a variant (e.g. for an insertion or deletion) spans between the end of one shard and the beginning of the next shard, it will be entered TWICE into the resultant vat table as there will be duplicate entries in each of the shards. To address this issue, we run a manual SQL statement in BigQuery in order to remove these duplicates.
**Note:** The purpose of this query is to create a NEW copy of the VAT table with the duplicates removed. So, when you run it, have BigQuery write the output to a new table (More -> Query Settings -> Destination )

>     select * except(row_number) from (
>      SELECT
>       *,
>       row_number()
>       over (partition by vid, transcript)
>       row_number
>       FROM
>        `<project_id>.<dataset>.<vat_table_name>`    
>       )
>     where row_number = 1;

### Run GvsValidateVAT

This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option. The `project_id` and `dataset_name` are the same as those used for `GvsCreateVATfromVDS`, and `vat_table_name` is `filter_set_name` + "_vat" (+ "_v" + `vat_version`, if used).

### Using the Nirvana reference image in Terra (AoU Echo and later callsets)

The Variants team created a Cromwell reference image containing all reference files used by Nirvana 3.18.1. This is
useful to avoid having to download tens of GiBs of Nirvana references in each shard of the scattered `AnnotateVCF` task.
In order to use this reference disk, the 'Use reference disk' option in Terra must be selected as shown below:

![Terra Use reference disks](Reference%20Disk%20Terra%20Opt%20In.png)

**Note:** this Nirvana reference image was not available in time for it to be used to create the Delta VAT. The
`GvsCreateVATfromVDS` WDL has been updated to take advantage of references loaded from an attached
reference image and should be ready to use starting with the AoU Echo callset.
