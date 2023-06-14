# Creating the Variant Annotations Table

### The VAT pipeline is a python script and a set of WDLs

- [hail_create_vat_inputs.py](https://raw.githubusercontent.com/broadinstitute/gatk/ah_var_store/scripts/variantstore/wdl/extract/hail_create_vat_inputs.py) is a script to be run in a notebook terminal to generate a sites-only VCF.
- [GvsCreateVATfromVDS.wdl](/scripts/variantstore/wdl/GvsCreateVATfromVDS.wdl) creates the variant annotations table from a sites only VCF and an ancestry file TSV.
- [GvsValidateVAT.wdl](/scripts/variantstore/variant_annotations_table/GvsValidateVAT.wdl) checks and validates the created VAT and prints a report of any failing validation.

The pipeline takes in a Hail Variant Dataset (VDS), creates a queryable table in BigQuery, and outputs a bgzipped TSV file containing the contents of that table. That TSV file is the sole AoU deliverable.

### Generate a Sites-Only VCF

You will want to set up a [notebook in a similar configuration](../docs/aou/vds/cluster/AoU%20Delta%20VDS%20Cluster%20Configuration.md) as you did for creating the VAT (minus the custom wheel) and then copy two python scripts for running locally in the notebook. The two scripts are:

* `create_vat_inputs.py` which you can get via `curl -O https://raw.githubusercontent.com/broadinstitute/gatk/ah_var_store/scripts/variantstore/wdl/extract/create_vat_inputs.py`
* The `hail_create_vat_inputs.py` output from the task `GenerateHailScripts` from the run of `GvsExtractAvroFilesForHail` workflow that was used to generate the VDS.

Once you have copied them, run `python3 hail_create_vat_inputs.py` with the following input(s):

* `--ancestry_input_path` Required, GCS path of the TSV file that maps `sample_name`s to subpopulations.
* `--vds_input_path` Optional, should be pre-populated by the `GvsExtractAvroFilesForHail` workflow. GCS path of the top-level directory of the VDS.
* `--sites_only_output_path` Optional, should be pre-populated by the `GvsExtractAvroFilesForHail` workflow. GCS path to write the sites-only VCF to; save this for the `GvsCreateVATfromVDS` workflow.

When the script is done running it will print out the path to where it has written the sites-only VCF.

### Run GvsCreateVATfromVDS

**Note:** in order for this workflow to run successfully the 'Use reference disks' option must be selected in Terra workflow
configuration. If this option is not selected the `AnnotateVCF` tasks will refuse to run.

Most of the inputs are specific to where the VAT will live or will be the same for the VDS creation, like the `project_id`, `dataset_name`, `filter_set_name`, and `output_path`. This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option. For the specific data being put in the VAT, two inputs need to be copied into a GCP bucket that this pipeline will have access to and passed to the WDL:

- `input_sites_only_vcf`: The GCS path to a sites-only VCF that we affectionately call it a Decorated Sites Only VCF because it also includes subpopulation level calculations for AC AN, AF and SC. This is the output of the `hail_create_vat_inputs.py` script.
- `ancestry_file`: The same file that was passed via `--ancestry_file` to `hail_create_vat_inputs.py`.
- `output_path`: The GCS path to store the intermediate files for the pipeline; should have a trailing slash.

Other input of note:

- `vat_version`: if you are creating multiple VATs for one callset, you can distinguish between them (and not overwrite others) by passing in increasing numbers

There are two temporary tables that are created in addition to the main VAT table: the Genes and VT tables. They have a time to live of 24 hours.  The VAT table is created by that query fresh each time so that there is no risk of duplicates.

Variants may be filtered out of the VAT (that were in the VDS) for the following reasons:

- they are hard-filtered out based on the initial soft filtering from the GVS extract (site- and GT-level filtering)
- they have excess alternate alleles, currently that cut off is 50 alternate alleles
- they are spanning deletions
- they are duplicate variants; they are tracked via the `GvsCreateVATfromVDS` workflow's scattered `RemoveDuplicatesFromSites` task and then merged into one file by the `MergeTsvs` task

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
