# Creating the Variant Annotations Table

### The VAT pipeline is a python script and a set of WDLs

- hail_create_vat_inputs.py is a script to be run in a notebook terminal to generate a sites-only VCF.
- [GvsCreateVATfromVDS.wdl](/scripts/variantstore/wdl/GvsCreateVATfromVDS.wdl) creates the variant annotations table from a sites only VCF and an acestry file TSV.
- [GvsValidateVAT.wdl](/scripts/variantstore/variant_annotations_table/GvsValidateVAT.wdl) checks and validates the created VAT and prints a report of any failing validation.

The pipeline takes in a Hail Variant Dataset (VDS), creates a queryable table in BigQuery, and outputs a bgzipped TSV file containing the contents of that table. That TSV file is the sole AoU deliverable.

### Generate a Sites-Only VCF

You will want to set up a [notebook in a similar configuration](../docs/aou/vds/cluster/AoU%20Delta%20VDS%20Cluster%20Configuration.md) as you did for creating the VAT (minus the custom wheel) and then copy two python scripts for running locally in the notebook. The two scripts are:

* create_vat_inputs.py which you can get via `curl https://raw.githubusercontent.com/broadinstitute/gatk/ah_var_store/scripts/variantstore/wdl/extract/create_vat_inputs.py -o create_vat_inputs.py`
* the `hail_create_vat_inputs_script` output from the task `GenerateHailScripts`from the run of`GvsExtractAvroFilesForHail` workflow that was used to generate the VDS.

Once you have copied them, run `python3 hail_create_vat_inputs.py` with the following input(s)

* `--ancestry_file` GCS pointer to the TSV file that maps `sample_name`s to sub-populations
* `--vds_path` GCS pointer to the top-level direcotory of the VDS
* `--sites_only_vcf` (optional, should be pre-populated by the `GvsExtractAvroFilesForHail` workflow) the location to write the sites-only VCF to; save this for the `GvsCreateVATfromVDS` workflow

When the script is done running, it will spit out the path to where it has written the sites-only VCF.

### Run GvsCreateVATfromVDS

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
