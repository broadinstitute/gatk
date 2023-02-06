# Creating the Variant Annotations Table

### The VAT pipeline is a set of WDLs

- [GvsCreateVATfromVDS.wdl](/scripts/variantstore/wdl/GvsCreateVATfromVDS.wdl) creates the variant annotations table from a sites only VCF and an acestry file TSV.
- [GvsValidateVAT.wdl](/scripts/variantstore/variant_annotations_table/GvsValidateVAT.wdl) checks and validates the created VAT and prints a report of any failing validation.

The pipeline takes in a Hail Variant Dataset (VDS), creates a queryable table in BigQuery, and outputs a bgzipped TSV file containing the contents of that table. That TSV file is the sole AoU deliverable.

### Run GvsCreateVATfromVDS

Most of the inputs are specific to where the VAT will live or will be the same for the VDS creation, like the `project_id`, `dataset_name`, `filter_set_name`, and `output_path`. This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option. For the specific data being put in the VAT, two inputs need to be copied into a GCP bucket that this pipeline will have access to and passed to the WDL:

- `input_sites_only_vcf`: A sites-only VCF that we affectionately call it a Decorated Sites Only VCF because it also includes subpopulation level calculations for AC AN, AF and SC. This is the `sites_only_vcf_output_path` output of the `GvsExtractAvroFilesForHail` workflow's `GenerateHailScripts` task for this callset.
- `ancestry_file`: An ancestry file from the ancestry pipeline which will be used to calculate AC, AN and AF for all subpopulations.

Other input of note:

- `vat_version`: if you are creating multiple VATs for one callset, you can distinguish between them (and not overwrite others) by passing in increasing numbers

### Run GvsValidateVAT

This workflow does not use the Terra Data Entity Model to run, so be sure to select the `Run workflow with inputs defined by file paths` workflow submission option. The `project_id` and `dataset_name` are the same as those used for `GvsCreateVATfromVDS`, and `vat_table_name` is `filter_set_name` + "_vat" (+ "_v" + `vat_version`, if used).

### Notes

There are two temporary tables that are created in addition to the main VAT table: the Genes and VT tables. They have a time to live of 24 hours.
The VAT table is created by that query fresh each time so that there is no risk of duplicates.

To check that all of the shards have successfully made it past the first of the subworkflow steps (the most complicated and likely to fail) make sure that the full number of expected files are here: `gsutil ls  [output_path]/annotations/  | wc -l`

And then once they have been transformed by the Python script and are ready to be loaded into BQ they will be here:

```
gsutil ls  [output_path]/genes/  | wc -l
gsutil ls  [output_path]/vt/  | wc -l
```

These numbers are cumulative. Also the names of these json files are retained from the original shard names so as to not cause collisions. If you run the same shards through the VAT twice, the second runs should overwrite the first and the total number of jsons should not change.
Once the shards have make it into the /genes/ and /vt/ directories, the majority of the expense and transformations needed for that shard are complete.
They are ready to be loaded into BQ. You will notice that past this step, all there is to do is create the BQ tables, load the BQ tables, run a join query and then the remaining steps are all validations or an export into TSV.

Sometimes shards will fail with a 503 from Google. The shards that are affected will need to be collected and put into a new File of File names and re-run.
In theory you could actually just re-run the pipeline entirely and since the 503s seem to be completely random and intermittent, the same shards likely wont fail twice and you’d get everything you need, but at low efficiency and high expense.
To grab any files not in the bucket, but in the file of file names:

```
gsutil ls [output_path]/genes/ | awk '{print substr($0,90)}' RS='.json.gz\n'  > successful_vcf_shards.txt
gsutil cat [inputFileofFileNames] | awk '{print substr($0,41)}' RS='.vcf.gz\n' > all_vcf_shards.txt
comm -3 all_vcf_shards.txt successful_vcf_shards.txt > diff_vcf_shards.txt
awk '{print g[bucket GVS output VCF shards are in]" $0 ".vcf.gz"}' diff_vcf_shards.txt > fofn_missing_shards.txt
```

_(dont forget to do this for the indices as well)_

Variants may be filtered out of the VAT (that were in the VDS) for the following reasons:

- they are hard-filtered out based on the initial soft filtering from the GVS extract (site- and GT-level filtering)
- they have excess alternate alleles, currently that cut off is 50 alternate alleles
- they are spanning deletions
- they are duplicate variants; they are tracked via the `GvsCreateVATfromVDS` workflow's scattered `RemoveDuplicatesFromSites` task and then merged into one file by the `MergeTsvs` task
