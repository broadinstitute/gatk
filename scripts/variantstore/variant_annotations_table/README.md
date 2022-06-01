# Creating the Variant Annotations Table

### The VAT pipeline is a set of WDLs
- [GvsCreateVAT.wdl](/scripts/variantstore/wdl/GvsCreateVAT.wdl)
- [GvsValidateVAT.wdl](/scripts/variantstore/variant_annotations_table/GvsValidateVAT.wdl)

The pipeline takes in a jointVCF and outputs a table in BigQuery.

**GvsCreateVAT** creates the table, while
**GvsValidateVAT** checks and validates the created VAT and prints a report of any failing validation


### Run GvsCreateVAT:

Most of the inputs are constants — like the reference, or a table schema — and don't require any additional work, or have defaults set in the WDL). However for the specific data being put in the VAT, three inputs need to be created.

The first two of these inputs are two files — one of the file/vcf/shards you want to use for the VAT, and their corresponding index files. These are labelled as `inputFileofFileNames` and `inputFileofIndexFileNames` and need to be copied into a GCP bucket that this pipeline will have access to (eg. this bucket: `gs://aou-genomics-curation-prod-processing/vat/`) for easy access during the workflow.
The third input is the ancestry file from the ancestry pipeline which will be used to calculate AC, AN and AF for all subpopulations. It needs to be copied into a GCP bucket that this pipeline will have access to. This input has been labelled as the `ancestry_file`.

Most of the other files are specific to where the VAT will live, like the project_id and dataset_name and the table_suffix which will name the VAT itself as vat_`table_suffix` as well as a GCP bucket location, the output_path, for the intermediary files and the VAT export in tsv form.

All optional inputs are provided with default values.


### Notes:

Note that there are two temporary tables that are created in addition to the main VAT table: the Genes and VT tables. They have a time to live of 24 hours.
The VAT table is created by that query fresh each time so that there is no risk of duplicates.  
HOWEVER the Genes and VT tables are not. They are cleaned up after 24 hours, but this code needs to be tweaked so that you can’t get into a state where duplicates are created here. The real question here is going to be, is there a use case that we might want to run where adding to a VAT that was created say weeks ago is beneficial, but given that calculations occur on a sample summing level, this seems unlikely.


To check that all of the shards have successfully made it past the first of the sub workflow steps (the most complicated and likely to fail) they will need to be annotated / transformed into json files and put here:
`gsutil ls  [output_path]/annotations/  | wc -l`

And then once they have been transformed by the python script and are ready to be loaded into BQ they will be here:
`gsutil ls  [output_path]/genes/  | wc -l`
`gsutil ls  [output_path]/vt/  | wc -l`

These numbers are cumulative. Also the names of these json files are retained from the original shard names so as to not cause collisions. If you run the same shards through the VAT twice, the second runs should overwrite the first and the total number of jsons should not change.
Once the shards have make it into the /genes/ and /vt/ directories, the majority of the expense and transformations needed for that shard are complete.
They are ready to be loaded into BQ. You will notice that past this step, all there is to do is create the BQ tables, load the BQ tables, run a join query and then the remaining steps are all validations or an export into tsv.


Sometimes shards will fail with a 503 from Google. The shards that are affected will need to be collected and put into a new File of File names and re-run.
In theory you could actually just re-run the pipeline entirely and since the 503s seem to be completely random and intermittent, the same shards likely wont fail twice and you’d get everything you need, but at low efficiency and high expense.
To grab any files not in the bucket, but in the file of file names:

`gsutil ls [output_path]/genes/ | awk '{print substr($0,90)}' RS='.json.gz\n'  > successful_vcf_shards.txt`  
`gsutil cat [inputFileofFileNames] | awk '{print substr($0,41)}' RS='.vcf.gz\n' > all_vcf_shards.txt`  
`comm -3 all_vcf_shards.txt successful_vcf_shards.txt > diff_vcf_shards.txt`  
`awk '{print g[bucket GVS output VCF shards are in]" $0 ".vcf.gz"}' diff_vcf_shards.txt > fofn_missing_shards.txt`  

_(dont forget to do this for the indices as well)_

All dropped variants are tracked













