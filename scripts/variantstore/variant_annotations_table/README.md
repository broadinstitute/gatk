
# Creating the Variant Annotations Table

### The VAT pipeline is a set of WDLs
- [GvsCreateVATfromVDS.wdl](/scripts/variantstore/wdl/GvsCreateVATfromVDS.wdl)
- [GvsValidateVAT.wdl](/scripts/variantstore/variant_annotations_table/GvsValidateVAT.wdl)

The pipeline takes in jointVCFs (and their index files), creates a queryable table in BigQuery, and outputs a bgzipped TSV file containing the contents of that table. That TSV file is the sole AoU deliverable.

**GvsCreateVATfromVDS** creates the variant annotations table from a sites only VCF and a TSV both created from a VDS
**GvsValidateVAT** checks and validates the created VAT and prints a report of any failing validation


### Run GvsCreateVATfromVDS:

Most of the inputs are constants — like the reference, or a table schema — and don't require any additional work, or have defaults set in the WDL). 
However for the specific data being put in the VAT, three inputs need to be created.

The first of these inputs is a Sites Only VCF
The second of these inputs is a TSV of custom annotations (AC AN and AF)
The third input is the ancestry file from the ancestry pipeline which will be used to calculate AC, AN and AF for all subpopulations. It needs to be copied into a GCP bucket that this pipeline will have access to. This input has been labelled as the `ancestry_file`.

Most of the other files are specific to where the VAT will live, like the project_id and dataset_name and the table_suffix which will name the VAT itself as vat_`table_suffix` as well as a GCP bucket location, the output_path, for the intermediary files and the VAT export in TSV form.

All optional inputs are provided with default values.

### Preparing to run GvsCreateVATfromVDS:

Two inputs need to be created from the VDS.
The first input can be created using the ancestry pipeline which is in another workspace. This input is also needed to create the second input.
The second input can be creating using this Python script:`scripts/variantstore/wdl/extract/hail_create_vat_inputs.py`
Making the second input:
`hail_create_vat_inputs.py` script gives us a Sites Only VCF which can then be used to create a custom annotations file for Nirvana.
This Python script needs the following inputs:
* VDS -- this can be created using the GVS callset creation pipeline
* Ancestry TSV -- this can be created using the Ancestry pipeline (and this is also the first input needed to run the WDL)
We recommend using a Terra Python Notebook with a small Hail cluster to run the Python script and make these VAT pipeline inputs.
Once those inputs have been created, they can be used by the GvsCreateVATfromVDS WDL

  These are the required parameters which must be supplied to the GvsCreateVATfromVDS workflow:

| Parameter            | Description                                           |
|----------------------|-------------------------------------------------------|
| ancestry_file        | the GCP path to the ancestry TSV                      |
| dataset_name         | the name of the dataset you created above             |
| filter_set_name      | the prefix of the VAT                                 |
| input_sites_only_vcf | the GCP path to the output of the above python script |
| output_path          | the path that the temp files will be stored           |
| project_id           | the name of the google project containing the dataset |

### Notes:

Note that there are two temporary tables that are created in addition to the main VAT table: the Genes and VT tables. They have a time to live of 24 hours.
The VAT table is created by that query fresh each time so that there is no risk of duplicates.  

To check that all of the shards have successfully made it past the first of the sub workflow steps (the most complicated and likely to fail) make sure that the full number of expected files are here:
`gsutil ls  [output_path]/annotations/  | wc -l`

And then once they have been transformed by the Python script and are ready to be loaded into BQ they will be here:
`gsutil ls  [output_path]/genes/  | wc -l`
`gsutil ls  [output_path]/vt/  | wc -l`

These numbers are cumulative. Also the names of these json files are retained from the original shard names so as to not cause collisions. If you run the same shards through the VAT twice, the second runs should overwrite the first and the total number of jsons should not change.
Once the shards have make it into the /genes/ and /vt/ directories, the majority of the expense and transformations needed for that shard are complete.
They are ready to be loaded into BQ. You will notice that past this step, all there is to do is create the BQ tables, load the BQ tables, run a join query and then the remaining steps are all validations or an export into TSV.


Sometimes shards will fail with a 503 from Google. The shards that are affected will need to be collected and put into a new File of File names and re-run.
In theory you could actually just re-run the pipeline entirely and since the 503s seem to be completely random and intermittent, the same shards likely wont fail twice and you’d get everything you need, but at low efficiency and high expense.
To grab any files not in the bucket, but in the file of file names:

`gsutil ls [output_path]/genes/ | awk '{print substr($0,90)}' RS='.json.gz\n'  > successful_vcf_shards.txt`  
`gsutil cat [inputFileofFileNames] | awk '{print substr($0,41)}' RS='.vcf.gz\n' > all_vcf_shards.txt`  
`comm -3 all_vcf_shards.txt successful_vcf_shards.txt > diff_vcf_shards.txt`  
`awk '{print g[bucket GVS output VCF shards are in]" $0 ".vcf.gz"}' diff_vcf_shards.txt > fofn_missing_shards.txt`  

_(dont forget to do this for the indices as well)_

Variants may be filtered out of the VAT (that were in the GVS extract) for the following reasons:
- they are hard-filtered out based on the initial soft filtering from the GVS extract (site and GT level filtering)
- they have excess alternate alleles--currently that cut off is 50 alternate alleles
- they are spanning deletions

## note that there are sometimes duplicate variants that are removed and tracked













