# Creating the Variant Annotations Table

### The VAT pipeline is a set of WDLs
- [GvsCreateVAT.wdl](scripts/variantstore/wdl/GvsCreateVAT.wdl)
- [GvsValidateVAT.wdl](scripts/variantstore/variant_annotations_table/GvsValidateVAT.wdl)

The pipeline takes in a jointVCF and outputs a table in BigQuery.

**GvsCreateVAT** creates the table, while
**GvsValidateVAT** checks and validates the VAT


### Run GvsCreateVAT:

Most of the inputs are constants---like the reference, or a table schema--and dont require any additional work (and paths can be found in the [example inputs json](scripts/variantstore/wdl/GvsCreateVAT.example.inputs.json)). However for the specific data being put in the VAT, three inputs need to be created.
The first two of these inputs are two files—one of the file/vcf/shards you want to use for the VAT, and one of their indices. These are labelled as `inputFileofFileNames` and `inputFileofIndexFileNames`.

You will need to know where all of the GVS extract output files live.
Example of creating the files of file names and indices:
`gsutil ls gs://prod-drc-broad/beta-release-99k-v3/beta_99k_*.vcf.gz > inputFileofFileNames.txt`
`gsutil ls gs://prod-drc-broad/beta-release-99k-v3/beta_99k_*.vcf.gz.tbi > inputFileofIndexFileNames.txt`

Those two created files will then need to get copied into a GCP bucket that this pipeline will have access to (eg. this bucket: `gs://aou-genomics-curation-prod-processing/vat/`)

The third input is the ancestry file that Lee has made and it needs to be copied into a GCP bucket that this pipeline will have access to (eg. this bucket: `gs://aou-genomics-curation-prod-processing/vat/`) for easy access during the workflow. This input has been labelled as the `ancestry_file`.
Most of the other files are specific to where the VAT will live, like the project_id and dataset_name and the table_suffix which will name the VAT itself as vat_`table_suffix` as well as a GCP bucket location, the output_path, for the intermediary files and the VAT export in tsv form.

The other inputs can be copied from where I’ve put several files (maybe we should find a better place?!? Since many are in my scratch dir for now—tho anything in my scratch dir is also living in github as it’s things like what BQ schema is needed)

The provided [example inputs json](scripts/variantstore/wdl/GvsCreateVAT.example.inputs.json) indicates all of the needed inputs.







### Brief overview of the pipeline itself:

- `MakeSubpopulationFiles` grabs Lee’s ancestry file and keeps only the columns it needs. It also grabs the File of shard names / paths. Each of those paths is then used to kick off one of the next step which is a subworkflow of three tasks. There will be a subworkflow for each vcf shard.

  - `ExtractAnAcAfFromVCF` is the first step of the subworkflow and filters out bad sites, sites with too many alt alleles, duplicate variants and calculates the AC, AN and AF values for each subpopulation
  - `AnnotateVCF` is the second step of the subworkflow and adds Nirvana annotations to each of the variants
  - `PrepAnnotationJson` is the final step and prepares the data to be loaded into a BQ table, and loads it into a GCP staging bucket, using the output_path

- `BigQueryLoadJson` creates BQ tables (two temporary and one that is the VAT) and then gets everything from the staging buckets (`genes` and `vt`) and loads it into the two temp tables in BQ. Those two tables are joined together to create the VAT.

- `BigQuerySmokeTest` checks that the number of variants in the new VAT is the same as the expected number from the many subworkflows. This will never fail the pipeline, but is helpful to validate results.

- `BigQueryExportVat` takes the VAT and exports it by chromosome into a series of tsv files and loads it into a GCP bucket in the directory created by the `output_path`






### Notes:

When running this pipeline, currently I try to stay around 5k-10k shards—I have yet to successfully run more than that, but I think it’s more about nerves than anything else. At 20k shards, the shards do tend to step on each other if there are too many—not that anything fails or is wrong, but I’d rather get results for a run that I then do multiple times than wait for queued jobs.

Note that there are two temporary tables that are created in addition to the main VAT table: the Genes and VT tables. They have a time to live of 24 hours.
The VAT table is created by that query fresh each time so that there is no risk of duplicates.  
HOWEVER the Genes and VT tables are not. They are cleaned up after 24 hours, but this code needs to be tweaked so that you can’t get into a state where duplicates are created here. The real question here is going to be, is there a use case that we might want to run where adding to a VAT that was created say weeks ago is beneficial, but given that calculations occur on a sample summing level, this seems unlikely   






To check that all of the shards have successfully made it past the first of the sub workflow steps (the most complicated and likely to fail) they will need to be annotated / transformed into json files and put here:
`gsutil ls  [output_path]/annotations/  | wc -l`

And then once they have been transformed by the python script and are ready to be loaded into BQ they will be here:
`gsutil ls  [output_path]/genes/  | wc -l`
`gsutil ls  [output_path]/vt/  | wc -l`

These numbers are cumulative. Also the names of these json files are retained from the original shard names so as to not cause collisions. If you run the same shards through the VAT twice, the second runs should overwrite the first and the total number of jsons should not change.
Once the shards have make it into the /genes/ and /vt/ directories, the majority of the expense and transformations needed for that shard are complete.
They are ready to be loaded into BQ. You will notice that past this step, all there is to do is create the BQ tables, load the BQ tables, run a join query and then the remaining steps are all validations or an export into tsv.


There are often a fair number of failures from google in this workflow—-so far they have all been 503s. Because of the re-working of the workflow, they should not interrupt unaffected shards, but the shards that are affected will need to be collected and put into a new File of File names and re-run.
In theory you could actually just re-run the pipeline entirely and since the 503s seem to be completely random and intermittent, the same shards likely wont fail twice and you’d get everything you need, but at low efficiency and high expense.
To grab any files not in the bucket, but in the file of file names:

`gsutil ls [output_path]/genes/ | awk '{print substr($0,90)}' RS='.json.gz\n'  > successful_vcf_shards.txt`  
`gsutil cat [inputFileofFileNames] | awk '{print substr($0,41)}' RS='.vcf.gz\n' > all_vcf_shards.txt`  
`comm -3 all_vcf_shards.txt successful_vcf_shards.txt > diff_vcf_shards.txt`  
`awk '{print g[bucket GVS output VCF shards are in]" $0 ".vcf.gz"}' diff_vcf_shards.txt > fofn_missing_shards.txt`  

_(dont forget to do this for the indices as well)_




There is a line of code in ExtractAnAcAfFromVCF (the most expensive in $ and time task in the workflow) that can be removed because it is used to track the variants that are dropped _TODO: Rori to make using it a parameter_















