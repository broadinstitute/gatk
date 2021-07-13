This directory has the WDLs for generating and interacting with a variant store in BigQuery. It contains
`ImportArrayManifest.wdl` which uploads the `probe_info` table from a Manifest file, `ImportArrays.wdl` which uploads
 the raw array data from a list of vcfs, and `raw_array_cohort_extract.wdl` 
with example inputs `raw_array_cohort_extract_inputs.json` which extracts cohorts from a fully uploaded dataset. The 
docker image required for the extract WDL can be generated from the Dockerfile in the `extract/` directory. There is no
WDL for generating the `genotype_counts` table, but this can be achieved from running the `raw_array_cohort_extract.py` 
script in the `extract` directory with the `--extract_genotype_counts_only` argument set to `true`.
