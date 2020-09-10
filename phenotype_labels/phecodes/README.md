# Phecodes

Phecodes are a mapping of ICD (9 and 10) codes to a hierarchy of diseases.

To put Phecodes on top of a UKBB dataset (ingested into bigquery according to the ml4h ingest script) you have to 
1. follow the instructions in load_phecodes.sh to load the raw phecode information (icd9 -> phenotype) into the database.
2. run map_phecodes.py to create a local csv file with mappings of the HESIN table to phecodes (the dataset is hardcoded for now)
3. load that csv file into bigquery manually, for example with 

```
gsutil cp ukbb_dev_phecode_mapping.csv.gz gs://ml4cvd/projects/pbatra/ukbb_dev/
bq load \
 --replace \
 --source_format=CSV \
 --schema  phecode_mapping.json \
 ukbb_dev.phecode_mapping gs://ml4cvd/projects/pbatra/ukbb_dev/ukbb_dev_phecode_mapping.csv.gz
```
