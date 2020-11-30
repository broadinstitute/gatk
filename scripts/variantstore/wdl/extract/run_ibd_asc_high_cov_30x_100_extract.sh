PROJECT="spec-ops-aou"
DATASET="ibd_asc_wgs_EXAMPLE"

python ngs_cohort_extract.py \
  --fq_petvet_dataset ${PROJECT}.${DATASET} \
  --fq_temp_table_dataset ${PROJECT}.temp_tables \
  --fq_destination_dataset ${PROJECT}.${DATASET} \
  --destination_table exported_cohort_100_test \
  --fq_cohort_sample_names ${PROJECT}.${DATASET}.cohort_100_of_1084 \
  --query_project ${PROJECT} \
  --fq_sample_mapping_table ${PROJECT}.${DATASET}.metadata