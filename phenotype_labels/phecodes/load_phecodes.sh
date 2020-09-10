#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

#RAW DATA LOCATOINS
#beta version of ICD10 (WHO) -> phecodes mapping
#https://phewascatalog.org/phecodes_icd10
#There is an ICD10CM (US) version as well
#Phecode definition map is here
#https://phewascatalog.org/phecodes

#SHARED_DATASET -- should already be created
SHARED_DATA=shared_data
#specific to this func
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#Phecode CSVs (beta 2.2 release) are located in 
ICD10_FILE="gs://ml4h/data/phecodes/Phecode_map_v1_2_icd10_beta.csv.gz"
PHECODE_DEF_FILE="gs://ml4h/data/phecodes/phecode_definitions1.2.csv.gz"

#phecode definitions
bq load \
 --replace \
 --source_format=CSV \
 --skip_leading_rows 1 \
 --schema  ${__dir}/phecode_dictionary.json \
 ${SHARED_DATA}.phecode_dictionary ${PHECODE_DEF_FILE}
        
#icd10
bq load \
 --replace \
 --source_format=CSV \
 --skip_leading_rows 1 \
 --format=prettyjson \
 --schema  ${__dir}/icd10.json \
 ${SHARED_DATA}.phecode_icd10 ${ICD10_FILE}  

