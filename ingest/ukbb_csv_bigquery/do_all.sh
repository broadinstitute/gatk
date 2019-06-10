#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

#to run
#bash do_all.sh 2>&1 | tee <name>.log
#can replace pigz with gzip -c if prefered.


#SET THESE ARGUMENTS -- not passing these in -- just gets too messy.
BUCKET_ROOT="" #gs://ml4cvd/projects/pbatra"
DATASET="" #pb_test_2" if dataset exists already, will error out.
PHENO_FOLDER="" #/Users/pbatra/ukbb_etl/downloads/ukbb_raw_data/pheno_test_files" #folder that contains multiple pheno.csv files
HESIN_FOLDER="" #/Users/pbatra/ukbb_etl/downloads/ukbb_raw_data/hesin_files" #folder that contains all hesin files
DEATH_CENSOR="" #2018-01-31" #lookup at https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=Data_providers_and_dates
PHENO_CENSOR="" #2017-03-31" #lookup at https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=Data_providers_and_dates

if [ "${BUCKET_ROOT}" == "" ] || [ "${DATASET}" == "" ] || [ "${PHENO_FOLDER}" == "" ] || \ 
    [ "${HESIN_FOLDER}" == "" ] || [ "${DEATH_CENSOR}" == "" ] || [ "${PHENO_CENSOR}" == "" ]; then
        echo "set all parameters, please"
        exit 1
fi

#print out all pheno files being loaded
pheno_files=$( ls ${PHENO_FOLDER})
echo "loading pheno csvs: 
${pheno_files}"

hesin_files="hesin.tsv.gz
hesin_diag10.tsv.gz
hesin_diag9.tsv.gz
hesin_oper.tsv.gz"

#check to make sure we have all the hesin files we expect
for file in ${hesin_files}
do
    if [ -f ${HESIN_FOLDER}/${file} ]; then
        echo "${file} exists"
    else
        echo "${file} missing" 
        exit 1
    fi
done

#STORAGE
BUCKET=${BUCKET_ROOT}/${DATASET} #gs target location for all generated data
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOCAL_DATA_FOLDER=${__dir}/${DATASET} #local storage location
bq mk --dataset broad-ml4cvd:${DATASET} #make dataset in bigquery, error if already there
mkdir ${LOCAL_DATA_FOLDER} #will error if already exists
echo "storing data locally in ${LOCAL_DATA_FOLDER}"


#dictionary
file=dictionary.tsv.gz
go run ${__dir}/convertdict/main.go | pigz > ${LOCAL_DATA_FOLDER}/${file}
gsutil cp ${LOCAL_DATA_FOLDER}/${file} ${BUCKET}/${file} 
bash ${__dir}/importdict/import.sh ${BUCKET} ${DATASET}

#coding
file=coding.tsv.gz
go run ${__dir}/convertcoding/main.go | pigz > ${LOCAL_DATA_FOLDER}/${file}
gsutil cp ${LOCAL_DATA_FOLDER}/${file} ${BUCKET}/${file} 
bash ${__dir}/importcoding/import.sh ${BUCKET} ${DATASET}

#pheno
pheno_args=""
for pheno_file in $pheno_files
do
    pheno_args="$pheno_args -pheno ${PHENO_FOLDER}/${pheno_file}"
done

go run ${__dir}/convertpheno/*.go \
    -bigquery ${DATASET} \
    ${pheno_args} \
    -ack | pigz > ${LOCAL_DATA_FOLDER}/phenotype.tsv.gz
gsutil cp ${LOCAL_DATA_FOLDER}/phenotype.tsv.gz ${BUCKET}/phenotype.tsv.gz
bash ${__dir}/importpheno/append.sh ${BUCKET} ${DATASET}

#hesin
for file in ${hesin_files}
do
    gsutil cp ${HESIN_FOLDER}/${file} ${BUCKET}/${file}
done
bash ${__dir}/importhesin/import.sh ${BUCKET} ${DATASET}

#sample
##SKIP

#censor
file=censor.tsv.gz
go run  ${__dir}/censor/*.go \
    -bigquery ${DATASET} \
    -death_censor ${DEATH_CENSOR} \
    -pheno_censor ${PHENO_CENSOR} | pigz > ${LOCAL_DATA_FOLDER}/${file}
gsutil cp ${LOCAL_DATA_FOLDER}/${file} ${BUCKET}/${file}
bash ${__dir}/importcensor/import.sh ${BUCKET} ${DATASET}
