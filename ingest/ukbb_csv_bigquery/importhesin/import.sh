#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset
#pass in
BUCKET=$1 #e.g. "gs://ml4cvd/projects/jamesp/bigquery/201903"
DATASET=$2 #e.g. "ukbb7089_201903"

#specific to this func
GEO="US"
__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


#for lubitz, must replace hesin.json with hesin_lubitz.json
for NAME in hesin hesin_diag10 hesin_diag9 hesin_oper
do
    bq --location=${GEO} load \
       --field_delimiter "\t" \
       --replace \
       --source_format=CSV \
       --skip_leading_rows 1 \
       ${DATASET}.${NAME} ${BUCKET}/${NAME}.tsv.gz \
       ${__dir}/${NAME}.json
done
