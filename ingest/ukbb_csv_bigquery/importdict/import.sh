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

#note hardcoded dictionary.json, dictionary.tsv.gz, dictionary table
# Special for dict: need to disable quotes
bq --location=${GEO} load \
 --field_delimiter "\t" \
 --replace \
 --quote "" \
 --source_format=CSV \
 --skip_leading_rows 1 \
 ${DATASET}.dictionary ${BUCKET}/dictionary.tsv.gz \
 ${__dir}/dictionary.json
