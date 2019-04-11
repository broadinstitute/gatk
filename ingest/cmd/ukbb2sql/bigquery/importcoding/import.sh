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

bq --location=${GEO} load \
 --field_delimiter "\t" \
 --replace \
 --quote "" \
 --source_format=CSV \
 --skip_leading_rows 1 \
 ${DATASET}.coding ${BUCKET}/coding.tsv.gz \
 ${__dir}/coding.json
