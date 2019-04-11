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

# For phenotypes, we expect to add repeatedly, so we don't replace here. Note:
# if you run append.sh with the same data twice, you'll just duplicate the
# contents of the table.
bq --location=${GEO} load \
 --field_delimiter "\t" \
 --quote "" \
 --null_marker "NULL" \
 --source_format=CSV \
 --skip_leading_rows 1 \
 ${DATASET}.phenotype ${BUCKET}/phenotype.tsv.gz \
 ${__dir}/phenotype.json
