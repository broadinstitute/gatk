#!/usr/bin/env bash
set -o errexit
set -o pipefail
set -o nounset

__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#RAW DATA LOCATIONS
#field file should already be in cloud before running this. looks like
# field_id
# val1
# val2
# etc
# upload using gsutil cp fieldids.csv gs://ml4h/data/
# depending on access pattern, we can make this less one-offy.
FIELDID_FILE="gs://ml4cvd/data/fieldids.csv"

#SHARED_DATASET -- should already be created, location of shared data across UKBB applications
SHARED_DATA="shared_data"


bq load \
    --replace \
    --source_format=CSV \
    --skip_leading_rows 1 \
    --schema  ${__dir}/fieldids.json \
    ${SHARED_DATA}.tensorization_fieldids ${FIELDID_FILE}
