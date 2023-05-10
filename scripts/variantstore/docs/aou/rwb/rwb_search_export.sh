#!/usr/bin/env bash

# Prepend date, time and pwd to xtrace log entries.
PS4='\D{+%F %T} \w $ '
set -o errexit -o pipefail -o xtrace

usage() {
  echo "Usage: $(basename "$0") --project-id <BigQuery project id> --dataset <BigQuery dataset> --extract-path <GCS extract path> --vat-table <VAT table name> [--chromosome <1-22, X or Y>]" 1>&2
  echo "All parameters except --chromosome are mandatory." 1>&2
  exit 1
}

VALID_ARGS=$(getopt --options p:d:e:v:c: --long project-id:dataset:extract-path:vat-table:chromosome: -- "$@")
if [[ $? -ne 0 ]]; then
  usage
fi

unset -v BIGQUERY_PROJECT_ID BIGQUERY_DATASET GCS_EXTRACT_PATH CHROMOSOME

eval set -- "$VALID_ARGS"
while true
do
  case "$1" in
    -p|--project-id)
      BIGQUERY_PROJECT_ID="$2"
      shift 2
      ;;
    -d|--dataset)
      BIGQUERY_DATASET="$2"
      shift 2
      ;;
    -e|--extract-path)
      GCS_EXTRACT_PATH="$2"
      shift 2
      ;;
    -v|--vat-table)
      VAT_TABLE="$2"
      shift 2
      ;;
    -c|--chromosome)
      CHROMOSOME="$2"
      shift 2
      ;;
    --) shift;
      break
      ;;
  esac
done

if [[ -z "${BIGQUERY_PROJECT_ID}" ]] || [[ -z "${BIGQUERY_DATASET}" ]] || [[ -z "${GCS_EXTRACT_PATH}" ]] || [[ -z "${VAT_TABLE}" ]]
then
  usage
fi

if ! [[ "${GCS_EXTRACT_PATH}" =~ "^gs://" ]]
then
  echo "--extract-path must be a GCS path beginning with 'gs://'." 1>&2
  usage
elif ! [[ "${GCS_EXTRACT_PATH}" =~ "/$" ]]
then
  echo "--extract-path must be a GCS path ending with a slash." 1>&2
  usage
fi

set -o nounset

BIGQUERY_EXTRACT_SCRIPT="/tmp/rwb_export.sql"

unset -v START_CHROMOSOME END_CHROMOSOME CHROMOSOME_WHERE_CLAUSE

if [[ -n "${CHROMOSOME}" ]]
then
  case "${CHROMOSOME}" in
    x|X)
      START_CHROMOSOME_NUMBER=23
      ;;
    y|Y)
      START_CHROMOSOME_NUMBER=24
      ;;
    1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22)
      START_CHROMOSOME_NUMBER="${CHROMOSOME}"
      ;;
    *)
      echo "Unrecognized --chromosome value '$CHROMOSOME', should be 1-22, X or Y." 1>&2
      usage
  esac
fi

CHROMOSOME_MULTIPLIER=1000000000000

if [[ -n "${START_CHROMOSOME_NUMBER}" ]]
then
  END_CHROMOSOME_NUMBER=$(("${START_CHROMOSOME_NUMBER}" + 1))
  read -r -d '' EXPORT_WHERE_CLAUSE <<FIN
      where aa.location > $((START_CHROMOSOME_NUMBER * CHROMOSOME_MULTIPLIER))
        and aa.location < $((END_CHROMOSOME_NUMBER * CHROMOSOME_MULTIPLIER))
FIN
else
  # No specified chromosome means we should export everything
  EXPORT_WHERE_CLAUSE=""
fi


cat > "${BIGQUERY_EXTRACT_SCRIPT}" <<FIN

EXPORT DATA OPTIONS (
    uri ='${GCS_EXTRACT_PATH}rwb_export_*.tsv', format ='CSV', header = false, field_delimiter ="\t") as
    SELECT (case cast(aa.location / ${CHROMOSOME_MULTIPLIER} as int64)
                when 23 then 'X'
                when 24 then 'Y'
                else cast(cast(aa.location / ${CHROMOSOME_MULTIPLIER} as int64) as string) end) as chromosome,
           MOD(aa.location, ${CHROMOSOME_MULTIPLIER})                                           as position,
           aa.ref,
           aa.allele,
           si.sample_name
    FROM \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.alt_allele\` as aa
             join \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.sample_info\` as si
                  on aa.sample_id = si.sample_id
             join
         (SELECT ((case SPLIT(vid, '-')[OFFSET(0)]
                       when 'X' then 23
                       when 'Y' then 24
                       else cast(SPLIT(vid, '-')[OFFSET(0)] AS int64) end) * ${CHROMOSOME_MULTIPLIER} +
                  CAST(SPLIT(vid, '-')[OFFSET(1)] AS int64)) as location,
                 SPLIT(vid, '-')[OFFSET(2)]                  as ref,
                 SPLIT(vid, '-')[OFFSET(3)]                  as alt
          FROM \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.${VAT_TABLE}\`)
             as vat
         on vat.ref = aa.ref and vat.alt = aa.allele and vat.location = aa.location
${EXPORT_WHERE_CLAUSE}
    order by aa.location

FIN
