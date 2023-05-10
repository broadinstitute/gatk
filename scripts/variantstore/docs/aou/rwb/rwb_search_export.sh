# Prepend date, time and pwd to xtrace log entries.
PS4='\D{+%F %T} \w $ '
set -o pipefail -o errexit

usage() {
  echo "Usage: $(basename "$0") --project-id <BigQuery project id> --dataset <BigQuery dataset> --export-path <GCS export path> --vat-table <VAT table name> [--chromosome <1-22, X or Y>] [--overwrite] [--xtrace]" 1>&2
  echo "All parameters except --chromosome, --overwrite, and --xtrace are mandatory." 1>&2
  exit 1
}

VALID_ARGS=$(getopt --options p:d:e:v:c:ox --longoptions project-id:,dataset:,export-path:,vat-table:,chromosome:,overwrite,xtrace -- "$@")
if [[ $? -ne 0 ]]
then
  usage
fi

unset -v BIGQUERY_PROJECT_ID BIGQUERY_DATASET BIGQUERY_EXPORT_PATH CHROMOSOME
BIGQUERY_EXPORT_OVERWRITE=false

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
    -e|--export-path)
      BIGQUERY_EXPORT_PATH="$2"
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
    -o|--overwrite)
      BIGQUERY_EXPORT_OVERWRITE=true
      shift
      ;;
    -x|--xtrace)
      set -o xtrace
      shift
      ;;
    --) shift;
      break
      ;;
  esac
done

if [[ -z "${BIGQUERY_PROJECT_ID}" ]] || [[ -z "${BIGQUERY_DATASET}" ]] || [[ -z "${BIGQUERY_EXPORT_PATH}" ]] || [[ -z "${VAT_TABLE}" ]]
then
  usage
fi

if ! [[ "${BIGQUERY_EXPORT_PATH}" =~ ^gs:// ]]
then
  echo "--export-path must be a GCS path beginning with 'gs://'." 1>&2
  usage
elif ! [[ "${BIGQUERY_EXPORT_PATH}" =~ /$ ]]
then
  echo "--extport-path must be a GCS path ending with a slash." 1>&2
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
  END_CHROMOSOME_NUMBER=$((${START_CHROMOSOME_NUMBER} + 1))
  # 'read' apparently returns non-zero
  set +o errexit
  read -r -d '' EXPORT_WHERE_CLAUSE <<FIN
      WHERE aa.location > $((${START_CHROMOSOME_NUMBER} * ${CHROMOSOME_MULTIPLIER}))
        AND aa.location < $((${END_CHROMOSOME_NUMBER} * ${CHROMOSOME_MULTIPLIER}))
FIN
  set -o errexit
else
  # No specified chromosome means we should export everything
  EXPORT_WHERE_CLAUSE=""
fi

BIGQUERY_EXPORT_SHARDS_PATH="${BIGQUERY_EXPORT_PATH}raw_bigquery_shards/rwb_export_*.tsv"

cat > "${BIGQUERY_EXTRACT_SCRIPT}" <<FIN

EXPORT DATA OPTIONS (
    uri ='${BIGQUERY_EXPORT_SHARDS_PATH}', format = 'CSV', header = true, field_delimiter = "\t", overwrite = ${BIGQUERY_EXPORT_OVERWRITE}) AS
    SELECT (CASE CAST(aa.location / ${CHROMOSOME_MULTIPLIER} AS int64)
                WHEN 23 THEN 'X'
                WHEN 24 THEN 'Y'
                ELSE CAST(CAST(aa.location / ${CHROMOSOME_MULTIPLIER} AS int64) AS string) END) AS chromosome,
           MOD(aa.location, ${CHROMOSOME_MULTIPLIER})                                           AS position,
           aa.ref AS ref,
           aa.allele AS allele,
           si.sample_name AS sample_name
    FROM \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.alt_allele\` AS aa
             JOIN \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.sample_info\` AS si
                  ON aa.sample_id = si.sample_id
             JOIN
         (SELECT ((CASE SPLIT(vid, '-')[OFFSET(0)]
                       WHEN 'X' THEN 23
                       WHEN 'Y' THEN 24
                       ELSE CAST(SPLIT(vid, '-')[OFFSET(0)] AS int64) END) * ${CHROMOSOME_MULTIPLIER} +
                 CAST(SPLIT(vid, '-')[OFFSET(1)] AS int64))  AS location,
                 SPLIT(vid, '-')[OFFSET(2)]                  AS ref,
                 SPLIT(vid, '-')[OFFSET(3)]                  AS alt
          FROM \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.${VAT_TABLE}\`)
             AS vat
         ON vat.ref = aa.ref AND vat.alt = aa.allele AND vat.location = aa.location
${EXPORT_WHERE_CLAUSE}
    ORDER BY aa.location, ref, allele, sample_name

FIN

bq query --nouse_legacy_sql --project_id "${BIGQUERY_PROJECT_ID}" < "${BIGQUERY_EXTRACT_SCRIPT}"

echo "Files exported to '${BIGQUERY_EXPORT_SHARDS_PATH}'." 1>&2
