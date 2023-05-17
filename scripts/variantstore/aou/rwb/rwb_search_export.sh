BIGQUERY_EXTRACT_SCRIPT="/tmp/rwb_export.sql"

usage() {

  echo "
Script to extract CSV files from a BigQuery GVS dataset to support search for Researcher Workbench.

Usage: $(basename "$0") --project-id <BigQuery project id> --dataset <BigQuery dataset>
                        --export-path <GCS export path> --vat-table <VAT table name>
                        [--chromosome <1-22, X or Y>] [--overwrite] [--dry-run] [--xtrace]

All parameters except --chromosome, --overwrite, and --xtrace are required.

If --chromosome is not specified, data will be extracted for all chromosomes.
If --overwrite is not specified and data exists at the specified --export-path, export will fail.
--dry-run will write a file containing the BigQuery query to '${BIGQUERY_EXTRACT_SCRIPT}' but will not actually run it.
--xtrace optionally turns on debug logging.

" 1>&2
  exit 1
}

# Prepend date, time and pwd to xtrace log entries.
PS4='\D{+%F %T} \w $ '
set -o pipefail -o errexit

VALID_ARGS=$(getopt --options p:d:e:v:c:oxy --longoptions project-id:,dataset:,export-path:,vat-table:,chromosome:,overwrite,xtrace,dry-run -- "$@")
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
    -y|--dry-run)
      DRY_RUN=true
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
  echo "--export-path must be a GCS path ending with a slash." 1>&2
  usage
fi

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
  # 'read' returns non-zero on EOF, so turn off errexit temporarily
  # https://stackoverflow.com/q/4165135
  # http://www.linuxcommand.org/lc3_man_pages/readh.html#:~:text=The%20return%20code%20is%20zero%2C%20unless%20end%2Dof%2Dfile%20is%20encountered
  set +o errexit
  read -r -d '' CHROMOSOME_FILTER_WHERE_CLAUSE <<FIN
      WHERE aa.location > $((${START_CHROMOSOME_NUMBER} * ${CHROMOSOME_MULTIPLIER}))
        AND aa.location < $((${END_CHROMOSOME_NUMBER} * ${CHROMOSOME_MULTIPLIER}))
FIN
  set -o errexit
else
  # No specified chromosome means we should export everything.
  CHROMOSOME_FILTER_WHERE_CLAUSE=""
fi

set -o nounset

BIGQUERY_EXPORT_SHARDS_PATH="${BIGQUERY_EXPORT_PATH}raw_bigquery_shards/rwb_export_*.tsv"

cat > "${BIGQUERY_EXTRACT_SCRIPT}" <<FIN

EXPORT DATA OPTIONS (
    uri ='${BIGQUERY_EXPORT_SHARDS_PATH}', format = 'CSV', header = true, field_delimiter = "\t", overwrite = ${BIGQUERY_EXPORT_OVERWRITE}) AS
    SELECT DISTINCT vat.vid as vid, si.sample_name AS sample_name
    FROM \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.alt_allele\` AS aa
             JOIN \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.sample_info\` AS si
                  ON aa.sample_id = si.sample_id
             JOIN
         (SELECT vid, (CASE SPLIT(vid, '-')[OFFSET(0)]
                         WHEN 'X' THEN 23
                         WHEN 'Y' THEN 24
                         ELSE CAST(SPLIT(vid, '-')[OFFSET(0)] AS int64) END) * ${CHROMOSOME_MULTIPLIER} +
                 CAST(SPLIT(vid, '-')[OFFSET(1)] AS int64) AS location,
                 SPLIT(vid, '-')[OFFSET(2)] AS ref_allele,
                 SPLIT(vid, '-')[OFFSET(3)] AS alt_allele
          FROM \`${BIGQUERY_PROJECT_ID}.${BIGQUERY_DATASET}.${VAT_TABLE}\`) AS vat
          ON
            vat.ref_allele = aa.ref AND
            vat.alt_allele = aa.allele AND
            vat.location = aa.location
${CHROMOSOME_FILTER_WHERE_CLAUSE}
    ORDER BY
      CAST((CASE SPLIT(vat.vid, '-')[OFFSET(0)]
        WHEN 'X' THEN '23'
        WHEN 'Y' THEN '24'
        ELSE SPLIT(vat.vid, '-')[OFFSET(0)] END) AS int64),
      CAST(SPLIT(vat.vid, '-')[OFFSET(1)] AS int64),
      SPLIT(vat.vid, '-')[OFFSET(2)],
      SPLIT(vat.vid, '-')[OFFSET(3)],
      sample_name

FIN

set +o nounset
if [[ -z "${DRY_RUN}" ]]
then
  bq query --nouse_legacy_sql --project_id "${BIGQUERY_PROJECT_ID}" < "${BIGQUERY_EXTRACT_SCRIPT}"
  echo "Files exported to '${BIGQUERY_EXPORT_SHARDS_PATH}'." 1>&2
else
  echo "Query file written to '${BIGQUERY_EXTRACT_SCRIPT}'." 1>&2
fi
