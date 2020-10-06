#!/usr/bin/env bash
set -e

if [ $# -lt 5 ]; then
  echo "usage: $0 <project-id> <dataset-name> <storage-location> <table-id> <load> <uuid>"
  exit 1
fi

PROJECT_ID=$1
DATASET_NAME=$2
STORAGE_LOCATION=$3
TABLE_ID=$4
if [ $5 == "true" ]; then
  LOAD=true
else
  LOAD=false
fi
if [ $# -eq 6 ]; then
  UUID_FOR_TABLE="${6}_"
else
  UUID_FOR_TABLE=""
fi
SAMPLE_DIR=$STORAGE_LOCATION/sample_tsvs/
RAW_DIR=$STORAGE_LOCATION/raw_tsvs/

let "PARTITION_START=($TABLE_ID-1)*4000+1"
let "PARTITION_END=$PARTITION_START+3999"
let "PARTITION_STEP=1"
PARTITION_FIELD="sample_id"
printf -v PADDED_TABLE_ID "%03d" $TABLE_ID

RAW_FILES="raw_${PADDED_TABLE_ID}_*"
METADATA_FILES="sample_${PADDED_TABLE_ID}_*"

NUM_RAW_FILES=$(gsutil ls $RAW_DIR${RAW_FILES} | wc -l)
NUM_METADATA_FILES=$(gsutil ls $SAMPLE_DIR${METADATA_FILES} | wc -l)

if [ $NUM_RAW_FILES -eq 0 -a $NUM_METADATA_FILES -eq 0 ]; then
  "no files for table ${PADDED_TABLE_ID} to process in $STORAGE_LOCATION; exiting"
  exit
fi

# schema and TSV header need to be the same order
RAW_SCHEMA="schemas/raw_array_schema.json"
SAMPLE_LIST_SCHEMA="schemas/arrays_sample_list_schema.json"

# create a metadata table and load
SAMPLE_LIST_TABLE="${DATASET_NAME}.${UUID_FOR_TABLE}sample_list"
if [ $NUM_METADATA_FILES -gt 0 ]; then
  set +e
  bq ls --project_id $PROJECT_ID $DATASET_NAME > /dev/null
  set -e
  if [ $? -ne 0 ]; then
    echo "making dataset $DATASET_NAME"
    bq mk --project_id=$PROJECT_ID $DATASET_NAME
  fi
  set +e
  bq show --project_id $PROJECT_ID $SAMPLE_LIST_TABLE > /dev/null
  set -e
  if [ $? -ne 0 ]; then
    echo "making table $SAMPLE_LIST_TABLE"
    bq --location=US mk --project_id=$PROJECT_ID $SAMPLE_LIST_TABLE $SAMPLE_LIST_SCHEMA
    #TODO: add a Google Storage Transfer for the table when we make it.
  fi
  if [ "$LOAD" = true ]; then
    bq load --location=US --project_id=$PROJECT_ID --skip_leading_rows=1 --null_marker="null" --source_format=CSV -F "\t" $SAMPLE_LIST_TABLE $SAMPLE_DIR$METADATA_FILES $SAMPLE_LIST_SCHEMA
    echo "ingested ${METADATA_FILES} file from $SAMPLE_DIR into table $SAMPLE_LIST_TABLE"
  else
    echo "${METADATA_FILES} will be ingested from $SAMPLE_DIR by Google Storage Transfer"
  fi
else
  echo "no metadata files to process"
fi

# create array table
TABLE="${DATASET_NAME}.${UUID_FOR_TABLE}arrays_${PADDED_TABLE_ID}"
if [ $NUM_RAW_FILES -gt 0 ]; then
  set +e
  bq show --project_id $PROJECT_ID $TABLE > /dev/null
  set -e
  if [ $? -ne 0 ]; then
    echo "making table $TABLE"
    bq --location=US mk --range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP \
      --project_id=$PROJECT_ID $TABLE $RAW_SCHEMA
    #TODO: add a Google Storage Transfer for the table when we make it.
  fi
  if [ "$LOAD" = true ]; then
    bq load --location=US --project_id=$PROJECT_ID --skip_leading_rows=1 --null_marker="null" --source_format=CSV -F "\t" $TABLE $RAW_DIR$RAW_FILES $RAW_SCHEMA
    echo "ingested ${RAW_FILES} files from $RAW_DIR into table $TABLE"
  else
    echo "${RAW_FILES} will be ingested from $RAW_DIR
     by Google Storage Transfer"
  fi
else
  echo "no raw data files to process"
fi

