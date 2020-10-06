#!/usr/bin/env bash
set -e

if [ $# -lt 4 ]; then
  echo "usage: $0 <project-id> <dataset-name> <storage-location> <table-id>"
  exit 1
fi

PROJECT_ID=$1
DATASET_NAME=$2
STORAGE_LOCATION=$3
TABLE_ID=$4
PROCESSING_DIR=$STORAGE_LOCATION/
DONE_DIR=$STORAGE_LOCATION/done/

let PARTITION_START=($TABLE_ID-1)*4000+1
let PARTITION_END=$PARTITION_START+3999
printf -v PADDED_TABLE_ID "%03d" $TABLE_ID

RAW_FILES="raw_${PADDED_TABLE_ID}_*"
METADATA_FILES="sample_${PADDED_TABLE_ID}_*"

NUM_RAW_FILES=$(gsutil ls ${PROCESSING_DIR}${RAW_FILES} | wc -l)
NUM_METADATA_FILES=$(gsutil ls $PROCESSING_DIR${METADATA_FILES} | wc -l)

if [ $NUM_RAW_FILES -eq 0 -a $NUM_METADATA_FILES -eq 0 ]; then
  "no files for table ${PADDED_TABLE_ID} to process in $PROCESSING_DIR; exiting"
  exit
fi

# schema and TSV header need to be the same order
RAW_SCHEMA="raw_array_schema.json"
SAMPLE_LIST_SCHEMA="arrays_sample_list_schema.json"

# create a metadata table and load
SAMPLE_LIST_TABLE=$DATASET_NAME.sample_list
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
  fi
  bq load --location=US --project_id=$PROJECT_ID --skip_leading_rows=1 --null_marker="null" --source_format=CSV -F "\t" $SAMPLE_LIST_TABLE $PROCESSING_DIR$METADATA_FILES $SAMPLE_LIST_SCHEMA
  echo "ingested ${METADATA_FILES} file from $PROCESSING_DIR into table $SAMPLE_LIST_TABLE"
else
  echo "no metadata files to process"
fi

# create array table
TABLE=$DATASET_NAME.arrays_$PADDED_TABLE_ID
if [ $NUM_RAW_FILES -gt 0 ]; then
  set +e
  bq show --project_id $PROJECT_ID $TABLE > /dev/null
  set -e
  if [ $? -ne 0 ]; then
    echo "making table $TABLE"
    bq --location=US mk --range_partitioning=$PARTITION_FIELD,$PARTITION_START,$PARTITION_END,$PARTITION_STEP \
      --project_id=$PROJECT_ID $TABLE $RAW_SCHEMA
  fi
  bq load --location=US --project_id=$PROJECT_ID --skip_leading_rows=1 --null_marker="null" --source_format=CSV -F "\t" $TABLE $PROCESSING_DIR$RAW_FILES $RAW_SCHEMA
  echo "ingested ${RAW_FILES} files from $PROCESSING_DIR into table $TABLE"
else
  echo "no raw data files to process"
fi
echo "moving files from processing to done"
gsutil -q -m mv $PROCESSING_DIR$METADATA_FILES $DONE_DIR 
gsutil -q -m mv $PROCESSING_DIR$RAW_FILES $DONE_DIR 
