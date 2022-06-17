# -*- coding: utf-8 -*-
import uuid
import datetime
import argparse
import re

from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
from google.oauth2 import service_account

import utils

JOBS = []

#
# CONSTANTS
#
REF_TABLE_PREFIX = "ref_ranges_"
VET_TABLE_PREFIX = "vet_"
SAMPLES_PER_PARTITION = 4000

FINAL_TABLE_TTL = ""
#FINAL_TABLE_TTL = " OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 72 HOUR)) "

# temp-table-uuid
output_table_prefix = str(uuid.uuid4()).split("-")[0]
print(f"running with prefix {output_table_prefix}")

REF_VET_TABLE_COUNT = -1
client = None
default_config = None

EXTRACT_SAMPLE_TABLE = f"{output_table_prefix}_sample_names"

def get_partition_range(i):
  if i < 1 or i > REF_VET_TABLE_COUNT:
    raise ValueError(f"out of partition range")

  return { 'start': (i-1)*SAMPLES_PER_PARTITION + 1, 'end': i*SAMPLES_PER_PARTITION }

def get_samples_for_partition(sample_ids, i):
  return [ s for s in sample_ids if s >= get_partition_range(i)['start'] and s <= get_partition_range(i)['end'] ]

def split_lists(samples, n):
  return [samples[i * n:(i + 1) * n] for i in range((len(samples) + n - 1) // n )]

def load_sample_names(sample_names_to_extract, fq_temp_table_dataset):
  schema = [ bigquery.SchemaField("sample_name", "STRING", mode="REQUIRED") ]
  fq_sample_table = f"{fq_temp_table_dataset}.{EXTRACT_SAMPLE_TABLE}"

  job_labels = client._default_query_job_config.labels
  job_labels["gvs_query_name"] = "load-sample-names"

  job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.CSV, skip_leading_rows=0, schema=schema, labels=job_labels)

  with open(sample_names_to_extract, "rb") as source_file:
    job = client.load_table_from_file(source_file, fq_sample_table, job_config=job_config)

  job.result()  # Waits for the job to complete.

  # setting the TTL needs to be done as a second API call
  table = bigquery.Table(fq_sample_table, schema=schema)
  expiration = datetime.datetime.utcnow() + datetime.timedelta(hours=TEMP_TABLE_TTL_HOURS)
  table.expires = expiration
  table = client.update_table(table, ["expires"])

  return fq_sample_table

def get_all_sample_ids(fq_destination_table_samples):
  sql = f"select sample_id from `{fq_destination_table_samples}`"

  query_return = utils.execute_with_retry(client, "read cohort sample table", sql)
  JOBS.append({'job': query_return['job'], 'label': query_return['label']})
  sample_ids = [row.sample_id for row in list(query_return['results'])]
  sample_ids.sort()
  return sample_ids

def create_extract_samples_table(control_samples, fq_destination_table_samples, fq_sample_name_table, fq_sample_mapping_table):
  sql = f"CREATE OR REPLACE TABLE `{fq_destination_table_samples}` AS (" \
        f"SELECT m.sample_id, m.sample_name, m.is_loaded, m.withdrawn, m.is_control FROM `{fq_sample_name_table}` s JOIN `{fq_sample_mapping_table}` m ON (s.sample_name = m.sample_name) " \
        f"WHERE m.is_loaded is TRUE AND m.is_control = {control_samples})"

  print(sql)

  query_return = utils.execute_with_retry(client, "create extract sample table", sql)
  JOBS.append({'job': query_return['job'], 'label': query_return['label']})
  return query_return['results']

def get_table_count(fq_pet_vet_dataset):
  sql = f"SELECT MAX(CAST(SPLIT(table_name, '_')[OFFSET(1)] AS INT64)) max_table_number " \
        f"FROM `{fq_pet_vet_dataset}.INFORMATION_SCHEMA.TABLES` " \
        f"WHERE REGEXP_CONTAINS(lower(table_name), r'^(pet_[0-9]+)$') "
  query_return = utils.execute_with_retry(client, "get max table", sql)
  JOBS.append({'job': query_return['job'], 'label': query_return['label']})
  return int([row.max_table_number for row in list(query_return['results'])][0])

def create_final_extract_vet_table(fq_destination_table_vet_data):
  # first, create the table structure

  sql = f"""
        CREATE OR REPLACE TABLE `{fq_destination_table_vet_data}` 
        (
              location      INT64,
              sample_id	    INT64,
              ref           STRING,
              alt           STRING,
              call_GT       STRING,
              call_GQ       INT64,
              call_AD       STRING,
              AS_QUALapprox STRING,
              QUALapprox    STRING,
              call_PL       STRING	
        )
          PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 26000000000000, 6500000000))
          CLUSTER BY location
          {FINAL_TABLE_TTL}        
        """
  print(sql)
  query_return = utils.execute_with_retry(client, "create final export vet table", sql)
  JOBS.append({'job': query_return['job'], 'label': query_return['label']})

def create_final_extract_ref_table(fq_destination_table_ref_data):
  # first, create the table structure
  sql = f"""
        CREATE OR REPLACE TABLE `{fq_destination_table_ref_data}` 
        (
              location      INT64,
              sample_id	    INT64,
              length        INT64,
              state	        STRING	
        )
          PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 26000000000000, 6500000000))
          CLUSTER BY location
          {FINAL_TABLE_TTL}        
        """
  print(sql)
  query_return = utils.execute_with_retry(client, "create final export ref table", sql)
  JOBS.append({'job': query_return['job'], 'label': query_return['label']})


def populate_final_extract_table_with_ref(fq_ranges_dataset, fq_destination_table_data, sample_ids):
  def get_ref_subselect(fq_ref_table, samples, id):
    sample_stanza = ','.join([str(s) for s in samples])
    sql = f"    q_{id} AS (SELECT location, sample_id, length, state FROM \n" \
          f" `{fq_ref_table}` WHERE sample_id IN ({sample_stanza})), "
    return sql

  for i in range(1, REF_VET_TABLE_COUNT+1):
    partition_samples = get_samples_for_partition(sample_ids, i)  #sample ids for the partition

    if len(partition_samples) > 0:
      subs = {}
      insert = f"\nINSERT INTO `{fq_destination_table_data}` (location, sample_id, length, state) \n WITH \n"
      fq_ref_table = f"{fq_ranges_dataset}.{REF_TABLE_PREFIX}{i:03}"
      j = 1

      for samples in split_lists(partition_samples, 1000):
        id = f"{i}_{j}"
        subs[id] = get_ref_subselect(fq_ref_table, samples, id)
        j = j + 1

      sql = insert + ("\n".join(subs.values())) + "\n" + \
            "q_all AS (" + (" union all ".join([ f"(SELECT * FROM q_{id})" for id in subs.keys()]))  + ")\n" + \
            f" (SELECT * FROM q_all)"
      print(sql)
      print(f"{fq_ref_table} query is {utils.utf8len(sql)/(1024*1024)} MB in length")
      query_return = utils.execute_with_retry(client, "populate destination table with reference data", sql)
      JOBS.append({'job': query_return['job'], 'label': query_return['label']})
  return

def populate_final_extract_table_with_vet(fq_ranges_dataset, fq_destination_table_data, sample_ids):
  def get_ref_subselect(fq_vet_table, samples, id):
    sample_stanza = ','.join([str(s) for s in samples])
    sql = f"    q_{id} AS (SELECT location, sample_id, ref, alt, call_GT, call_GQ, call_AD, AS_QUALapprox, QUALapprox, CALL_PL FROM \n" \
          f" `{fq_vet_table}` WHERE sample_id IN ({sample_stanza})), "
    return sql

  for i in range(1, REF_VET_TABLE_COUNT+1):
    partition_samples = get_samples_for_partition(sample_ids, i)  #sample ids for the partition

    if len(partition_samples) > 0:
      subs = {}
      insert = f"\nINSERT INTO `{fq_destination_table_data}` (location, sample_id, ref, alt, call_GT, call_GQ, call_AD, AS_QUALapprox, QUALapprox, CALL_PL) \n WITH \n"
      fq_vet_table = f"{fq_ranges_dataset}.{VET_TABLE_PREFIX}{i:03}"
      j = 1

      for samples in split_lists(partition_samples, 1000):
        id = f"{i}_{j}"
        subs[id] = get_ref_subselect(fq_vet_table, samples, id)
        j = j + 1

      sql = insert + ("\n".join(subs.values())) + "\n" + \
            "q_all AS (" + (" union all ".join([ f"(SELECT * FROM q_{id})" for id in subs.keys()]))  + ")\n" + \
            f" (SELECT * FROM q_all)"
      print(sql)
      print(f"{fq_vet_table} query is {utils.utf8len(sql)/(1024*1024)} MB in length")
      query_return = utils.execute_with_retry(client, "populate destination table with variant data", sql)
      JOBS.append({'job': query_return['job'], 'label': query_return['label']})
  return

def make_extract_table(call_set_identifier,
                       control_samples,
                       fq_ranges_dataset,
                       max_tables,
                       sample_names_to_extract,
                       fq_cohort_sample_names,
                       query_project,
                       query_labels,
                       fq_temp_table_dataset,
                       fq_destination_dataset,
                       destination_table_prefix,
                       fq_sample_mapping_table,
                       sa_key_path,
                       temp_table_ttl_hours
                       ):
  try:
    fq_destination_table_ref_data = f"{fq_destination_dataset}.{destination_table_prefix}__REF_DATA"
    fq_destination_table_vet_data = f"{fq_destination_dataset}.{destination_table_prefix}__VET_DATA"
    fq_destination_table_samples = f"{fq_destination_dataset}.{destination_table_prefix}__SAMPLES"

    global client
    global default_config
    # this is where a set of labels are being created for the cohort extract
    query_labels_map = {}
    query_labels_map["id"]= output_table_prefix
    query_labels_map["gvs_tool_name"]= "gvs_prepare_ranges_callset"

    # query_labels is string that looks like 'key1=val1, key2=val2'
    if query_labels is not None and len(query_labels) != 0:
      for query_label in query_labels:
        kv = query_label.split("=", 2)
        key = kv[0].strip().lower()
        value = kv[1].strip().lower()
        query_labels_map[key] = value

        if not (bool(re.match(r"[a-z0-9_-]+$", key)) & bool(re.match(r"[a-z0-9_-]+$", value))):
          raise ValueError(f"label key or value did not pass validation--format should be 'key1=val1, key2=val2'")

    # add labels for DSP Cloud Cost Control Labeling and Reporting
    query_labels.update({'service':'gvs','team':'variants','managedby':'prepare_ranges_callset'})

    #Default QueryJobConfig will be merged into job configs passed in
    #but if a specific default config is being updated (eg labels), new config must be added
    #to the client._default_query_job_config that already exists
    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE", use_query_cache=True)

    if sa_key_path:
      credentials = service_account.Credentials.from_service_account_file(
        sa_key_path, scopes=["https://www.googleapis.com/auth/cloud-platform"],
      )

      client = bigquery.Client(credentials=credentials,
                               project=query_project,
                               default_query_job_config=default_config)
    else:
      client = bigquery.Client(project=query_project,
                               default_query_job_config=default_config)

    ## TODO -- provide a cmdline arg to override this (so we can simulate smaller datasets)

    global REF_VET_TABLE_COUNT
    REF_VET_TABLE_COUNT = max_tables

    global TEMP_TABLE_TTL_HOURS
    TEMP_TABLE_TTL_HOURS = temp_table_ttl_hours

    global TEMP_TABLE_TTL
    TEMP_TABLE_TTL = f" OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL {TEMP_TABLE_TTL_HOURS} HOUR)) "

    print(f"Using {REF_VET_TABLE_COUNT} tables in {fq_ranges_dataset}...")

    # if we have a file of sample names, load it into a temporary table
    if (sample_names_to_extract):
      fq_sample_name_table = load_sample_names(sample_names_to_extract, fq_temp_table_dataset)
    else:
      fq_sample_name_table = fq_cohort_sample_names

    # At this point one way or the other we have a table of sample names in BQ,
    # join it to the sample_info table to drive the extract
    create_extract_samples_table(control_samples, fq_destination_table_samples, fq_sample_name_table, fq_sample_mapping_table)

    # pull the sample ids back down
    sample_ids = get_all_sample_ids(fq_destination_table_samples)
    print(f"Discovered {len(sample_ids)} samples in {fq_destination_table_samples}...")

    # create the tables for extract data
    create_final_extract_ref_table(fq_destination_table_ref_data)
    create_final_extract_vet_table(fq_destination_table_vet_data)

    populate_final_extract_table_with_ref(fq_ranges_dataset, fq_destination_table_ref_data, sample_ids)
    populate_final_extract_table_with_vet(fq_ranges_dataset, fq_destination_table_vet_data, sample_ids)

  finally:
    utils.dump_job_stats(JOBS, client, f"{fq_destination_dataset}", call_set_identifier, 'GvsPrepareRanges', 'PrepareRangesCallsetTask', output_table_prefix)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract a cohort from BigQuery Variant Store ')
  parser.add_argument('--call_set_identifier',type=str, help='callset identifier used to track costs in cost_observability table', default='false')
  parser.add_argument('--control_samples',type=str, help='true for control samples only, false for participant samples only', default='false')
  parser.add_argument('--fq_ranges_dataset',type=str, help='project.dataset location of ranges/vet data', required=True)
  parser.add_argument('--fq_temp_table_dataset',type=str, help='project.dataset location where results should be stored', required=True)
  parser.add_argument('--fq_destination_dataset',type=str, help='project.dataset location where results should be stored', required=True)
  parser.add_argument('--destination_cohort_table_prefix',type=str, help='prefix used for destination cohort extract tables (e.g. my_fantastic_cohort)', required=True)
  parser.add_argument('--query_project',type=str, help='Google project where query should be executed', required=True)
  parser.add_argument('--query_labels',type=str, action='append', help='Labels to put on the BQ query that will show up in the billing. Ex: --query_labels key1=value1 --query_labels key2=value2', required=False)
  parser.add_argument('--fq_sample_mapping_table',type=str, help='Mapping table from sample_id to sample_name', required=True)
  parser.add_argument('--sa_key_path',type=str, help='Path to json key file for SA', required=False)
  parser.add_argument('--max_tables',type=int, help='Maximum number of PET/VET tables to consider', required=False, default=250)
  parser.add_argument('--ttl',type=int, help='Temp table TTL in hours', required=False, default=72)

  sample_args = parser.add_mutually_exclusive_group(required=True)
  sample_args.add_argument('--sample_names_to_extract',type=str, help='File containing list of samples to extract, 1 per line')
  sample_args.add_argument('--fq_cohort_sample_names',type=str, help='FQN of cohort table to extract, contains "sample_name" column')


  # Execute the parse_args() method
  args = parser.parse_args()

  make_extract_table(args.call_set_identifier,
                     args.control_samples,
                     args.fq_ranges_dataset,
                     args.max_tables,
                     args.sample_names_to_extract,
                     args.fq_cohort_sample_names,
                     args.query_project,
                     args.query_labels,
                     args.fq_temp_table_dataset,
                     args.fq_destination_dataset,
                     args.destination_cohort_table_prefix,
                     args.fq_sample_mapping_table,
                     args.sa_key_path,
                     args.ttl)
