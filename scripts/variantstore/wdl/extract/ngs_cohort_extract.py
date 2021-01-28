# -*- coding: utf-8 -*-
import sys
import uuid
import time

from concurrent.futures import ThreadPoolExecutor, as_completed
from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig

import argparse

JOB_IDS = set()

#
# CONSTANTS
#
PET_TABLE_PREFIX = "pet_"
VET_TABLE_PREFIX = "vet_"
SAMPLES_PER_PARTITION = 4000

TEMP_TABLE_TTL = " OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 1 HOUR)) "
FINAL_TABLE_TTL = ""
#FINAL_TABLE_TTL = " OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL 72 HOUR)) "

# temp-table-uuid
output_table_prefix = str(uuid.uuid4()).split("-")[0]
print(f"running with prefix {output_table_prefix}")

PET_VET_TABLE_COUNT = -1
client = None

VET_DISTINCT_POS_TABLE = f"{output_table_prefix}_vet_distinct_pos"
PET_NEW_TABLE = f"{output_table_prefix}_pet_new"
VET_NEW_TABLE = f"{output_table_prefix}_vet_new"
COHORT_EXTRACT_TABLE = f"{output_table_prefix}_cohort_extract"

def utf8len(s):
    return len(s.encode('utf-8'))

def dump_job_stats():
  total = 0

  for jobid in JOB_IDS:
    job = client.get_job(jobid[1])

    bytes_billed = int(0 if job.total_bytes_billed is None else job.total_bytes_billed)
    total = total + bytes_billed

    print(jobid[0], " <====> Cache Hit:", job.cache_hit, bytes_billed/(1024 * 1024), " MBs")

  print(" Total GBs billed ", total/(1024 * 1024 * 1024), " GBs")

def execute_with_retry(label, sql):
  retry_delay = [30, 60, 90] # 3 retries with incremental backoff

  start = time.time()
  while len(retry_delay) > 0:
    try:
      query = client.query(sql)
      
      print(f"STARTING - {label}")
      JOB_IDS.add((label, query.job_id))
      results = query.result()
      print(f"COMPLETED ({time.time() - start} s, {3-len(retry_delay)} retries) - {label}")
      return results
    except Exception as err:

      # if there are no retries left... raise
      if (len(retry_delay) == 0):
        raise err
      else:
        t = retry_delay.pop(0)
        print(f"Error {err} running query {label}, sleeping for {t}")
        time.sleep(t)

def get_partition_range(i):
  if i < 1 or i > PET_VET_TABLE_COUNT:
    raise ValueError(f"out of partition range")

  return { 'start': (i-1)*SAMPLES_PER_PARTITION + 1, 'end': i*SAMPLES_PER_PARTITION }

def get_samples_for_partition(cohort, i):
  return [ s for s in cohort if s >= get_partition_range(i)['start'] and s <= get_partition_range(i)['end'] ]

def split_lists(samples, n):
  return [samples[i * n:(i + 1) * n] for i in range((len(samples) + n - 1) // n )]

def get_all_samples(fq_cohort_sample_names, fq_sample_mapping_table):
  sql = f"select m.sample_id from `{fq_cohort_sample_names}` c JOIN `{fq_sample_mapping_table}` m ON (m.sample_name = c.sample_name)"
      
  results = execute_with_retry("read cohort table", sql)    
  cohort = [row.sample_id for row in list(results)]
  cohort.sort()
  return cohort

def get_table_count(fq_pet_vet_dataset):
  sql = f"SELECT MAX(CAST(SPLIT(table_name, '_')[OFFSET(1)] AS INT64)) max_table_number " \
        f"FROM `{fq_pet_vet_dataset}.INFORMATION_SCHEMA.TABLES` " \
        f"WHERE REGEXP_CONTAINS(lower(table_name), r'^(pet_[0-9]+)$') "
  results = execute_with_retry("get max table", sql)
  return int([row.max_table_number for row in list(results)][0])

def make_new_vet_union_all(fq_pet_vet_dataset, fq_temp_table_dataset, cohort):
  def get_subselect(fq_vet_table, samples, id):
    sample_stanza = ','.join([str(s) for s in samples])
    sql = f"    q_{id} AS (SELECT location, sample_id, ref, alt, call_GT, call_GQ, call_pl, AS_QUALapprox from `{fq_vet_table}` WHERE sample_id IN ({sample_stanza})), "
    return sql
   
  subs = {}
  for i in range(1, PET_VET_TABLE_COUNT+1):
    partition_samples = get_samples_for_partition(cohort, i)

      # KCIBUL -- grr, should be fixed width
    fq_vet_table = f"{fq_pet_vet_dataset}.{VET_TABLE_PREFIX}{i:03}"
    if len(partition_samples) > 0:
      j = 1
      for samples in split_lists(partition_samples, 1000):
        id = f"{i}_{j}"
        subs[id] = get_subselect(fq_vet_table, samples, id)
        j = j + 1

  sql = f"create or replace table `{fq_temp_table_dataset}.{VET_NEW_TABLE}` {TEMP_TABLE_TTL} AS \n" + \
        "with\n" + \
        ("\n".join(subs.values())) + "\n" \
        "q_all AS (" + (" union all ".join([ f"(SELECT * FROM q_{id})" for id in subs.keys()]))  + ")\n" + \
        f" (SELECT * FROM q_all)"

  print(sql) 
  print(f"VET Query is {utf8len(sql)/(1024*1024)} MB in length")  
  results = execute_with_retry("insert vet new table", sql)    
  return results



def create_position_table(fq_temp_table_dataset, min_variant_samples):
  dest = f"{fq_temp_table_dataset}.{VET_DISTINCT_POS_TABLE}"
  
  # only create this clause if min_variant_samples > 0, becuase if 
  # it is == 0 then we don't need to touch the sample_id column (which doubles the cost of this query)
  min_sample_clause = ""
  if min_variant_samples > 0:
      min_sample_clause = f"HAVING count(distinct sample_id) >= {min_variant_samples}"

  create_vet_distinct_pos_query = client.query(
        f"""
          create or replace table `{dest}` {TEMP_TABLE_TTL}
          as (
            select location from `{fq_temp_table_dataset}.{VET_NEW_TABLE}` WHERE alt != '*' GROUP BY location {min_sample_clause}
          )
        """
    )

  create_vet_distinct_pos_query.result()
  JOB_IDS.add((f"create positions table {dest}", create_vet_distinct_pos_query.job_id))
  return

def make_new_pet_union_all(fq_pet_vet_dataset, fq_temp_table_dataset, cohort):
  def get_pet_subselect(fq_pet_table, samples, id):
    sample_stanza = ','.join([str(s) for s in samples])
    sql = f"    q_{id} AS (SELECT p.location, p.sample_id, p.state from {fq_pet_table} p " \
          f"    join `{fq_temp_table_dataset}.{VET_DISTINCT_POS_TABLE}` v on (p.location = v.location) WHERE p.sample_id IN ({sample_stanza})), "
    return sql

  subs = {}
  for i in range(1, PET_VET_TABLE_COUNT+1):
    partition_samples = get_samples_for_partition(cohort, i)

      # KCIBUL -- grr, should be fixed width
    fq_pet_table = f"{fq_pet_vet_dataset}.{PET_TABLE_PREFIX}{i:03}"
    if len(partition_samples) > 0:
      j = 1
      for samples in split_lists(partition_samples, 1000):
        id = f"{i}_{j}"
        subs[id] = get_pet_subselect(fq_pet_table, samples, id)
        j = j + 1

  sql = f"create or replace table `{fq_temp_table_dataset}.{PET_NEW_TABLE}` {TEMP_TABLE_TTL} AS \n" + \
        "with\n" + \
        ("\n".join(subs.values())) + "\n" \
        "q_all AS (" + (" union all ".join([ f"(SELECT * FROM q_{id})" for id in subs.keys()]))  + ")\n" + \
        f" (SELECT * FROM q_all)"

  #print(sql)      
  print(f"PET Query is {utf8len(sql)/(1024*1024)} MB in length")
  results = execute_with_retry("insert pet new table", sql)    
  return results


def populate_final_extract_table(fq_temp_table_dataset, fq_destination_dataset, destination_table, fq_sample_mapping_table):
  dest = f"{fq_destination_dataset}.{destination_table}"
  sql = f"""
        CREATE OR REPLACE TABLE `{dest}` 
        PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 26000000000000, 6500000000))
        CLUSTER BY location
        {FINAL_TABLE_TTL}
        as (SELECT
            new_pet.location,
            s.sample_name as sample_name,
            new_pet.state,
            new_vet.ref,
            REPLACE(new_vet.alt,",<NON_REF>","") alt,
            new_vet.call_GT,
            new_vet.call_GQ,
            cast(SPLIT(new_vet.call_pl,",")[OFFSET(0)] as int64) as call_RGQ,
            new_vet.AS_QUALapprox,
            new_vet.call_PL
          FROM
            `{fq_temp_table_dataset}.{PET_NEW_TABLE}` new_pet
          LEFT OUTER JOIN
            `{fq_temp_table_dataset}.{VET_NEW_TABLE}`  new_vet
          ON (new_pet.location = new_vet.location AND new_pet.sample_id = new_vet.sample_id)
          LEFT OUTER JOIN
            `{fq_sample_mapping_table}` s ON (new_pet.sample_id = s.sample_id))
      """
  #print(sql)
  cohort_extract_final_query_job = client.query(sql)

  cohort_extract_final_query_job.result()
  JOB_IDS.add((f"insert final cohort table {dest}", cohort_extract_final_query_job.job_id))  
  return

def do_extract(fq_pet_vet_dataset,
               max_tables,
               fq_cohort_sample_names,
               query_project,
               fq_temp_table_dataset,
               fq_destination_dataset,
               destination_table,
               min_variant_samples,
               fq_sample_mapping_table
              ):
  try:  
      
    global client
    client = bigquery.Client(project=query_project, 
                             default_query_job_config=QueryJobConfig(labels={ "id" : f"test_cohort_export_{output_table_prefix}"}, priority="INTERACTIVE", use_query_cache=False ))

    ## TODO -- provide a cmdline arg to override this (so we can simulat smaller datasets)
    global PET_VET_TABLE_COUNT
    PET_VET_TABLE_COUNT = max_tables
    print(f"Using {PET_VET_TABLE_COUNT} PET tables in {fq_pet_vet_dataset}...")

    cohort = get_all_samples(fq_cohort_sample_names, fq_sample_mapping_table)
    print(f"Discovered {len(cohort)} samples in {fq_cohort_sample_names}...")

    make_new_vet_union_all(fq_pet_vet_dataset, fq_temp_table_dataset, cohort)

    create_position_table(fq_temp_table_dataset, min_variant_samples)  
    make_new_pet_union_all(fq_pet_vet_dataset, fq_temp_table_dataset, cohort)
    populate_final_extract_table(fq_temp_table_dataset, fq_destination_dataset, destination_table, fq_sample_mapping_table)
  except Exception as err:
    print(err)

  dump_job_stats()
  print(f"\nFinal cohort extract written to {fq_destination_dataset}.{destination_table}\n")

if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract a cohort from BigQuery Variant Store ')
  
  parser.add_argument('--fq_petvet_dataset',type=str, help='project.dataset location of pet/vet data', required=True)
  parser.add_argument('--fq_temp_table_dataset',type=str, help='project.dataset location where results should be stored', required=True)
  parser.add_argument('--fq_destination_dataset',type=str, help='project.dataset location where results should be stored', required=True)
  parser.add_argument('--destination_table',type=str, help='destination table', required=True)
  parser.add_argument('--fq_cohort_sample_names',type=str, help='FQN of cohort table to extract, contains "sample_name" column', required=True)
  parser.add_argument('--query_project',type=str, help='Google project where query should be executed', required=True)
  parser.add_argument('--min_variant_samples',type=int, help='Minimum variant samples at a site required to be emitted', required=False, default=0)
  parser.add_argument('--fq_sample_mapping_table',type=str, help='Mapping table from sample_id to sample_name', required=True)

  parser.add_argument('--max_tables',type=int, help='Maximum number of PET/VET tables to consider', required=False, default=250)


  # Execute the parse_args() method
  args = parser.parse_args()

  do_extract(args.fq_petvet_dataset,
             args.max_tables,
             args.fq_cohort_sample_names,
             args.query_project,
             args.fq_temp_table_dataset,
             args.fq_destination_dataset,
             args.destination_table,
             args.min_variant_samples,
             args.fq_sample_mapping_table)