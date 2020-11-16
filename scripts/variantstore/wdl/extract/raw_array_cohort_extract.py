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
RAW_ARRAY_TABLE_PREFIX = "arrays_"
SAMPLES_PER_PARTITION = 4000

RAW_ARRAY_TABLE_COUNT = -1
client = None

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
  while True:
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
        print(f"Error {err} running query {label}, but no retries left.  RAISING!")
        raise err
      else:
        t = retry_delay.pop(0)
        print(f"Error {err} running query {label}, sleeping for {t}")
        time.sleep(t)
  
def get_partition_range(i):
  if i < 1 or i > RAW_ARRAY_TABLE_COUNT:
    raise ValueError(f"out of partition range")

  return { 'start': (i-1)*SAMPLES_PER_PARTITION + 1, 'end': i*SAMPLES_PER_PARTITION }

def get_samples_for_partition(cohort, i):
  return [ s for s in cohort if s >= get_partition_range(i)['start'] and s <= get_partition_range(i)['end'] ]

def split_lists(samples, n):
  return [samples[i * n:(i + 1) * n] for i in range((len(samples) + n - 1) // n )]

def get_all_samples(fq_sample_mapping_table, cohort_sample_names_file, sample_map_outfile):
  sample_names = [line.strip() for line in open(cohort_sample_names_file).readlines()]
  joined_sample_names = ",".join('"' + s + '"' for s in sample_names)

  sql = f"select sample_id,sample_name from `{fq_sample_mapping_table}` WHERE sample_name IN ({joined_sample_names})"
      
  results = execute_with_retry("read cohort table", sql)    

  cohort = []
  csv_str = ""

  for row in list(results):
    cohort.append(row.sample_id)
    csv_str = csv_str + str(row.sample_id) + "," + row.sample_name + "\n"
  
  cohort.sort()

  if sample_map_outfile is not None:
    with open(sample_map_outfile, 'w') as outfile:
      outfile.write(csv_str)

  return cohort

def populate_extract_table(fq_dataset, cohort, fq_destination_table, ttl, number_of_partitions, probes_per_partition, extract_genotype_counts_only):
  def get_subselect(fq_array_table, samples, id, extract_genotype_counts_only):
    fields_to_extract = "sample_id, probe_id, GT_encoded" if extract_genotype_counts_only else "sample_id, probe_id, GT_encoded, NORMX, NORMY, BAF, LRR"
    sample_stanza = ','.join([str(s) for s in samples])
    sql = f"    q_{id} AS (SELECT {fields_to_extract} from `{fq_array_table}` WHERE sample_id IN ({sample_stanza})), "
    return sql
   
  subs = {}
  for i in range(1, RAW_ARRAY_TABLE_COUNT+1):
    partition_samples = get_samples_for_partition(cohort, i)

    fq_array_table = f"{fq_dataset}.{RAW_ARRAY_TABLE_PREFIX}{i:03}"
    if len(partition_samples) > 0:
      j = 1
      # subset the query to 1000 samples because if you have more than that BigQuery doesn't optimize the query correctly
      # TODO: test the 1000 items in a list limit, maybe BigQuery has fixed this by now
      for samples in split_lists(partition_samples, 1000):
        id = f"{i}_{j}"
        subs[id] = get_subselect(fq_array_table, samples, id, extract_genotype_counts_only)
        j = j + 1

  # TODO: make genotype flexible to older encodings
  select_sql = (
                f" (SELECT probe_id, " +
                f"SUM(IF(GT_encoded is null OR GT_encoded = 'R', 1, 0)) hom_ref, \n" +
                f"SUM(IF(GT_encoded = 'X', 1, 0)) het, \n" +
                f"SUM(IF(GT_encoded = 'A', 1, 0)) hom_var, \n" +
                f"SUM(IF(GT_encoded = 'Y', 1, 0)) het_1_2, \n" +
                f"SUM(IF(GT_encoded = 'B', 1, 0)) hom_var_2_2, \n" +
                f"SUM(IF(GT_encoded = '.', 1, 0)) no_call \n" +
                f"FROM q_all \n" +
                f"GROUP BY probe_id)"
  ) if extract_genotype_counts_only else f" (SELECT * FROM q_all)"

  sql = (
        f"CREATE OR REPLACE TABLE `{fq_destination_table}` \n"
        f"PARTITION BY RANGE_BUCKET(probe_id, GENERATE_ARRAY(0, {number_of_partitions * probes_per_partition}, {probes_per_partition})) \n"
        f"OPTIONS( expiration_timestamp=TIMESTAMP_ADD(CURRENT_TIMESTAMP(), INTERVAL {ttl} HOUR)) "
        f"AS \n" +
        f"with\n" +
        ("\n".join(subs.values())) + "\n"
        "q_all AS (" + (" union all ".join([ f"(SELECT * FROM q_{id})" for id in subs.keys()])) + ")\n" +
        f"{select_sql}"
        )

  print(sql) 
  print(f"Extract Query is {utf8len(sql)/(1024*1024)} MB in length")  
  results = execute_with_retry("create extract table table", sql)    
  return results

def do_extract(fq_dataset,
               max_tables,
               query_project,
               fq_destination_table,
               fq_sample_mapping_table,
               cohort_sample_names_file,
               sample_map_outfile,
               ttl,
               number_of_partitions,
               probes_per_partition,
               extract_genotype_counts_only
              ):
  try:  
    global client
    client = bigquery.Client(project=query_project, 
                             default_query_job_config=QueryJobConfig(priority="INTERACTIVE", use_query_cache=False ))

    global RAW_ARRAY_TABLE_COUNT
    RAW_ARRAY_TABLE_COUNT = max_tables
    print(f"Using {RAW_ARRAY_TABLE_COUNT} tables in {fq_dataset}...")

    cohort = get_all_samples(fq_sample_mapping_table, cohort_sample_names_file, sample_map_outfile)
    print(f"Discovered {len(cohort)} samples in {fq_sample_mapping_table}...")

    populate_extract_table(fq_dataset, cohort, fq_destination_table, ttl, number_of_partitions, probes_per_partition, extract_genotype_counts_only)

    print(f"\nFinal cohort extract written to {fq_destination_table}\n")
  except Exception as err:
      print(f"Unexpected error! {err}")
      raise
  finally:
      dump_job_stats()

if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Extract a raw array cohort from BigQuery Variant Store ')
  
  parser.add_argument('--dataset',type=str, help='project.dataset location of raw array data', required=True)
  parser.add_argument('--fq_destination_table',type=str, help='fully qualified destination table', required=True)
  parser.add_argument('--query_project',type=str, help='Google project where query should be executed', required=True)
  parser.add_argument('--fq_sample_mapping_table',type=str, help='Mapping table from sample_id to sample_name', required=True)
  parser.add_argument('--cohort_sample_names_file',type=str, help='File containing newline separated sample_names in the cohort', required=True)
  parser.add_argument('--max_tables',type=int, help='Maximum number of array_xxx tables to consider', required=False, default=250)
  parser.add_argument('--ttl',type=int, help='how long should the destination table be kept before expiring (in hours)', required=False, default=24)
  parser.add_argument('--number_of_partitions',type=int, help='how many partitions to create', required=False, default=1)
  parser.add_argument('--probes_per_partition',type=int, help='how many probes in each partition', required=False, default=2000000)
  parser.add_argument('--extract_genotype_counts_only', type=bool, help='Extract only genoype counts for QC metric calculations', required=False, default=False)
  parser.add_argument('--sample_map_outfile', type=str, help='Write out sample_id,sample_name map as a CSV', required=False, default=False)

  # Execute the parse_args() method
  args = parser.parse_args()

  do_extract(args.dataset,
             args.max_tables,
             args.query_project,
             args.fq_destination_table,
             args.fq_sample_mapping_table,
             args.cohort_sample_names_file,
             args.sample_map_outfile,
             args.ttl,
             args.number_of_partitions, 
             args.probes_per_partition,
             args.extract_genotype_counts_only)
             
