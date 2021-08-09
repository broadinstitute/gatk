from google.cloud import bigquery
from google.cloud.bigquery.job import QueryJobConfig
from google.oauth2 import service_account

import argparse

client = None
THREE_QUOTES = "\"\"\""
alt_allele_temp_function = r"""
                            CREATE TEMPORARY FUNCTION minimize(ref STRING, allele STRING)
                            RETURNS STRING 
                              LANGUAGE js AS """ + THREE_QUOTES + r"""
                                let done = false
                                while (!done && ref.length !== 1) {
                                  if (ref.slice(-1) === allele.slice(-1)) {
                                    ref = ref.slice(0, -1)
                                    allele = allele.slice(0,-1)
                                  } else {
                                    done = true
                                  }
                                }
                                return ref+','+allele
                            """ + THREE_QUOTES + ";"
alt_allele_query_end = r"""
select location, sample_id, 
SPLIT(minimize(ref, SPLIT(alt,',')[OFFSET(0)]))[OFFSET(0)] as ref,
SPLIT(minimize(ref, SPLIT(alt,',')[OFFSET(0)]))[OFFSET(1)] as allele,
1 as allele_pos, call_GT, call_GQ,
as_raw_mq,
cast(SPLIT(as_raw_mq,'|')[OFFSET(1)] as int64) raw_mq,
as_raw_mqranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_mqranksum,',')[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_mqranksum_x_10, 
as_qualapprox,
qualapprox,
cast(SPLIT(as_qualapprox,'|')[OFFSET(0)] as int64) as qual,
as_raw_readposranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_readposranksum,',')[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_readposranksum_x_10,
as_sb_table,
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(0)],',')[OFFSET(0)] as int64) as sb_ref_plus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(0)],',')[OFFSET(1)] as int64) as sb_ref_minus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(1)],',')[OFFSET(0)] as int64) as sb_alt_plus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(1)],',')[OFFSET(1)] as int64) as sb_alt_minus, 
call_AD,
cast(SPLIT(call_AD,',')[OFFSET(0)] as int64) as ref_ad, 
cast(SPLIT(call_AD,',')[OFFSET(1)] as int64) as ad
from position1

union all

select location, sample_id, 
SPLIT(minimize(ref, SPLIT(alt,',')[OFFSET(0)]))[OFFSET(0)] as ref,
SPLIT(minimize(ref, SPLIT(alt,',')[OFFSET(0)]))[OFFSET(1)] as allele,
1 as allele_pos, call_GT, call_GQ,
as_raw_mq,
cast(SPLIT(as_raw_mq,'|')[OFFSET(1)] as int64) raw_mq,
as_raw_mqranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_mqranksum,',')[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_mqranksum_x_10,
as_qualapprox,
qualapprox,
cast(SPLIT(as_qualapprox,'|')[OFFSET(0)] as int64) as qual,
as_raw_readposranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_readposranksum,',')[SAFE_OFFSET(0)] as float64) * 10.0 as int64) as raw_readposranksum_x_10,
as_sb_table,
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(0)],',')[OFFSET(0)] as int64) as sb_ref_plus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(0)],',')[OFFSET(1)] as int64) as sb_ref_minus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(1)],',')[OFFSET(0)] as int64) as sb_alt_plus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(1)],',')[OFFSET(1)] as int64) as sb_alt_minus, 
call_AD,
cast(SPLIT(call_AD,',')[OFFSET(0)] as int64) as ref_ad, 
cast(SPLIT(call_AD,',')[OFFSET(1)] as int64) as ad
from position2

union all

select location, sample_id,
SPLIT(minimize(ref, SPLIT(alt,',')[OFFSET(1)]))[OFFSET(0)] as ref,
SPLIT(minimize(ref, SPLIT(alt,',')[OFFSET(1)]))[OFFSET(1)] as allele,
2 as allele_pos, call_GT, call_GQ,
as_raw_mq,
cast(SPLIT(as_raw_mq,'|')[OFFSET(2)] as int64) raw_mq,
as_raw_mqranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_mqranksum,',')[SAFE_OFFSET(1)] as float64) * 10.0 as int64) as raw_mqranksum_x_10,
as_qualapprox,
qualapprox,
cast(SPLIT(as_qualapprox,'|')[OFFSET(1)] as int64) as qual,
as_raw_readposranksum,
SAFE_cast(SAFE_cast(SPLIT(as_raw_readposranksum,',')[SAFE_OFFSET(1)] as float64) * 10.0 as int64) as raw_readposranksum_x_10,
as_sb_table,
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(0)],',')[OFFSET(0)] as int64) as sb_ref_plus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(0)],',')[OFFSET(1)] as int64) as sb_ref_minus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(2)],',')[OFFSET(0)] as int64) as sb_alt_plus, 
cast(SPLIT(SPLIT(as_sb_table,'|')[OFFSET(2)],',')[OFFSET(1)] as int64) as sb_alt_minus, 
call_AD,
cast(SPLIT(call_AD,',')[OFFSET(0)] as int64) as ref_ad, 
cast(SPLIT(call_AD,',')[OFFSET(2)] as int64) as ad
from position2;
"""

def get_vet_table_names(fq_pet_vet_dataset):
  sql = f"SELECT table_name FROM {fq_pet_vet_dataset}.INFORMATION_SCHEMA.TABLES WHERE table_name LIKE 'vet_%' ORDER BY table_name"

  query = client.query(sql)
  results = query.result()
  return results


def make_alt_allele_table(fq_pet_vet_dataset,
                          query_project,
                          sa_key_path
                          ):
  global client
  default_config = QueryJobConfig(priority="INTERACTIVE", use_query_cache=True)

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

  vet_table_names = get_vet_table_names(fq_pet_vet_dataset)

  first = True
  for vet_table in vet_table_names:
    query_beginning = f"CREATE OR REPLACE TABLE {fq_pet_vet_dataset}.alt_allele \
                        PARTITION BY \
                        RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000, 1000000000000)) \
                        CLUSTER BY location, sample_id \
                        AS \n"
    if (not first):
      query_beginning = f"INSERT INTO {fq_pet_vet_dataset}.alt_allele \m"
      first = False
    fq_vet_table = f"{fq_pet_vet_dataset}.{vet_table.table_name}"
    query_with = f"""WITH 
                      position1 as (select * from {fq_vet_table} 
                       WHERE call_GT IN ('0/1', '1/0', '1/1', '0|1', '1|0', '1|1', '0/2', '0|2','2/0', '2|0')), 
                      position2 as (select * from {fq_vet_table} 
                       WHERE call_GT IN ('1/2', '1|2', '2/1', '2|1'))
                  """

    sql = alt_allele_temp_function + query_beginning + query_with + alt_allele_query_end
    query = client.query(sql)
    results = query.result()
    print(f"data loaded into alt_allele from {vet_table.table_name}\n")


if __name__ == '__main__':
  parser = argparse.ArgumentParser(allow_abbrev=False, description='Generate an alt_allele table for the BigQuery Variant Store ')

  parser.add_argument('--fq_petvet_dataset',type=str, help='project.dataset location of pet/vet data', required=True)
  parser.add_argument('--query_project',type=str, help='Google project where query should be executed', required=True)
  parser.add_argument('--sa_key_path',type=str, help='Path to json key file for SA', required=False)


  # Execute the parse_args() method
  args = parser.parse_args()

  make_alt_allele_table(args.fq_petvet_dataset,
                     args.query_project,
                     args.sa_key_path)
