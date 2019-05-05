""" Creates a ukbb_dev dataset which is a miniature version of a typical ukbbs.

NOTE: dataset names, tables to be copied are currently hardcoded,
will fail if dev dataset already
exists.
"""
from google.cloud import bigquery

# CONSTANTS

ORIGINAL_DATASET = 'ukbb7089_201904'
DEV_DATASET = 'ukbb_dev'
EXACT_TABLES = ['censor', 'coding', 'dictionary']
LIMITED_SAMPLE_TABLES = ['hesin', 'hesin_diag10', 'hesin_diag9', 'hesin_oper']
client = bigquery.Client()  # should already be set to default project


if __name__ == '__main__':
    PROJECT = client.project
    FULL_DEV_DATASET = f"{PROJECT}.{DEV_DATASET}"
    print(FULL_DEV_DATASET)
    dataset = bigquery.Dataset.from_string(FULL_DEV_DATASET)
    dataset = client.create_dataset(dataset)
    print('Dataset {} created.'.format(dataset.dataset_id))

    print(f"working on {DEV_DATASET}.phenotype")
    # create 1/1000th size dataset
    query_job = client.query(f"""
                CREATE TABLE {DEV_DATASET}.phenotype
                AS SELECT * FROM {ORIGINAL_DATASET}.phenotype
                WHERE MOD(sample_id,1000)=4""")

    rows = query_job.result()

    # copy some tables exactly
    for table in EXACT_TABLES:
        print(f"working on {DEV_DATASET}.{table}")
        query_job = client.query(f"""
                CREATE TABLE {DEV_DATASET}.{table}
                AS SELECT * FROM {ORIGINAL_DATASET}.{table}
                """)
        rows = query_job.result()

    # copy hesin tables by limiting eids within sample_ids
    for table in LIMITED_SAMPLE_TABLES:
        print(f"working on {DEV_DATASET}.{table}")
        query_job = client.query(f"""
                CREATE TABLE {DEV_DATASET}.{table}
                AS SELECT * FROM {ORIGINAL_DATASET}.{table}
                where eid in (select sample_id from
                {DEV_DATASET}.phenotype)
                """)
        rows = query_job.result()
