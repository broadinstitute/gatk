from abc import abstractmethod, ABC

import sqlite3
from google.cloud.bigquery import Client


class DatabaseClient(ABC):
    def __init__(self, client):
        self.client = client
        super(DatabaseClient, self).__init__()

    @abstractmethod
    def execute(self, query: str):
        pass


class BigQueryDatabaseClient(DatabaseClient):
    """ If running locally, run the following commandline to authenticate yourself:
    gcloud auth application-default login
    """

    def __init__(self, client=None, credentials_file=None):
        if client is not None:
            super(BigQueryDatabaseClient, self).__init__(client)
        else:
            if credentials_file is not None:
                bigquery_client = Client.from_service_account_json(credentials_file)
            else:
                raise ValueError("BigQueryDatabaseClient requires a client or a credentials_file.")
            super(BigQueryDatabaseClient, self).__init__(bigquery_client)

    def execute(self, query: str):
        query_job = self.client.query(query)  # API request
        rows = query_job.result()  # Waits for query to finish
        return rows


class SqLiteDatabaseClient(DatabaseClient):
    def __init__(self, client=None, db_file=None):
        if client is not None:
            super(SqLiteDatabaseClient, self).__init__(client)
        else:
            if db_file is not None:
                super(SqLiteDatabaseClient, self).__init__(sqlite3.connect(db_file).cursor())
            else:
                raise ValueError("SqLiteDatabaseClient requires a client or a db_file.")

    def execute(self, query: str):
        return self.client.execute(query)


if '__main__' == __name__:
    credentials_file = '/Users/kyuksel/ml4h/bigquery-viewer-credentials.json'
    db_client = BigQueryDatabaseClient(credentials_file=credentials_file)

    dataset = 'broad-ml4cvd.ukbb7089_r10data'

    dictionary_table = f"`{dataset}.dictionary`"
    phenotype_table = f"`{dataset}.phenotype`"
    coding_table = f"`{dataset}.coding`"

    fid = 20001
    fids = [3143, 3144]
    sample_id = 2907043

    job_title_field_id = 22600
    icd10_field = 41202
    query = \
        f"SELECT value FROM {phenotype_table} WHERE fieldid={icd10_field} AND sample_id={sample_id}"

    rows = db_client.execute(query.format())
    for row in rows:
        print(row)
