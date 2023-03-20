from azure.cosmos import CosmosClient, DatabaseProxy
from azure.identity import DefaultAzureCredential

import argparse
import pandas


"""
Python-based Cosmos data loader, good enough to load 10 sample_info records but definitely not recommended for loading
"real" variants data in the tens of millions of rows per sample.
"""


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Say Hello to Cosmos DB from Python')
    parser.add_argument('--endpoint', type=str, help='Cosmos DB NoSQL Endpoint URI', required=True)
    parser.add_argument('--database', type=str, help='Cosmos DB database name', required=True)
    parser.add_argument('--data-table-tsv', type=str,
                        help='TSV representing the data table with sample, VCF, and index data', required=True)
    parser.add_argument('--data-table-sample-name', type=str, help='Header representing the sample name', required=True)

    args = parser.parse_args()

    credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)
    client = CosmosClient(url=args.endpoint, credential=credential)

    database_proxy = DatabaseProxy(client.client_connection, args.database)

    container_client = database_proxy.get_container_client('sample_info')

    dataframe = pandas.read_csv(args.data_table_tsv, sep='\t')
    item_id = 0
    for idx, sample_name in enumerate(dataframe[args.data_table_sample_name]):
        item_id = idx + 1
        item = {
            "id": sample_name,
            "sampleId": item_id,
            "sampleName": sample_name,
            "isLoaded": False,
            "isControl": False,
            "withdrawn": False
        }
        container_client.create_item(item)
        if item_id % 10 == 0:
            print(f"{item_id} items created")

    if item_id % 10 != 0:
        print(f"{item_id} items created")
