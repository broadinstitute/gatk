from azure.cosmos import CosmosClient, DatabaseProxy
from azure.identity import DefaultAzureCredential

import argparse
import pandas

if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Say Hello to Cosmos DB from Python')
    parser.add_argument('--endpoint', type=str, help='Cosmos DB NoSQL Endpoint URI', required=True)
    parser.add_argument('--database', type=str, help='Cosmos DB database name', required=True)
    parser.add_argument('--data-table-tsv', type=str,
                        help='TSV representing the data table with sample, VCF, and index data', required=True)
    parser.add_argument('--data-table-sample-name', type=str, help='Header representing the sample name', required=True)
    # parser.add_argument('--data-table-vcf', type=str, help='Header representing the VCF path')
    # parser.add_argument('--data-table-vcf-index', type=str, help='Header representing the VCF index')

    args = parser.parse_args()

    credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)
    client = CosmosClient(url=args.endpoint, credential=credential)

    database_proxy = DatabaseProxy(client.client_connection, args.database)

    container_client = database_proxy.get_container_client('sample_info')

    dataframe = pandas.read_csv(args.data_table_tsv, sep='\t')
    for idx, sample_name in enumerate(dataframe[args.data_table_sample_name]):
        id = idx + 1
        item = {
            "id": sample_name,
            "sampleId": id,
            "sampleName": sample_name,
            "isLoaded": False,
            "isControl": False,
            "withdrawn": False
        }
        container_client.create_item(item)
        if id % 10 == 0:
            print(f"{id} items created")
