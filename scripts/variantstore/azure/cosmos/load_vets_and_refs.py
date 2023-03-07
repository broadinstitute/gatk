from azure.cosmos import CosmosClient, DatabaseProxy
from azure.identity import DefaultAzureCredential
from fastavro import reader

import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Load basic Variants data into Cosmos')
    parser.add_argument('--endpoint', type=str, help='Cosmos DB NoSQL Endpoint URI', required=True)
    parser.add_argument('--database', type=str, help='Cosmos DB database name', required=True)
    parser.add_argument('--datatype', type=str, help='Data type to load, must be "vets" or "refs"', required=True)

    args, files = parser.parse_known_args()

    container_name = None
    if args.datatype == 'vets':
        container_name = 'vets'
    elif args.datatype == 'refs':
        container_name = 'ref_ranges'
    else:
        raise ValueError(f'--datatype must be "refs" or "vets"')

    credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)
    client = CosmosClient(url=args.endpoint, credential=credential)

    database_proxy = DatabaseProxy(client.client_connection, args.database)
    container_client = database_proxy.get_container_client(container_name)

    for file in files:
        print(f"Processing '{file}'...")
        with open(file, 'rb') as fo:
            avro_reader = reader(fo)
            count = 0
            for record in avro_reader:
                item = {'id': str(count)}
                for k, v in record.items():
                    item[k] = v

                container_client.create_item(item)
                count = count + 1

                if count % 10 == 0:
                    print(f"{count} items created")