import argparse
from azure.cosmos import CosmosClient, DatabaseProxy
from azure.identity import DefaultAzureCredential


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Say Hello to Cosmos DB from Python')
    parser.add_argument('--endpoint', type=str, help='Cosmos DB NoSQL Endpoint URI', required=True)
    parser.add_argument('--database', type=str, help='Cosmos DB database name', required=True)
    args = parser.parse_args()

    credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)
    client = CosmosClient(url=args.endpoint, credential=credential)

    database_proxy = DatabaseProxy(client.client_connection, args.database)

    container_client = database_proxy.get_container_client('sample_info')

    for container in database_proxy.list_containers():
        print(container['id'])
