import argparse
from azure.cosmos import CosmosClient
from azure.identity import DefaultAzureCredential


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Say Hello to Cosmos DB from Python')
    parser.add_argument('--cosmos-endpoint', type=str, help='Cosmos DB NoSQL Endpoint URI', required=True)
    args = parser.parse_args()

    credential = DefaultAzureCredential(exclude_shared_token_cache_credential=True)
    client = CosmosClient(url=args.cosmos_endpoint, credential=credential)

    print(client)
