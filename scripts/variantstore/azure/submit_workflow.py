import argparse
import os
import re

from azure.identity import DefaultAzureCredential
from azure.mgmt.resource import ResourceManagementClient
from azure.mgmt.storage import StorageManagementClient
from azure.mgmt.subscription import SubscriptionClient
from azure.storage.blob import BlobClient, BlobServiceClient


def exactly_one_or_die(result, kind, filter=None, describer=None):
    found = list(result)

    if filter:
        found = [f for f in found if filter(f)]

    if len(found) == 0:
        raise ValueError(f"Could not find a {kind}!")

    if len(found) > 1:
        message = f"Found multiple {kind}!"
        if describer:
            message = message + " : " + ", ".join([describer(f) for f in found])
        raise ValueError(message)

    return found[0]


def get_subscription(credentials):
    subscription_client = SubscriptionClient(credentials)
    return exactly_one_or_die(subscription_client.subscriptions.list(),
                              "subscription",
                              describer=lambda s: s.id)


def get_resource_group(credentials, subscription):
    resource_client = ResourceManagementClient(credentials, subscription.subscription_id)

    if args.resource_group:
        resource_group_filter = lambda g: g.name == args.resource_group
        resource_group_descriptor = f"resource group '{args.resource_group}'"
    else:
        pattern = f"{os.environ['USER']}-[a-f0-9]+$"
        resource_group_descriptor = f"resource group matching pattern '{pattern}'"
        resource_group_filter = lambda g: re.match(pattern, g.name)

    return exactly_one_or_die(resource_client.resource_groups.list(),
                              resource_group_descriptor,
                              filter=resource_group_filter)


def get_storage_account(credentials, subscription, resource_group):
    storage_client = StorageManagementClient(credentials, subscription.subscription_id)
    # `az` returns storage account JSONs with a `resourceGroup` attribute, but the objects returned by the Python API do
    # not have this attribute. However the `id`s of these Python objects do contain the resource group in a predictable
    # pattern, so look for that instead.
    id_prefix = f"/subscriptions/{subscription.subscription_id}/resourceGroups/{resource_group.name}"
    return exactly_one_or_die(storage_client.storage_accounts.list(), "storage account",
                              filter=lambda a: a.id.startswith(id_prefix),
                              describer=lambda a: a.name)


def write_input_file(storage_account, path, data):
    blob_url = f"{storage_account.primary_endpoints.blob}/inputs/{path}"

    blob_client = BlobClient.from_blob_url(
        blob_url=blob_url,
        credential=credentials
    )

    blob_client.upload_blob(data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Submit workflow to Cromwell on Azure')
    parser.add_argument('--workflow', type=str, help='Workflow WDL source', required=True)
    parser.add_argument('--inputs', type=str, help='Workflow inputs', required=False)
    parser.add_argument('--resource-group', type=str, help='Azure Resource Group name', required=False)
    args = parser.parse_args()

    # https://github.com/Azure/azure-sdk-for-python/issues/22822#issuecomment-1024668507
    credentials = DefaultAzureCredential(exclude_shared_token_cache_credential=True)

    subscription = get_subscription(credentials)
    resource_group = get_resource_group(credentials, subscription)
    storage_account = get_storage_account(credentials, subscription, resource_group)

    with open(args.workflow, "rb") as workflow:
        write_input_file(storage_account, "hello/HelloAzure.wdl", workflow)
