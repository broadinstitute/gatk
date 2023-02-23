import argparse
import sys

import inflection
import os
import re

from azure.identity import DefaultAzureCredential
from azure.mgmt.resource import ResourceManagementClient
from azure.mgmt.sql import SqlManagementClient
from azure.mgmt.storage import StorageManagementClient
from azure.mgmt.subscription import SubscriptionClient
from azure.storage.blob import BlobClient, BlobServiceClient

from pathlib import Path
from uuid import uuid4


def exactly_one_or_die(result, kind, filter=None, describer=None):
    """
    Pull out exactly one item from the result or die trying.
    """
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


def get_resource_group_filter(resource_group_name):
    if resource_group_name:
        resource_group_filter = lambda g: g.name.startswith(resource_group_name)
        resource_group_descriptor = f"resource group '{resource_group_name}'"
    else:
        pattern = f"{os.environ['USER']}-[a-f0-9]+$"
        resource_group_descriptor = f"resource group matching pattern '{pattern}'"
        resource_group_filter = lambda g: re.match(pattern, g.name)
    return resource_group_filter, resource_group_descriptor


def get_resource_group(credentials, subscription, resource_group_name=None):
    resource_group_filter, resource_group_descriptor = get_resource_group_filter(resource_group_name)
    resource_client = ResourceManagementClient(credentials, subscription.subscription_id)
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


def get_sql_server(credentials, subscription):
    sql_management_client = SqlManagementClient(credentials, subscription.subscription_id)
    return exactly_one_or_die(sql_management_client.servers.list(), 'Azure SQL Server')


def get_sql_database(credentials, subscription, resource_group, server):
    sql_management_client = SqlManagementClient(credentials, subscription.subscription_id)
    resource_group_filter, _ = get_resource_group_filter(resource_group.name)
    return exactly_one_or_die(sql_management_client.databases.list_by_server(resource_group.name, server.name),
                              'Azure SQL Database',
                              resource_group_filter)


def get_blob_service_client():
    return BlobServiceClient.from_connection_string(os.getenv('AZURE_CONNECTION_STRING'))


def generate_inputs_json(workflow_name, access_token_storage_path, python_script_storage_path,
                         ammonite_script_storage_path, server_name, database_name):
    return f"""
{{
  "{workflow_name}.utf8_token_file": "{access_token_storage_path}",
  "{workflow_name}.python_script": "{python_script_storage_path}",
  "{workflow_name}.ammonite_script": "{ammonite_script_storage_path}",
  "{workflow_name}.sql_server": "{server_name}",
  "{workflow_name}.sql_database": "{database_name}"
}}
"""


def generate_trigger_json(workflow_storage_path, inputs_storage_path):
    """
    Creates a trigger JSON of the form accepted by CromwellOnAzure. This is conceptually similar to the JSON payload
    of a workflow submission POST that would normally go to Cromwell's REST interface.
    """
    return f"""
{{
  "WorkflowUrl": "{workflow_storage_path}",
  "WorkflowInputsUrl": "{inputs_storage_path}",
  "WorkflowInputsUrls": null,
  "WorkflowOptionsUrl": null,
  "WorkflowDependenciesUrl": null
}}
    """.strip()


if __name__ == '__main__':
    description = """

    Cromwell on Azure (CoA) "Hello Azure!" workflow submission script does the following:

    1. Stages the workflow and its `File` inputs (scripts, database access token) to the CoA inputs container.
    2. Generates an inputs JSON corresponding to the inputs file from 1. and stages this to the CoA inputs container.
    3. Generates a trigger JSON for this workflow + inputs and stages this to the CoA workflows container under 'new'.

    The script does *not* attempt to poll the submitted workflow for status, this is simple "fire and forget".
    """

    parser = argparse.ArgumentParser(allow_abbrev=False, description=description)
    parser.add_argument('--workflow', type=str, help='Workflow WDL source', required=True)
    parser.add_argument('--python-script', type=str, help="Hello World Python script", required=True)
    parser.add_argument('--ammonite-script', type=str, help="Hello World Ammonite script", required=True)
    parser.add_argument('--sql-server', type=str, help='Azure SQL Server name', required=True)
    parser.add_argument('--sql-database', type=str, help='Azure SQL Server database', required=True)
    parser.add_argument('--utf8-access-token', type=str, help='UTF-8 encoded Azure SQL Database access token',
                        required=True)
    parser.add_argument('--resource-group', type=str, help='Azure Resource Group name', required=False)
    args = parser.parse_args()

    if not os.getenv('AZURE_CONNECTION_STRING'):
        raise ValueError("Must define 'AZURE_CONNECTION_STRING' as a SAS token with write permissions to the CoA storage account, see https://learn.microsoft.com/en-us/azure/storage/common/storage-configure-connection-string#store-a-connection-string")

    # The shared token cache was causing issues with attempts to use expired tokens so disable that.
    # `DefaultAzureCredential` appears to fall back to Azure CLI credentials which works fine.
    # https://github.com/Azure/azure-sdk-for-python/issues/22822#issuecomment-1024668507
    credentials = DefaultAzureCredential(exclude_shared_token_cache_credential=True)

    subscription = get_subscription(credentials)
    resource_group = get_resource_group(credentials, subscription, args.resource_group)
    storage_account = get_storage_account(credentials, subscription, resource_group)
    sql_server = get_sql_server(credentials, subscription)
    sql_database = get_sql_database(credentials, subscription, resource_group, sql_server)

    blob_service_client = get_blob_service_client()
    inputs_client = blob_service_client.get_container_client('inputs')
    # `name` is the filename without leading directory components.
    # e.g. name for /path/to/Hello.wdl is Hello.wdl
    # `stem` is the filename without leading directory components and without an extension.
    # e.g. stem for /path/to/Hello.wdl is Hello
    workflow_path = Path(args.workflow)

    # Stage the workflow into /<storage container>/inputs/<snake cased workflow name>/<workflow file>.
    with open(args.workflow, "rb") as workflow:
        # `inflection.underscore` snake-cases Pascal-cased workflow names. Not strictly required here but nice.
        # e.g. "HelloAzure" ==> "hello_azure"
        blob_address = f"{inflection.underscore(workflow_path.stem)}/{workflow_path.name}"
        blob_client = inputs_client.get_blob_client(blob_address)
        blob_client.upload_blob(workflow, overwrite=True)
        workflow_storage_path = f'/{storage_account.name}/inputs/{blob_address}'

    # Stage the Python script into /<storage container>/inputs/<snake cased workflow name>/<Python script>.
    python_script_path = Path(args.python_script)
    with open(args.python_script, "rb") as python_script:
        blob_address = f"{inflection.underscore(workflow_path.stem)}/{python_script_path.name}"
        blob_client = inputs_client.get_blob_client(blob_address)
        blob_client.upload_blob(python_script, overwrite=True)
        python_script_storage_path = f'/{storage_account.name}/inputs/{blob_address}'

    # Stage the Ammonite script into /<storage container>/inputs/<snake cased workflow name>/<Ammonite script>.
    ammonite_script_path = Path(args.ammonite_script)
    with open(args.ammonite_script, "rb") as ammonite_script:
        blob_address = f"{inflection.underscore(workflow_path.stem)}/{ammonite_script_path.name}"
        blob_client = inputs_client.get_blob_client(blob_address)
        blob_client.upload_blob(ammonite_script, overwrite=True)
        ammonite_script_storage_path = f'/{storage_account.name}/inputs/{blob_address}'

    # Stage the access token into /<storage container>/inputs/<user name>/db_access_token.txt
    with open(args.utf8_access_token, "rb") as token:
        blob_address = f"{os.environ['USER']}/db_access_token.txt"
        blob_client = inputs_client.get_blob_client(blob_address)
        blob_client.upload_blob(token, overwrite=True)
        access_token_storage_path = f'/{storage_account.name}/inputs/{blob_address}'

    inputs_json = generate_inputs_json(
        workflow_path.stem, access_token_storage_path, python_script_storage_path, ammonite_script_storage_path,
        args.sql_server, args.sql_database)

    # Stage the workflow inputs into /<storage container>/inputs/<snake cased workflow name>/<workflow name>.inputs.json
    blob_address = f"{inflection.underscore(workflow_path.stem)}/{workflow_path.stem}.inputs.json"
    blob_client = inputs_client.get_blob_client(blob_address)
    blob_client.upload_blob(bytes(inputs_json, 'utf8'), overwrite=True)
    inputs_storage_path = f'/{storage_account.name}/inputs/{blob_address}'

    # Create the trigger JSON and stage into /<storage account name>/workflows/new.
    trigger_json = generate_trigger_json(workflow_storage_path, inputs_storage_path)
    workflows_client = blob_service_client.get_container_client('workflows')

    blob_address = f'new/{workflow_path.stem}-{uuid4()}.json'
    blob_client = workflows_client.get_blob_client(blob_address)
    blob_client.upload_blob(bytes(trigger_json, 'utf8'))

    print(f"Trigger JSON staged to /{storage_account.name}/workflows/{blob_address}.")
