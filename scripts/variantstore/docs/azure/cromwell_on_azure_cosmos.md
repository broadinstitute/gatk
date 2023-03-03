# **Note**

The following describes setup for Microsoft's [CromwellOnAzure](https://github.com/microsoft/CromwellOnAzure) which is
*not* the same or necessarily even similar to Cromwell on Azure in Terra. Much of what follows below may not be relevant
to running Cromwell in Terra.

# Setup for Mac

* Install the [Azure CLI locally](https://learn.microsoft.com/en-us/cli/azure/install-azure-cli-macos). 
* Authenticate locally with `az login`.
* Install .NET 6.0 from [here](https://dotnet.microsoft.com/en-us/download/dotnet). Note `CromwellOnAzure` wants .NET
  version 6.0 and will not run with .NET 7.0.
* Clone the `CromwellOnAzure` git repo from [here](https://github.com/microsoft/CromwellOnAzure).
* Build `CromwellOnAzure` with

```
cd CromwellOnAzure
dotnet build
```

* Create your deployment of `CromwellOnAzure` in the Variants Azure subscription with:

```
SUBSCRIPTION_ID=$(az account list | jq -r '.[] | select(.name | test("variants")) | .id')

dotnet src/deploy-cromwell-on-azure/bin/Debug/net6.0/deploy-cromwell-on-azure.dll --SubscriptionId ${SUBSCRIPTION_ID} --RegionName eastus --MainIdentifierPrefix ${USER}
```

The CoA deployment can take more than 20 minutes to complete. Once that's done:

```
# Fish out the CromwellOnAzure resource group name that was named after your Broad username.
# Intro to Azure resource groups:
# https://learn.microsoft.com/en-us/azure/azure-resource-manager/management/manage-resource-groups-portal#what-is-a-resource-group
RESOURCE_GROUP=$(az group list | jq -r ".[] | .name | select(test(\"^${USER}-[0-9a-f]+$\"))")
RESOURCE_GROUP_LOCATION=$(az group list | jq -r ".[] | select(.name | test(\"^${USER}-[0-9a-f]+$\")) | .location")


# Set default resource group and location
# https://learn.microsoft.com/en-us/cli/azure/azure-cli-configuration#configure-settings-using-az-config 
az configure --defaults group="${RESOURCE_GROUP}" location=${RESOURCE_GROUP_LOCATION}
```


# Azure Cosmos DB

Generally following https://learn.microsoft.com/en-us/azure/cosmos-db/nosql/quickstart-python?tabs=azure-cli%2Cpasswordless%2Clinux%2Csign-in-azure-cli%2Csync#create-an-azure-cosmos-db-account

```
COSMOS_DB_NAME=${RESOURCE_GROUP}-cosmos-db

az cosmosdb create --name ${COSMOS_DB_NAME} --locations regionName=${RESOURCE_GROUP_LOCATION}

COSMOS_DB_DOCUMENT_URI="$(az cosmosdb show --name ${COSMOS_DB_NAME} --query documentEndpoint | jq -r)"
COSMOS_DB_PRIMARY_KEY="$(az cosmosdb keys list --name ${COSMOS_DB_NAME} --type keys --query primaryMasterKey | jq -r)"

COSMOS_DB_DATA_CONTRIBUTOR_GROUP_ID=$(az ad group create --display-name "${RESOURCE_GROUP} Cosmos DB Data Contributors Group" --mail-nickname "${RESOURCE_GROUP}-cosmos-data-contributors" | jq -r .id)
# Add members of the Variants team plus the Cromwell on Azure User Assigned Managed Identity to the Cosmos DB Data
# Contributors group
VARIANTS_GROUP_ID=$(az ad group list | jq -r '.[] | select(.displayName == "DSP Variants Team") | .id')
COA_UAMI_PRINCIPAL_ID=$(az identity list | jq -r ".[] | select(.name == \"${RESOURCE_GROUP}-identity\") | .principalId")

# Get the id for the Data Contributor role
COSMOS_DB_DATA_CONTRIBUTOR_ROLE_ID=$(az cosmosdb sql role definition list --account-name $COSMOS_DB_NAME | jq -r '.[] | select(.roleName == "Cosmos DB Built-in Data Contributor") | .id')

# Assign the Cosmos DB Data Contributor role to the Cosmos DB Data Contributor group.
az cosmosdb sql role assignment create --account-name "${COSMOS_DB_NAME}" --scope "/" \
 --principal-id "${COSMOS_DB_DATA_CONTRIBUTOR_GROUP_ID}" --role-definition-id "${COSMOS_DB_DATA_CONTRIBUTOR_ROLE_ID}"

```
