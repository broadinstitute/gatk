# **Note**

This Cosmos DB exploration spike is not going to follow the pattern of the Azure SQL Database spike of trying to run in
Microsoft's Cromwell on Azure (CoA) system. We already learned in the Azure SQL Database spike that the Azure Batch VMs
spun up in CoA are not authed nor have any assigned identity. We believe from our discussions with the Workflows teams
that the Azure Batch VMs spun up in Terra on Azure will be authed as the user, which is a very significant difference
with regards to being able to do passwordless auth to Cosmos DB. For that reason the instructions that follow will
assume an environment where the user is authenticated via the Azure CLI. i.e. the following should return a JSON blob
that looks like your identity in the DSP Dev Azure Tenant:

```
az ad signed-in-user show
```

For now these instructions leverage some of the useful items that were created as part of the CoA deployment for
Azure SQL Database, specifically the resource group and User Assigned Managed Identity. If the ideas from this spike
live on it would be nice to not have dependencies on these artifacts from the previous spike.

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

Generally
following https://learn.microsoft.com/en-us/azure/cosmos-db/nosql/quickstart-python?tabs=azure-cli%2Cpasswordless%2Clinux%2Csign-in-azure-cli%2Csync#create-an-azure-cosmos-db-account

```
COSMOS_DB_NAME=${RESOURCE_GROUP}-cosmos-db

az cosmosdb create --name ${COSMOS_DB_NAME} --locations regionName=${RESOURCE_GROUP_LOCATION}

COSMOS_DB_DOCUMENT_URI="$(az cosmosdb show --name ${COSMOS_DB_NAME} --query documentEndpoint | jq -r)"
COSMOS_DB_PRIMARY_KEY="$(az cosmosdb keys list --name ${COSMOS_DB_NAME} --type keys --query primaryMasterKey | jq -r)"

# Add members of the Variants team plus the Cromwell on Azure User Assigned Managed Identity to the Cosmos DB Data
# Contributors group
VARIANTS_GROUP_ID=$(az ad group list | jq -r '.[] | select(.displayName == "DSP Variants Team") | .id')
COA_UAMI_PRINCIPAL_ID=$(az identity list | jq -r ".[] | select(.name == \"${RESOURCE_GROUP}-identity\") | .principalId")

COSMOS_DB_DATA_CONTRIBUTOR_GROUP_ID=$(az ad group create --display-name "${RESOURCE_GROUP} Cosmos DB Data Contributors Group" --mail-nickname "${RESOURCE_GROUP}-cosmos-data-contributors" | jq -r .id)
for id in ${VARIANTS_GROUP_ID} ${COA_UAMI_PRINCIPAL_ID}
do
  az ad group member add --group ${COSMOS_DB_DATA_CONTRIBUTOR_GROUP_ID} --member-id ${id}
done

# Get the id for the Data Contributor role
COSMOS_DB_DATA_CONTRIBUTOR_ROLE_ID=$(az cosmosdb sql role definition list --account-name $COSMOS_DB_NAME | jq -r '.[] | select(.roleName == "Cosmos DB Built-in Data Contributor") | .id')

# Assign the Cosmos DB Data Contributor role to the Cosmos DB Data Contributor group.
# Unfortunately this does not seem to be working. If we give the Variants Team group this role everything seems to work.
# Are nested groups really not supported for role assignments?
# for id in ${COSMOS_DB_DATA_CONTRIBUTOR_GROUP_ID}

# Assign the Cosmos DB Data Contributor role to the Variants team and the CoA UAMI.
for id in ${VARIANTS_GROUP_ID} ${COA_UAMI_PRINCIPAL_ID}
do
  az cosmosdb sql role assignment create --account-name "${COSMOS_DB_NAME}" --scope "/" \
    --principal-id "${id}" --role-definition-id "${COSMOS_DB_DATA_CONTRIBUTOR_ROLE_ID}"
done

az cosmosdb sql database create --account-name ${COSMOS_DB_NAME} --name cosmos_gvs

az cosmosdb sql container create --account-name ${COSMOS_DB_NAME} --database-name cosmos_gvs --partition-key-path "/sampleId" --name sample_info

```

Download the samples table from the quickstart, process into a FOFN with VCFs and indexes:

```
gtail --lines=+2 sample.tsv | awk '{print $2,"\n", $3}' | sed 's/ *//g' > files.txt
```

`gtail` is GNU tail as Homebrew likes to call it. Next copy all of these VCFs from their scattered locations to a single
location to facilitate the transfer to Azure:

```
cat files.txt | gcloud storage cp -I gs://<some bucket>/<assembled quickstart vcf path>
```

`gs://<some bucket>/<assembled quickstart vcf path>` will be the source for the `azcopy copy` command. 

Download Microsoft's `azcopy` utility from [here](https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10)
and put it in your `PATH`.

Create a new storage container in the Azure portal to hold these VCFs and indexes. 

From 'Shared access tokens' for this container on the navigation pane of the portal, select the appropriate permissions,
IP address, and time range and press the big blue 'Generate SAS token and URL' button. Copy the value in the 'Blob SAS
URL' field; this will be the destination for the `azcopy copy` command.

Somewhat annoyingly, `azcopy copy` insists on using a service account key to do this copying and will not use your
personal auth. Go to the Google Cloud Console for the project containing your source bucket, then navigate to IAM ->
Service Accounts and select the 'Compute Engine default service account'. Select the 'KEYS' blade and add a key. This
will download a service account JSON to your local machine. Set the Google credentials environment variable like so:

```
export GOOGLE_APPLICATION_CREDENTIALS=<path to service account json>
```

And finally you can run `azcopy`:

```
azcopy cp 'https://storage.cloud.google.com/<some bucket>/<assembled quickstart vcf path>/*' '<Blob SAS URL from Azure Portal>'
```

Make sure to enclose both the source and destination paths in single quotes to escape metacharacters.
