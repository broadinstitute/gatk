# **Note**

This Cosmos DB exploration spike builds on the learnings from the Azure SQL Database spike, most notably in not trying
to run inside Microsoft's Cromwell on Azure (CoA) system. The Azure SQL Database spike revealed that the Azure Batch VMs
spun up in CoA are not authed nor have any assigned identity. However from our discussions with the Workflows teams we
believe that the Azure Batch VMs spun up in Terra on Azure will be authed as the user. For that reason the instructions
that follow assume an environment where the user is authenticated via the Azure CLI. i.e. the following should return a
JSON blob that looks like your identity in the DSP Dev Azure Tenant:

```
az ad signed-in-user show
```

For now these instructions leverage some of the useful items that were created as part of the CoA deployment for
Azure SQL Database, specifically the resource group and User Assigned Managed Identity. If the ideas from this spike
live on it would be nice to not have dependencies on these artifacts from the previous Azure SQL Database spike.

```
# Fish out the CromwellOnAzure resource group name that was named after your Broad username.
# Intro to Azure resource groups:
# https://learn.microsoft.com/en-us/azure/azure-resource-manager/management/manage-resource-groups-portal#what-is-a-resource-group
export RESOURCE_GROUP=$(az group list | jq -r ".[] | .name | select(test(\"^${USER}-[0-9a-f]+$\"))")
export RESOURCE_GROUP_LOCATION=$(az group list | jq -r ".[] | select(.name | test(\"^${USER}-[0-9a-f]+$\")) | .location")


# Set default resource group and location
# https://learn.microsoft.com/en-us/cli/azure/azure-cli-configuration#configure-settings-using-az-config 
az configure --defaults group="${RESOURCE_GROUP}" location=${RESOURCE_GROUP_LOCATION}
```

# Azure Cosmos DB

## Basics

Generally
following [this](https://learn.microsoft.com/en-us/azure/cosmos-db/nosql/quickstart-python?tabs=azure-cli%2Cpasswordless%2Clinux%2Csign-in-azure-cli%2Csync#create-an-azure-cosmos-db-account):

```
export COSMOS_DB_NAME=${RESOURCE_GROUP}-cosmos-db

az cosmosdb create --name ${COSMOS_DB_NAME} --locations regionName=${RESOURCE_GROUP_LOCATION}

export COSMOS_ENDPOINT="$(az cosmosdb show --name ${COSMOS_DB_NAME} --query documentEndpoint | jq -r)"
export COSMOS_KEY="$(az cosmosdb keys list --name ${COSMOS_DB_NAME} --type keys --query primaryMasterKey | jq -r)"

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
export COSMOS_DB_DATA_CONTRIBUTOR_ROLE_ID=$(az cosmosdb sql role definition list --account-name $COSMOS_DB_NAME | jq -r '.[] | select(.roleName == "Cosmos DB Built-in Data Contributor") | .id')

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

# How many documents are in the container?
export CONTAINER_NAME=vets
az cosmosdb sql container show --database-name cosmos_gvs --name ${CONTAINER_NAME} --account-name $COSMOS_DB_NAME  | jq -r '..|.documentCount? //empty'

```

## VCFs

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

Download Microsoft's `azcopy` utility
from [here](https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10)
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

Note: the code below creates Cosmos containers with a fixed, default bandwidth of 400 RU/s which is way too slow for
variant and reference data loading. The containers used for this spike were created through the Azure portal with and
autoscale bandwidth from 10K to 100K RU/s. It should be possible to modify the code below to create containers with
similar parameters.

```
for container in vets ref_ranges
do
  az cosmosdb sql container create --account-name ${COSMOS_DB_NAME} --database-name cosmos_gvs --partition-key-path "/sample_id" --name $container
done
```

## Avros

Extract variant data from a GCP Quickstart Integration run. This is what was actually used to drive the
`cosmos_load_data.sc` script:

```
EXPORT DATA OPTIONS(
                uri='gs://someplace/scratch/avro/quickstart-raw/vets/vet_001/vet_001_*.avro', format='AVRO', compression='SNAPPY') AS
                SELECT v.sample_id as sample_id, location, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, QUALapprox, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, call_GT, call_AD, call_GQ, call_PGT, call_PID, call_PL
                FROM `gvs-internal.quickit_vs_639_hail_testing_spike_ed64917_hail.vet_001` v
                INNER JOIN `gvs-internal.quickit_vs_639_hail_testing_spike_ed64917_hail.sample_info` s ON s.sample_id = v.sample_id
                WHERE withdrawn IS NULL AND
                is_control = false
                ORDER BY location
```

Extract reference data:

```
EXPORT DATA OPTIONS(
                uri='gs://someplace/scratch/avro/quickstart-raw/ref_ranges/ref_ranges_001/ref_ranges_001_*.avro', format='AVRO', compression='SNAPPY') AS
                SELECT r.sample_id as sample_id, location, length, state
                FROM `gvs-internal.quickit_vs_639_hail_testing_spike_ed64917_hail.ref_ranges_001` r
                INNER JOIN `gvs-internal.quickit_vs_639_hail_testing_spike_ed64917_hail.sample_info` s ON s.sample_id = r.sample_id
                WHERE withdrawn IS NULL AND
                is_control = false
                ORDER BY location
```

# Learnings

## Loading is slow, error prone, and expensive

I made `vets` and `ref_ranges` Cosmos containers for this spike and set their bandwidths to autoscale from 10K to 100K
RU/s (up from the default of 400 to 4000 RU/s), the maximum before the UI shows a scary checkbox that I would have to
click to acknowledge how much I would be charged for the requested bandwidth. However even 100K RU/s was painfully slow
for loading references data. With all indexing turned off, loading the references data for the 10 quickstart samples
would take about 18 hours and cost tens of dollars per sample. That cost is absolutely nuts, orders of magnitude more
than we pay per sample for the full variant calling pipeline in GCP. While it's possible we might get around the load
time bottleneck with a different container layout (e.g. one per sample), the cost problem is independent of that and
absolutely needs to be solved if we're seriously going to consider Cosmos as a storage solution for GVS on Azure. And
all of this is assuming the data can even be loaded correctly, something I never managed to do with the reference data
due to spurious 403 errors.

To see RU consumption (narrow this to a single container if appropriate):

```
Azure Portal -> CosmosDB -> Metrics
Metric Namespace = 'Cosmos DB standard metrics'
Metric = 'Normalized RU Consumption'
```

To see the number of documents written to a Cosmos container:

```
export CONTAINER_NAME=vets
az cosmosdb sql container show --database-name cosmos_gvs --name ${CONTAINER_NAME} --account-name ${COSMOS_DB_NAME} |
  jq -r '..|.documentCount? //empty' |
  paste --serial --delimiters=+ |
  bc
```

Picking this apart:

* The `az` command returns an informative JSON blob for our Cosmos container of interest.
* The `jq` command recursively picks out values for all fields keyed by `documentCount`. There can be multiple of these
  fields, one for each physical (?) Cosmos partition.
* The `paste` command joins all the `documentCount` lines together on one line, delimited by '+'.
* The `bc` line evaluates this as an arithmetic expression (a sum).

## What goes up does not come (all the way) down

Cosmos allows for container bandwidth to be edited after creation. After I loaded the vets data and (most of) the
references data, I wanted to scale down my Cosmos container bandwidth to save money. However it seems that Cosmos does
not allow *arbitrary* downsizing of containers, but only up to a maximum of 10x smaller. So setting container bandwidth
very high to get decent load performance has the very undesirable consequence of setting a high floor for container
costs after loading is complete.

## Reactor not Reacting

The V4 Cosmos Java client library makes extensive use of
the [Reactor library](https://projectreactor.io/docs/core/release/reference/) in
its [Bulk Executor APIs](https://learn.microsoft.com/en-us/azure/cosmos-db/bulk-executor-overview).
This library aims to "[provide efficient demand management](https://projectreactor.io/)". Unfortunately in the
context of this spike, the loader code I wrote did not appear to behave in a reactive manner. The Avro "source"
generated records far faster than the Cosmos "sink" could process them, causing unflushed records to pile up inside the
Cosmos client code. The original version of this code created a more elegant `flatMap`ped stream of Avro records which
were transformed into Jackson `ObjectNode`s and `Flux`ed into the Cosmos client library. However this invariably led to
memory pressure, thrashing from constant garbage collection, and eventually a crash with a 408 error. Included among the
very detailed error messages at the time of the crash was the fact that the Cosmos client library was holding
millions of unwritten records:

```
03:40:03.554 [bulk-executor-bounded-elastic-265] INFO com.azure.cosmos.implementation.batch.BulkExecutor - BulkExecutor.execute flux terminated - Signal: cancel - # left items 6514672, Context: BulkExecutor-32[n/a], Thread[Name: bulk-executor-bounded-elastic-265,Group: main, isDaemon: true, Id: 328]
```

The code committed as part of this spike replaced the `flatMap` with an iterator over the Avro files. While
aesthetically unfortunate, this structure did enable loading variant data cleanly, although the reference data ran into
a slew of unexplained 403 failures after several hours of running.

Per Microsoft documentation the "solution" for High CPU utilization-driven 408 errors it to
"[scale the client up and out](https://learn.microsoft.com/en-us/azure/cosmos-db/nosql/troubleshoot-java-sdk-request-timeout#solution)".
At least in this case that advice doesn't seem consistent with the intent of the Reactor library, though conceivably
there could be cases where the CPU load was due to a source that computationally struggled to make items quickly enough
to fill a fast-draining sink.

This is my first time programming to the Reactor API and I while it's certainly possible I have done something
horrifically wrong in setting up the spike program, I did closely follow the example of the trivially
small [sample code](https://github.com/Azure-Samples/azure-cosmos-java-sql-api-samples/blob/main/src/main/java/com/azure/cosmos/examples/bulk/async/SampleBulkQuickStartAsync.java).
From the limited reading I have done on Reactor, it does seem that it is
the [responsibility of the consumer](https://stackoverflow.com/a/57298393/21269164) to request
only as much data as it can actually handle from the producer. See
also [this piece](https://projectreactor.io/docs/core/release/reference/#reactor.hotCold) on hot versus cold sources.
