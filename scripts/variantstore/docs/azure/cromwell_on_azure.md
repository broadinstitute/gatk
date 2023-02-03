# Setup

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
# Get the Variants subscription ID from Azure Portal.
SUBSCRIPTION_ID=...

dotnet src/deploy-cromwell-on-azure/bin/Debug/net6.0/deploy-cromwell-on-azure.dll --SubscriptionId ${SUBSCRIPTION_ID} --RegionName eastus --MainIdentifierPrefix ${USER}
```

This command can take about 20 minutes to complete.

Next, launch a [Cloud Shell](https://portal.azure.com/#cloudshell/) from the Azure Portal. In the Cloud Shell:

```
# Specify the same Azure subscription id as we used in the local terminal above.
SUBSCRIPTION_ID=...

# Set our subscription id as the default so we don't need to include it as an argument to every command.
az account set --subscription ${SUBSCRIPTION_ID}

# Fish out our Broad username from this JSON blob stored in an environment variable (in Azure Cloud Shell $USER is set
# to our first name).
BROAD_USER=$(echo -n $ACC_STORAGE_PROFILE | jq -r '.fileShareName | split("-") | .[1]')

# Fish out the resource group name that was named after your Broad username.
RESOURCE_GROUP=$(az group list | jq -r ".[] | .name | select(test(\"^${BROAD_USER}-[0-9a-f]+$\"))")

# Make a strong random password (save this someplace).
ADMIN_PASSWORD=$(dd if=/dev/urandom bs=30 count=1 2> /dev/null | base64)

SQL_SERVER=${RESOURCE_GROUP}-azure-sql-server
SQL_DATABASE=${RESOURCE_GROUP}-azure-sql-database

# Mostly taken from here https://learn.microsoft.com/en-us/azure/azure-sql/database/scripts/create-and-configure-database-cli?view=azuresql

# Set resource group and sql server defaults to avoid having to include them repeatedly in `az` commands.
az configure --defaults group="${RESOURCE_GROUP}" sql-server="${SQL_SERVER}"

# Create the server for "serverless" Azure SQL Database
az sql server create --name ${SQL_SERVER} --location "eastus" --admin-user ${BROAD_USER} --admin-password "${ADMIN_PASSWORD}"
az sql db create --name ${SQL_DATABASE} --edition GeneralPurpose --family Gen5 --capacity 2 --zone-redundant true

# Get the external IP of this machine and create a firewall rule that allows it to connect to the Azure SQL Database. 
# https://dev.to/0xbf/get-your-public-ip-address-from-command-line-57l4
EXTERNAL_IP=$(curl --silent ifconfig.me)

az sql server firewall-rule create -n AllowYourIp --start-ip-address $EXTERNAL_IP --end-ip-address $EXTERNAL_IP

sqlcmd -S tcp:${SQL_SERVER}.database.windows.net,1433 -d ${SQL_DATABASE} -U ${BROAD_USER} -P ${ADMIN_PASSWORD} -N -l 30

1> select "Hello Azure SQL Database!"
2> go

-------------------------
Hello Azure SQL Database!


# TODO: Figure out AAD, username/password is ick

```
