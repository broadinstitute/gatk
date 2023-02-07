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
RESOURCE_GROUP=$(az group list | jq -r ".[] | .name | select(test(\"^${USER}-[0-9a-f]+$\"))")

# Make a strong random password. You should probably save this someplace even though it's currenty not used in these
# instructions beyond creating the Azure SQL Server.
ADMIN_PASSWORD=$(dd if=/dev/urandom bs=30 count=1 2> /dev/null | base64)

SQL_SERVER=${RESOURCE_GROUP}-azure-sql-server
SQL_DATABASE=${RESOURCE_GROUP}-azure-sql-database

# Mostly taken from here https://learn.microsoft.com/en-us/azure/azure-sql/database/scripts/create-and-configure-database-cli?view=azuresql

# Set resource group and sql server defaults to avoid having to include them repeatedly in `az` commands.
az configure --defaults group="${RESOURCE_GROUP}" sql-server="${SQL_SERVER}"

# Create the server for "serverless" Azure SQL Database. These two steps take a minute or so each.
az sql server create --name ${SQL_SERVER} --location "eastus" --admin-user ${USER} --admin-password "${ADMIN_PASSWORD}"
az sql db create --name ${SQL_DATABASE} --edition GeneralPurpose --family Gen5 --capacity 2 --zone-redundant true

# Get the external IP of this machine and create a firewall rule that allows it to connect to the Azure SQL Database. 
# https://dev.to/0xbf/get-your-public-ip-address-from-command-line-57l4
EXTERNAL_IP=$(curl --silent ifconfig.me)

az sql server firewall-rule create -n AllowYourIp --start-ip-address $EXTERNAL_IP --end-ip-address $EXTERNAL_IP

# Create an AD admin group for the Azure SQL Database.
AZ_SQLDB_AD_ADMIN_GROUP_ID=$(az ad group create --display-name "${RESOURCE_GROUP} Azure SQL Database AD Admin Group" --mail-nickname "${RESOURCE_GROUP}-ad-admin" | jq -r .id)

# Add the whole Variants team to this group.
VARIANTS_GROUP_ID=$(az ad group list | jq -r '.[] | select(.displayName == "DSP Variants Team") | .id')
az ad group member add --group $AZ_SQLDB_AD_ADMIN_GROUP_ID --member-id ${VARIANTS_GROUP_ID}

# Also add the User Assigned Managed Identity created by the CromwellOnAzure deployer to this group.
COA_UAMI_ID=$(az identity list | jq -r ".[] | select(.name == \"${RESOURCE_GROUP}-identity\") | .principalId")
az ad group member add --group $AZ_SQLDB_AD_ADMIN_GROUP_ID --member-id ${COA_UAMI_ID}

# Make the AD Admin group the AD Admin for the Azure SQL Server. All members of this group will be able to act as
# Azure AD Admins for this server.
az sql server ad-admin create --object-id ${AZ_SQLDB_AD_ADMIN_GROUP_ID} --display-name "${RESOURCE_GROUP} Azure SQL Database AD Admin"

# Get a database access token.
# https://learn.microsoft.com/en-us/sql/connect/odbc/linux-mac/connecting-with-sqlcmd?view=azuresqldb-current
SQLCMDPASSWORD=$(az account get-access-token --resource https://database.windows.net --output tsv | cut -f 1 | tr -d '\n' | iconv -f ascii -t UTF-16LE)

# Say hello to Azure SQL Database!
sqlcmd -S tcp:${SQL_SERVER}.database.windows.net,1433 -d ${SQL_DATABASE} -G -Q 'select @@version as "Hello Azure SQL Database!"'
                         
Hello Azure SQL Database!                                                                                                                                                                                                                                                                                   
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Microsoft SQL Azure (RTM) - 12.0.2000.8 
	Jan 12 2023 05:25:39 
	Copyright (C) 2022 Microsoft Corporation


(1 rows affected)

```
