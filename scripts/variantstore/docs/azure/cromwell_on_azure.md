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

# Add an allow rule for all Azure traffic so the Azure Batch VMs spun up by Cromwell on Azure can connect.
# Hopefully this can be done more narrowly in the future, but likely we won't know the Batch VMs exact IP addresses in advance. 
# https://learn.microsoft.com/en-us/azure/azure-sql/database/firewall-configure?view=azuresql#connections-from-inside-azure
az sql server firewall-rule create -n AllowAllWindowsAzureIps --start-ip-address 0.0.0.0 --end-ip-address 0.0.0.0 

# Create an AD admin group for the Azure SQL Database.
AZ_SQLDB_AD_ADMIN_GROUP_ID=$(az ad group create --display-name "${RESOURCE_GROUP} Azure SQL Database AD Admin Group" --mail-nickname "${RESOURCE_GROUP}-ad-admin" | jq -r .id)

# Add the whole Variants team to this group.
VARIANTS_GROUP_ID=$(az ad group list | jq -r '.[] | select(.displayName == "DSP Variants Team") | .id')
az ad group member add --group $AZ_SQLDB_AD_ADMIN_GROUP_ID --member-id ${VARIANTS_GROUP_ID}

# Also add the User Assigned Managed Identity created by the CromwellOnAzure deployer to this group.
COA_UAMI_PRINCIPAL_ID=$(az identity list | jq -r ".[] | select(.name == \"${RESOURCE_GROUP}-identity\") | .principalId")
az ad group member add --group $AZ_SQLDB_AD_ADMIN_GROUP_ID --member-id ${COA_UAMI_PRINCIPAL_ID}

# Grab the server ID as we want to the SQL Security Manager role to the UAMI to be able to add the VM to the server
# firewall's allowed IP addresses.
SQL_SERVER_ID=$(az sql server list | jq -r ".[] | select(test(\"${RESOURCE_GROUP}\")) | .id"

az role assignment create --role "SQL Security Manager" --assignee "${COA_UAMI_PRINCIPAL_ID}" --scope "${SQL_SERVER_ID}"

# Make the AD Admin group the AD Admin for the Azure SQL Server. All members of this group will be able to act as
# Azure AD Admins for this server.
az sql server ad-admin create --object-id ${AZ_SQLDB_AD_ADMIN_GROUP_ID} --display-name "${RESOURCE_GROUP} Azure SQL Database AD Admin"

get_db_token() {
    # Fetches a database access token with ~1 hour TTL, ASCII / UTF-8 encoded and newline terminated
    # https://learn.microsoft.com/en-us/sql/connect/odbc/linux-mac/connecting-with-sqlcmd?view=azuresqldb-current
    az account get-access-token --resource https://database.windows.net --query accessToken --output tsv
}

# Be aware that there are at least two versions of `sqlcmd` circulating and they are not compatible with respect to the
# handling of Azure Active Directory credentials. The instructions below are for the `mssql-tools` version (as the
# package is called in both Ubuntu and Homebrew). Note that the `sqlcmd` formula in Homebrew installs a Golang-based
# version of `sqlcmd` which uses an environment variable SQLCMDPASSWORD to hold the access token rather than a file
# specified with the -P option.

# Say hello to Azure SQL Database!
sqlcmd -S tcp:${SQL_SERVER}.database.windows.net,1433 -d ${SQL_DATABASE} -G -Q 'select @@version as "Hello Azure SQL Database!"' -P =(get_db_token | tr -d '\n' | iconv -f ascii -t UTF-16LE)

                         
Hello Azure SQL Database!                                                                                                                                                                                                                                                                                   
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Microsoft SQL Azure (RTM) - 12.0.2000.8 
	Jan 12 2023 05:25:39 
	Copyright (C) 2022 Microsoft Corporation


(1 rows affected)

# Yes sqlcmd (or more specifically the Microsoft ODBC Driver) really does require newline-stripped UTF-16 Little
# Endian encoded tokens and will fail to log in with no useful diagnostics if it gets anything else.
# =() is some zsh trickery that uses temporary files; the usual bash <() construct does not work here.
