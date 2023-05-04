# Importing vet and ref_ranges data to Azure SQL Database, serverless style!

Similar to and building on the [earlier VS-879 spike](azure_sql_database_ingest.md), this spike imports vet and ref
ranges data from Azure Blob Storage to an Azure SQL Database. However this time the server will be... serverless and
costs will be measured.

All steps below run from an Azure VM in which I was authed as myself.

## Create a new resource group

```
export RESOURCE_GROUP="azsql-serverless-group"
az group create --name "${RESOURCE_GROUP}" --location eastus
```

Make this resource group the default:

```
az configure --defaults group="${RESOURCE_GROUP}"
```


## Create the server

Yes we need an Azure SQL Database server even when we use a serverless configuration...

Create the server:

```
$ export SQL_SERVER="azsql-server-serverless"

$ az sql server create \
    --enable-ad-only-auth \
    --enable-public-network \
    --location eastus \
    --external-admin-principal-type User \
    --external-admin-sid "$(az ad signed-in-user show | jq -r .id)" \
    --external-admin-name "$(az ad signed-in-user show | jq -r .mail)" \
    --name "${SQL_SERVER}"
```

Make this server the default:

```
$ az configure --defaults sql-server="${SQL_SERVER}"
```

Set up a firewall rule so we can connect to the server from this Azure VM:

```
# Get the external IP of this machine and create a firewall rule that allows it to connect to the Azure SQL Database. 
# https://dev.to/0xbf/get-your-public-ip-address-from-command-line-57l4
$ EXTERNAL_IP=$(curl --silent ifconfig.me)

$ az sql server firewall-rule create -n AllowYourIp --start-ip-address $EXTERNAL_IP --end-ip-address $EXTERNAL_IP
```

## Create the database

Note the database is specified as 'Serverless', not the server.

```
$ export SQL_DATABASE="azsql-database-serverless"

$ az sql db create \
    --name "${SQL_DATABASE}" \
    -e GeneralPurpose \
    -f Gen5 \
    -c 2 \
    --auto-pause-delay 10 \
    --compute-model Serverless
```

## Generate a database token

```
$ get_db_token() {
    # Fetches a database access token with ~1 hour TTL, ASCII / UTF-8 encoded and newline terminated
    # https://learn.microsoft.com/en-us/sql/connect/odbc/linux-mac/connecting-with-sqlcmd?view=azuresqldb-current
    az account get-access-token --resource https://database.windows.net --query accessToken --output tsv
}

$ get_db_token | tr -d '\n' | iconv -f ascii -t UTF-16LE > token.txt
```

## Ingest

Largely copying from the previous [VS-879 ingest spike](azure_sql_database_ingest.md) and presuming everything that
spike put in blob storage is still there:

### Set required environment variables

```
$ export MASTER_KEY_PASSWORD='...'
$ export CSV_CONTAINER_SAS_TOKEN='...'
$ export STORAGE_ACCOUNT_NAME='...'
$ export STORAGE_CONTAINER_NAME='...'
```

### Generate loader script prologue

```
$ cat > loader_script.sql <<FIN


CREATE MASTER KEY ENCRYPTION BY PASSWORD = '${MASTER_KEY_PASSWORD}';


CREATE DATABASE SCOPED CREDENTIAL MyAzureBlobStorageCredential
WITH IDENTITY = 'SHARED ACCESS SIGNATURE',
SECRET = '${CSV_CONTAINER_SAS_TOKEN}';


CREATE EXTERNAL DATA SOURCE MyAzureBlobStorage
WITH ( 
    TYPE = BLOB_STORAGE,
    LOCATION = 'https://${STORAGE_ACCOUNT_NAME}.blob.core.windows.net/${STORAGE_CONTAINER_NAME}',
    CREDENTIAL= MyAzureBlobStorageCredential
);


CREATE TABLE ref_ranges (
    location BIGINT,
    sample_id INT,
    length SMALLINT,
    state TINYINT
);


CREATE TABLE vets (
    location BIGINT,
    sample_id INT,
    ref VARCHAR(1024),
    alt VARCHAR(4096),
    AS_RAW_MQ VARCHAR(255),
    AS_RAW_MQRankSum VARCHAR(255),
    QUALapprox VARCHAR(255),
    AS_QUALapprox VARCHAR(255),
    AS_RAW_ReadPosRankSum VARCHAR(255),
    AS_SB_TABLE VARCHAR(255),
    AS_VarDP VARCHAR(255),
    call_GT VARCHAR(255),
    call_AD VARCHAR(255),
    call_GQ VARCHAR(255),
    call_PGT VARCHAR(255),
    call_PID VARCHAR(1024),
    call_PL VARCHAR(255)
);


FIN

```

```

for kind in vet ref_ranges
do
    list_file="${kind}_list.txt"
    az storage blob list \
        --account-name "${STORAGE_ACCOUNT_NAME}" \
        --container-name "${STORAGE_CONTAINER_NAME}" \
        --prefix "${kind}" \
        --sas-token "$SAS_TOKEN" | \ 
        jq -r '.[] | .name' | \
        grep -E '.csv$' > \
        "${list_file}"
        
    for file in $(cat "${list_file}")
    do
        if [[ "${kind}" -eq "vet" ]]
        then
            cat >> loader_script.sql <<FIN
INSERT INTO vets with (TABLOCK)
  SELECT 
      CONVERT(BIGINT, location) AS location,
      CONVERT(INT, sample_id) AS sample_id,
      ref,
      alt,
      AS_RAW_MQ,
      AS_RAW_MQRankSum,
      QUALapprox,
      AS_QUALapprox,
      AS_RAW_ReadPosRankSum,
      AS_SB_TABLE,
      AS_VarDP,
      call_GT,
      call_AD,
      call_GQ,
      call_PGT,
      call_PID,
      call_PL
  FROM OPENROWSET(
              BULK '${file}',
              DATA_SOURCE = 'MyAzureBlobStorage',
              FORMAT ='CSV',
              FORMATFILE='vets_format.txt',
              FORMATFILE_DATA_SOURCE = 'MyAzureBlobStorage'
  ) as DataFile;
FIN
        else
            cat >> loader_script.sql <<FIN
  INSERT INTO ref_ranges with (TABLOCK)
  SELECT 
      CONVERT(BIGINT, location) AS location,
      CONVERT(INT, sample_id) AS sample_id,
      CONVERT(SMALLINT, length) AS length,
      CONVERT(TINYINT, state)
  FROM OPENROWSET(
              BULK '${file}',
              DATA_SOURCE = 'MyAzureBlobStorage',
              FORMAT ='CSV',
              FORMATFILE='ref_ranges_format.txt',
              FORMATFILE_DATA_SOURCE = 'MyAzureBlobStorage'
  ) as DataFile;
FIN
        fi
    done
done

```

# Costs

Described [here](https://learn.microsoft.com/en-au/azure/azure-sql/database/serverless-tier-overview?view=azuresql&tabs=general-purpose#billing).