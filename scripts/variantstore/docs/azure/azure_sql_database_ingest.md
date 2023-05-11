# Importing vet and ref_ranges data to Azure SQL Database

## Extract CSVs from BigQuery to Google Cloud Storage

### Variant data

```

EXPORT DATA OPTIONS(
                uri='gs://bucket/path/to/vets/vet_001/vet_001_*.csv', format='CSV', header=false) AS
                SELECT v.sample_id as sample_id, location, ref, alt, AS_RAW_MQ, AS_RAW_MQRankSum, QUALapprox, AS_QUALapprox, AS_RAW_ReadPosRankSum, AS_SB_TABLE, AS_VarDP, call_GT, call_AD, call_GQ, call_PGT, call_PID, call_PL
                FROM `gvs-internal.quickit_dataset.vet_001` v
                INNER JOIN `gvs-internal.quickit_dataset.sample_info` s ON s.sample_id = v.sample_id
                WHERE withdrawn IS NULL AND
                is_control = false
                ORDER BY location
```


### Reference data

```

EXPORT DATA OPTIONS (
                uri='gs://bucket/path/to/ref_ranges/ref_ranges_001/ref_ranges_001_*.csv', format='CSV', header=false) AS
                SELECT r.sample_id as sample_id, location, length, state
                FROM `gvs-internal.quickit_dataset.ref_ranges_001` r
                INNER JOIN `gvs-internal.quickit_dataset.sample_info` s ON s.sample_id = r.sample_id
                WHERE withdrawn IS NULL AND
                is_control = false
                ORDER BY location
```


## Transfer CSVs from GCP to Azure

Create a storage container in Azure, use `azcopy` to transfer files. Example below for vets, reference ranges are similar:


```
azcopy cp "https://storage.cloud.google.com/bucket/path/to/vets/" 'https://mystorageaccount.blob.core.windows.net/container?<SAS TOKEN>' --recursive
```

##  Configure Azure SQL Database to read from Azure Blob Storage

### Get a database token

```
$ get_db_token() {
    # Fetches a database access token with ~1 hour TTL, ASCII / UTF-8 encoded and newline terminated
    # https://learn.microsoft.com/en-us/sql/connect/odbc/linux-mac/connecting-with-sqlcmd?view=azuresqldb-current
    az account get-access-token --resource https://database.windows.net --query accessToken --output tsv
}

$ get_db_token | tr -d '\n' | iconv -f ascii -t UTF-16LE > token.txt
```

### One-time Azure SQL Database Setup

Interactive connection via `sqlcmd`. This can be done locally or on an Azure VM. `sqlcmd` leaves a lot to be desired as
a SQL client, if we can find another way to connect to Azure SQL database with this Azure token we should look into that.

```
% sqlcmd -S tcp:${SQL_SERVER}.database.windows.net,1433 -d ${SQL_DATABASE} -G -P token.txt
```

Within `sqlcmd`
(https://learn.microsoft.com/en-us/sql/t-sql/functions/openrowset-transact-sql?view=sql-server-ver16):

```
CREATE MASTER KEY ENCRYPTION BY PASSWORD = '<Password>';
go
```

The SAS token referenced below should have at least `Read` and `List` on the Azure storage container that contains the
CSVs. Define a credential referencing the storage container and a datasource referencing the credential:

```
CREATE DATABASE SCOPED CREDENTIAL MyAzureBlobStorageCredential
WITH IDENTITY = 'SHARED ACCESS SIGNATURE',
SECRET = '<SAS Token>';
go

```

```
CREATE EXTERNAL DATA SOURCE MyAzureBlobStorage
WITH ( TYPE = BLOB_STORAGE,
LOCATION = '<Storage Container URL>',
CREDENTIAL= MyAzureBlobStorageCredential
);
go
```

## Create `ref_ranges` and `vets` tables

```
CREATE TABLE ref_ranges (
    location BIGINT,
    sample_id INT,
    length SMALLINT,
    state TINYINT
);
go

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
go
```

## Upload vets and ref_ranges format files

Upload the `ref_ranges_format.txt` and `vets_format.txt` files in `src/script/variantstore/azure/cosmos` to the storage
container referenced above.

## Generate script to import from each of the CSVs in Azure Blob Storage

In an Azure VM or locally:

```
$ read -r -d '' TEMPLATE << FIN
  INSERT INTO ref_ranges with (TABLOCK)
  SELECT 
      CONVERT(BIGINT, location) AS location,
      CONVERT(INT, sample_id) AS sample_id,
      CONVERT(SMALLINT, length) AS length,
      CONVERT(TINYINT, state)
  FROM OPENROWSET(
              BULK 'ref_ranges_001/ref_ranges_001_%012d.csv',
              DATA_SOURCE = 'MyAzureBlobStorage',
              FORMAT ='CSV',
              FORMATFILE='ref_ranges_format.txt',
              FORMATFILE_DATA_SOURCE = 'MyAzureBlobStorage'
  ) as DataFile;
  
FIN

# 52 is the magic number for the last ref_ranges CSV for a Quickstart integration run.
$ for i in $(seq 0 52)
do
  printf "${TEMPLATE}\n\n" $i >> insert_ref_ranges.sql
done

sqlcmd -S tcp:${SQL_SERVER}.database.windows.net,1433 -d ${SQL_DATABASE} -G -P token.txt -Q "$(cat insert_ref_ranges.sql)"
```


```
read -r -d '' TEMPLATE << FIN
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
              BULK 'vet_001/vet_001_%012d.csv',
              DATA_SOURCE = 'MyAzureBlobStorage',
              FORMAT ='CSV',
              FORMATFILE='vets_format.txt',
              FORMATFILE_DATA_SOURCE = 'MyAzureBlobStorage'
  ) as DataFile;

FIN

# 30 is the magic number for the last vets CSV for a Quickstart integration run.
$ for i in $(seq 0 30)
do
  printf "${TEMPLATE}\n\n" $i >> insert_vets.sql
done

$ sqlcmd -S tcp:${SQL_SERVER}.database.windows.net,1433 -d ${SQL_DATABASE} -G -P token.txt -Q "$(cat insert_vets.sql)"

```

Done!