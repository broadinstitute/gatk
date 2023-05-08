# Importing vet and ref_ranges data to Azure SQL Database, serverless style!

Similar to and building on the [earlier VS-879 spike](azure_sql_database_ingest.md), this spike imports vet and ref
ranges data from Azure Blob Storage to an Azure SQL Database. However this time the server will be... serverless and
costs will be measured.

# Loader script

See `load_azure_sql_database.sh` elsewhere in this repo. This script automates the following:

* Resource group creation
* Azure SQL Server creation
* Azure SQL Server firewall rule creation
* Serverless Azure SQL Database creation
* Blob storage connectivity
* `vets` and `ref_ranges` table creation
* `vets` and `ref_ranges` data loading from blob storage

# Costs

Described [here](https://learn.microsoft.com/en-au/azure/azure-sql/database/serverless-tier-overview?view=azuresql&tabs=general-purpose#billing).

Actual vCore seconds consumed during the load was 8.22K, instance priced at $0.000145 per vCore second so

```
8220 * $0.000145 = $1.19
```

or about $0.12 per sample. This is without doing any sort of indexing which would have increased costs, and it is
waiting a full hour for the database to autopause during which it was using 0.5 vCore. The wait for timeout cost

```
1 hour * 60 min / hour * 60 sec / min = 3600 sec * $0.000145 / vCore second * 0.5 vCore = $0.26
```

This "wait for autopause" should be a constant term regardless of how many samples are being loaded and should therefore
decrease as the number of samples being loaded increases.
