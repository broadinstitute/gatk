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

All trials below with Family = General, Edition = Gen5 resulting in compute priced at $0.000145 per vCore second in all
cases. Cost per vCore second is visible at

```
Portal -> SQL Database -> <specific server> -> Compute + storage
```

![Cost per vCore second](./cost%20per%20vCore%20second.png)

Actual vCore * second usage is visible at
```
Portal -> SQL Database -> <specific server> -> Metrics
```

![Actual vCore * second usage](./App%20CPU%20billed.png)

using the metric "App CPU billed".


## Trial 1: Capacity = 2

Actual vCore seconds consumed during the load was 8.22K so

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

## Trial 2: Capacity = 8

Instance again priced at $0.000145 per vCore second.

```
22K * $0.000145 = $3.19
```

or about $0.32 per sample, clearly worse than the first trial. One obvious problem with this configuration is that at a
Capacity of 8 the vCore minimum is 1 versus the vCore minimum of 0.5 with the default Capacity of 2. However it seems
from the metrics that we're consistently above 50% of the log rate percentage and log rates appear to
be [a function of Capacity / max vCores](https://techcommunity.microsoft.com/t5/azure-sql-blog/raising-log-rate-limits-for-general-purpose-service-tier-in/ba-p/1784622).
So this 4x Capacity should have allowed us to write significantly faster than we could with a Capacity of 2, but we
don't seem to be able to use all of our newfound log rate bandwidth. Repeat the experiment with a Capacity of 4 to see
if we can get the benefits of improved log rates with a lower minimum vCore.

## Trial 3: Capacity = 4

```
12.19K * $0.000145 = $1.77
```

So $1.77 per sample. Not as good as the default.

## Trial 4: Capacity = 1

```
6.5K * $0.000145 = $0.94
```

or around $0.09 per sample. This configuration was clearly the cost winner for this toy 10 sample Quickstart dataset,
but real-world sample set sizes may require the use of a larger instance to load within a reasonable timeframe.
