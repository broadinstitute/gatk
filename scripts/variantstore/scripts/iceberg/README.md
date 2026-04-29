# GVS Parquet to Iceberg BigLake Pipeline

This guide outlines the steps to ingest raw GVS Parquet files (`vet`, `ref_ranges`, `sample_chromosome_ploidy`) into an Apache Iceberg warehouse on GCS, and expose them as highly optimized external tables in Google BigQuery.

The ingestion process applies specific partitioning (`sample_id`), physical clustering (linear sort by `location`), and custom Parquet block sizing (16MB row groups) to enable millisecond, byte-level file pruning in BigQuery.

## Prerequisites
1. **Google Cloud Storage (GCS):** A bucket to host the raw Parquet and the Iceberg warehouse.
2. **Dataproc Serverless:** API enabled and a valid subnet configured.
3. **BigLake Connection:** A pre-existing BigQuery external data connection (e.g., `iceberg-poc-connection-us`).

## Step 1: Set Environment Variables
Update these variables for your specific environment before running the commands below.

```bash
export GCP_REGION="us-central1"
export SUBNET_NAME="gvs-network-subnet"
export BUCKET="gs://iceberg-poc-2026-04-14"

export RAW_DATA_PATH="${BUCKET}/raw_parquet"
export ICEBERG_WAREHOUSE="${BUCKET}/iceberg-warehouse"
export CATALOG_NAME="iceberg_catalog_2026_04_14"
export NAMESPACE="iceberg_2026_04_14"

# Spark properties for Iceberg Hadoop Catalog
export SPARK_PROPS="spark.jars.packages=org.apache.iceberg:iceberg-spark-runtime-3.5_2.13:1.4.3,\
spark.sql.extensions=org.apache.iceberg.spark.extensions.IcebergSparkSessionExtensions,\
spark.sql.catalog.${CATALOG_NAME}=org.apache.iceberg.spark.SparkCatalog,\
spark.sql.catalog.${CATALOG_NAME}.type=hadoop,\
spark.sql.catalog.${CATALOG_NAME}.warehouse=${ICEBERG_WAREHOUSE}/"
```

## Step 2: Rewrite & Cluster Parquet to Iceberg

Submit a Dataproc Serverless batch job for each table type using the `register_iceberg_rewrite.py` script. This script creates the table, applies the `sample_id` partitioning, sets the `location` sort order, and tuning the Parquet block sizes for optimal BigQuery pruning.

### 1. Vet Table
```bash
gcloud dataproc batches submit pyspark register_iceberg_rewrite.py \
    --region=${GCP_REGION} \
    --version="2.2" \
    --subnet=${SUBNET_NAME} \
    --properties=${SPARK_PROPS} \
    -- \
    --parquet_path="${RAW_DATA_PATH}/vet/001/" \
    --table_name="${CATALOG_NAME}.${NAMESPACE}.vet_001_iceberg" \
    --row_group_size_bytes="16777216" \
    --page_size_bytes="1048576"
```

### 2. Ref Ranges Table
```bash
gcloud dataproc batches submit pyspark register_iceberg_rewrite.py \
    --region=${GCP_REGION} \
    --version="2.2" \
    --subnet=${SUBNET_NAME} \
    --properties=${SPARK_PROPS} \
    -- \
    --parquet_path="${RAW_DATA_PATH}/ref_ranges/001/" \
    --table_name="${CATALOG_NAME}.${NAMESPACE}.ref_ranges_001_iceberg" \
    --row_group_size_bytes="16777216"
```

### 3. Sample Chromosome Ploidy Table
```bash
gcloud dataproc batches submit pyspark register_iceberg_rewrite.py \
    --region=${GCP_REGION} \
    --version="2.2" \
    --subnet=${SUBNET_NAME} \
    --properties=${SPARK_PROPS} \
    -- \
    --parquet_path="${RAW_DATA_PATH}/sample_chromosome_ploidy/001/" \
    --table_name="${CATALOG_NAME}.${NAMESPACE}.sample_chromosome_ploidy_iceberg" \
    --row_group_size_bytes="16777216"
```

## Step 3: Automate BigQuery External Table Creation

This script automatically finds the highest versioned `.metadata.json` file for each table in GCS and uses the `bq` command-line tool to create the BigLake external tables.

```bash
# 1. Define your BigQuery project, dataset, and connection
export BQ_TARGET="your_project.your_dataset"
export BQ_CONNECTION="gvs-internal.us.iceberg-poc-connection-us"

# 2. Define a helper function to find the latest metadata file
get_latest_metadata() {
    local table_dir=$1
    # Lists the files, uses version sort (-V) to handle v2 vs v10 correctly, and grabs the last one
    gcloud storage ls "${ICEBERG_WAREHOUSE}/${NAMESPACE}/${table_dir}/metadata/v*.metadata.json" | sort -V | tail -n 1
}

# 3. Fetch the latest URIs
echo "Finding latest metadata files..."
VET_URI=$(get_latest_metadata "vet_001_iceberg")
REF_URI=$(get_latest_metadata "ref_ranges_001_iceberg")
PLOIDY_URI=$(get_latest_metadata "sample_chromosome_ploidy_iceberg")

echo "Vet URI: ${VET_URI}"
echo "Ref URI: ${REF_URI}"
echo "Ploidy URI: ${PLOIDY_URI}"

# 4. Execute the DDL in BigQuery
echo "Creating BigQuery external tables..."

bq query --nouse_legacy_sql "
CREATE OR REPLACE EXTERNAL TABLE \`${BQ_TARGET}.vet_001_iceberg\`
WITH CONNECTION \`${BQ_CONNECTION}\`
OPTIONS (format = 'ICEBERG', uris = ['${VET_URI}']);
"

bq query --nouse_legacy_sql "
CREATE OR REPLACE EXTERNAL TABLE \`${BQ_TARGET}.ref_ranges_001_iceberg\`
WITH CONNECTION \`${BQ_CONNECTION}\`
OPTIONS (format = 'ICEBERG', uris = ['${REF_URI}']);
"

bq query --nouse_legacy_sql "
CREATE OR REPLACE EXTERNAL TABLE \`${BQ_TARGET}.sample_chromosome_ploidy_iceberg\`
WITH CONNECTION \`${BQ_CONNECTION}\`
OPTIONS (format = 'ICEBERG', uris = ['${PLOIDY_URI}']);
"

echo "Successfully linked Iceberg tables to BigQuery!"
```

## Verification
To verify that Iceberg file pruning and Parquet block-level pruning are working, run a needle-in-a-haystack query with the BigQuery cache **disabled**:

```sql
SELECT SUM(MOD(location, 2))
FROM `your_project.your_dataset.vet_001_iceberg`
WHERE sample_id = 3
  AND location BETWEEN 1000000000000 AND 1000000000010;
```
Check the **Job Information -> Bytes processed**. It should be in the low megabytes (~1-3 MB) rather than scanning the full partition.