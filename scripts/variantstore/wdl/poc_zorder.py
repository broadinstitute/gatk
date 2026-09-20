from pyspark.sql import SparkSession

# --- CONFIGURATION ---
PROJECT_ID = "gvs-internal"
DATASET_ID = "NHGRI_AnVIL_ZERO"
BASELINE_TABLE = "poc_baseline_alt_allele"
ICEBERG_TABLE_NAME = "poc_iceberg_zordered"

# Initialize Spark Session (Configurations are passed via the submit command later)
spark = SparkSession.builder \
    .appName("BigQuery-to-Iceberg-ZOrder") \
    .getOrCreate()

# 1. Read the baseline data from BigQuery
print("Reading baseline data from BigQuery...")
df = spark.read \
    .format("bigquery") \
    .option("table", f"{PROJECT_ID}.{DATASET_ID}.{BASELINE_TABLE}") \
    .load()

# --- NEW: Create the namespace in the Iceberg catalog ---
print(f"Creating namespace local_catalog.{DATASET_ID}...")
spark.sql(f"CREATE NAMESPACE IF NOT EXISTS local_catalog.{DATASET_ID}")

iceberg_target = f"local_catalog.{DATASET_ID}.{ICEBERG_TABLE_NAME}"

# 2. Write data to GCS as an Iceberg table using the V2 API
print(f"Writing data to Iceberg table: {iceberg_target}...")
df.writeTo(iceberg_target) \
    .createOrReplace()

# 3. Execute the Z-Order compaction procedure
print("Executing Z-Order compaction on (location, sample_id)...")
spark.sql(f"""
    CALL local_catalog.system.rewrite_data_files(
        table => '{DATASET_ID}.{ICEBERG_TABLE_NAME}', 
        strategy => 'sort', 
        sort_order => 'zorder(location, sample_id)'
    )
""")

print("Job completed successfully!")
spark.stop()