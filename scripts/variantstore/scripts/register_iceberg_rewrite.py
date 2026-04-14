import argparse
from pyspark.sql import SparkSession

def main():
    parser = argparse.ArgumentParser(description="Rewrite Parquet files into a partitioned and clustered Iceberg table.")
    parser.add_argument("--parquet_path", required=True, help="GCS path to the raw Parquet files")
    parser.add_argument("--table_name", required=True, help="Full Iceberg table identifier (catalog.namespace.table)")
    args = parser.parse_args()

    parquet_path = args.parquet_path
    table_name = args.table_name

    spark = SparkSession.builder.appName(f"ParquetToIcebergRewrite_{table_name}").getOrCreate()

    print(f"Starting rewrite process for table: {table_name}")
    print(f"Reading full dataset from: {parquet_path}")

    # 1. Read the ACTUAL data (Notice we removed .limit(0))
    # Spark will still pick up the 'sample_id' partition from your folder structure
    df = spark.read.parquet(parquet_path)

    print(f"Creating empty partitioned Iceberg table: {table_name}")

    # 2. Create the empty Iceberg table first using just the schema
    df.limit(0).writeTo(table_name).partitionedBy("sample_id").create()

    # 3. Apply the clustering rule to the table
    if "location" in df.columns:
        print(f"Setting strict write order to 'location' for {table_name}...")
        spark.sql(f"ALTER TABLE {table_name} WRITE ORDERED BY location")

    print(f"Rewriting and clustering data into Iceberg warehouse... (This may take a while)")

    # 4. Append the data.
    # Because of the ALTER TABLE command above, Iceberg will automatically
    # shuffle and sort the data by 'location' as it writes the new files!
    df.writeTo(table_name).append()

    print(f"Iceberg rewrite and clustering complete for {table_name}!")

if __name__ == "__main__":
    main()