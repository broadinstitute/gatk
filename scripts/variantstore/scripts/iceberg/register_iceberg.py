import argparse
from pyspark.sql import SparkSession

def main():
    # 1. Set up argument parsing to accept parameters
    parser = argparse.ArgumentParser(description="Register Parquet files to an Iceberg table in-place.")
    parser.add_argument("--parquet_path", required=True, help="GCS path to the raw Parquet files")
    parser.add_argument("--table_name", required=True, help="Name of the Iceberg table to create")
    args = parser.parse_args()

    parquet_path = args.parquet_path
    table_name = args.table_name

    # 2. Initialize the Spark Session (Naming the app dynamically based on the table)
    spark = SparkSession.builder.appName(f"ParquetToIceberg_{table_name}").getOrCreate()

    print(f"Starting registration for table: {table_name}")
    print(f"Reading schema from: {parquet_path}")

    # 3. Read the schema
    df = spark.read.parquet(parquet_path).limit(0)

    print(f"Creating empty Iceberg table: {table_name}")

    # 4. Create the empty Iceberg table
    df.writeTo(table_name).create()

    # Extract the catalog name from the first part of the table_name
    # e.g., splits "iceberg_catalog_2026_04_13.db.table" and grabs "iceberg_catalog_2026_04_13"
    catalog_name = table_name.split('.')[0]

    print(f"Registering Parquet files to Iceberg metadata using catalog: {catalog_name}...")

    # Use the dynamic catalog_name in the CALL procedure
    spark.sql(f"""
        CALL {catalog_name}.system.add_files(
            table => '{table_name}',
            source_table => '`parquet`.`{parquet_path}`'
        )
    """)

    print(f"Iceberg registration complete for {table_name}!")

if __name__ == "__main__":
    main()