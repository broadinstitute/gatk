import argparse
from pyspark.sql import SparkSession

def main():
    parser = argparse.ArgumentParser(description="Rewrite Parquet files into a partitioned and clustered Iceberg table.")
    parser.add_argument("--parquet_path", required=True, help="GCS path to the raw Parquet files")
    parser.add_argument("--table_name", required=True, help="Full Iceberg table identifier (catalog.namespace.table)")

    # New optional arguments for Parquet tuning
    parser.add_argument("--row_group_size_bytes", required=False, default=None,
                        help="Optional: Set the Parquet row group size in bytes (e.g., 16777216 for 16MB)")
    parser.add_argument("--page_size_bytes", required=False, default=None,
                        help="Optional: Set the Parquet page size in bytes (e.g., 1048576 for 1MB)")

    args = parser.parse_args()

    parquet_path = args.parquet_path
    table_name = args.table_name
    row_group_size = args.row_group_size_bytes
    page_size = args.page_size_bytes

    spark = SparkSession.builder.appName(f"ParquetToIcebergRewrite_{table_name}").getOrCreate()

    print(f"Starting rewrite process for table: {table_name}")
    print(f"Reading full dataset from: {parquet_path}")

    # 1. Read the ACTUAL data
    df = spark.read.parquet(parquet_path)

    print(f"Creating empty partitioned Iceberg table: {table_name}")

    # 2. Create the empty Iceberg table
    df.limit(0).writeTo(table_name).partitionedBy("sample_id").create()

    # 3. Apply custom Parquet properties if the user provided them
    tbl_properties = []
    if row_group_size:
        tbl_properties.append(f"'write.parquet.row-group-size-bytes'='{row_group_size}'")
    if page_size:
        tbl_properties.append(f"'write.parquet.page-size-bytes'='{page_size}'")

    if tbl_properties:
        properties_sql = ", ".join(tbl_properties)
        print(f"Applying custom Parquet block sizes: {properties_sql}")
        spark.sql(f"ALTER TABLE {table_name} SET TBLPROPERTIES ({properties_sql})")

    # 4. Apply the clustering rule
    if "location" in df.columns:
        print(f"Setting strict write order to 'location' for {table_name}...")
        spark.sql(f"ALTER TABLE {table_name} WRITE ORDERED BY location")

    print(f"Rewriting and clustering data into Iceberg warehouse... (This may take a while)")

    # 5. Append the data using the new properties and sort order
    df.writeTo(table_name).append()

    print(f"Iceberg rewrite and clustering complete for {table_name}!")

if __name__ == "__main__":
    main()