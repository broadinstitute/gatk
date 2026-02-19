#!/usr/bin/env python3
"""
Create the BigQuery tracking table for Parquet load status.
"""

import argparse
import sys

from google.cloud import bigquery


def create_tracking_table(project_id, dataset_name):
    """
    Create the parquet_load_status tracking table if it doesn't exist.
    
    Args:
        project_id: BigQuery project ID
        dataset_name: BigQuery dataset name
    """
    client = bigquery.Client(project=project_id)
    table_id = f"{project_id}.{dataset_name}.parquet_load_status"
    
    print(f"Creating tracking table: {table_id}")
    
    schema = [
        bigquery.SchemaField("file_path", "STRING", mode="REQUIRED"),
        bigquery.SchemaField("table_name", "STRING", mode="REQUIRED"),
        bigquery.SchemaField("sample_id", "INT64", mode="REQUIRED"),
        bigquery.SchemaField("load_timestamp", "TIMESTAMP", mode="REQUIRED"),
        bigquery.SchemaField("load_job_id", "STRING", mode="REQUIRED"),
        bigquery.SchemaField("file_size_bytes", "INT64", mode="NULLABLE"),
        bigquery.SchemaField("rows_loaded", "INT64", mode="NULLABLE"),
    ]
    
    table = bigquery.Table(table_id, schema=schema)
    
    # Note: BigQuery doesn't enforce PRIMARY KEY, but we document it
    table.description = (
        "Tracks which Parquet files have been loaded into BigQuery tables. "
        "file_path is the primary key (enforced via MERGE logic, not by BigQuery)."
    )
    
    try:
        table = client.create_table(table, exists_ok=True)
        print(f"✓ Table {table_id} created successfully")
        print(f"  Schema: {len(schema)} columns")
        return True
    except Exception as e:
        print(f"✗ Error creating table: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Create BigQuery tracking table for Parquet loads"
    )
    parser.add_argument("--project-id", required=True, help="BigQuery project ID")
    parser.add_argument("--dataset-name", required=True, help="BigQuery dataset name")
    
    args = parser.parse_args()
    
    success = create_tracking_table(args.project_id, args.dataset_name)
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
