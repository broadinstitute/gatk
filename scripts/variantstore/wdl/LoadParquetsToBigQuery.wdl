version 1.0

## LoadParquetsToBigQuery
##
## Loads Parquet files from GCS into BigQuery tables with fault tolerance and idempotency.
##
## This workflow:
## 1. Creates a tracking table to record loaded files
## 2. Discovers Parquet files in GCS and groups them by target BigQuery table
## 3. Loads files for each table in parallel with preemption recovery
## 4. Verifies all files were successfully loaded
##
## Design document: docs/parquet_to_bigquery_loading_design.md

workflow LoadParquetsToBigQuery {
  input {
    String output_gcs_dir
    String project_id
    String dataset_name
    Array[String] table_prefixes = ["vet", "ref_ranges"]
    String? billing_project_id
    Int max_parallel_loads = 50
    Int batch_size = 10000
    
    # Docker images
    String cloud_sdk_docker = "google/cloud-sdk:latest"
    String python_bq_docker = "gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }
  
  # Ensure tracking table exists
  call CreateTrackingTable {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      docker = python_bq_docker
  }
  
  # Discover all files and group by table
  call DiscoverParquetFiles {
    input:
      output_gcs_dir = output_gcs_dir,
      project_id = project_id,
      dataset_name = dataset_name,
      table_prefixes = table_prefixes,
      tracking_table_ready = CreateTrackingTable.done,
      billing_project_id = billing_project_id,
      docker = cloud_sdk_docker
  }
  
  # Load each table in parallel
  scatter (i in range(length(DiscoverParquetFiles.table_names))) {
    call LoadParquetFilesToBQ {
      input:
        project_id = project_id,
        dataset_name = dataset_name,
        table_name = DiscoverParquetFiles.table_names[i],
        files_to_load = DiscoverParquetFiles.file_fofns[i],
        schema_path = DiscoverParquetFiles.schema_paths[i],
        batch_size = batch_size,
        billing_project_id = billing_project_id,
        docker = python_bq_docker
    }
  }
  
  # Verify everything loaded successfully
  call VerifyAllLoaded {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      gcs_files_list = DiscoverParquetFiles.all_files_list,
      load_outputs = LoadParquetFilesToBQ.completion_status,
      docker = python_bq_docker
  }
  
  output {
    Boolean all_loaded = VerifyAllLoaded.all_loaded
    Int total_files = VerifyAllLoaded.total_files
    Int loaded_files = VerifyAllLoaded.loaded_files
    Int missing_files = VerifyAllLoaded.missing_files
    Array[String] tables_loaded = DiscoverParquetFiles.table_names
    File? missing_files_list = VerifyAllLoaded.missing_files_list
    File verification_results = VerifyAllLoaded.results_json
    Array[File] load_stats = LoadParquetFilesToBQ.stats_json
  }
}

task CreateTrackingTable {
  input {
    String project_id
    String dataset_name
    String docker
  }
  
  command <<<
    set -euo pipefail
    
    python3 /app/create_tracking_table.py \
      --project-id ~{project_id} \
      --dataset-name ~{dataset_name}
  >>>
  
  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
  
  output {
    Boolean done = true
  }
}

task DiscoverParquetFiles {
  input {
    String output_gcs_dir
    String project_id
    String dataset_name
    Array[String] table_prefixes
    Boolean tracking_table_ready
    String? billing_project_id
    String docker
  }
  
  meta {
    volatile: true
  }
  
  command <<<
    set -euo pipefail
    
    # List all objects, filter for parquet files
    echo "Listing files in ~{output_gcs_dir}..."
    gcloud storage ls --recursive ~{"--billing-project " + billing_project_id} \
      "~{output_gcs_dir}/" > all_objects.txt || true
    
    grep '\.parquet$' all_objects.txt > all_files.txt || touch all_files.txt
    
    FILE_COUNT=$(wc -l < all_files.txt)
    echo "Found $FILE_COUNT Parquet files"
    
    # Parse and group by table
    python3 /app/parse_and_group_files.py \
      --input all_files.txt \
      --output-dir grouped_files \
      --project-id ~{project_id} \
      --dataset ~{dataset_name} \
      --table-prefixes ~{sep=" " table_prefixes}
    
    # For each table, create a dummy schema file
    # In production, these would be actual schema files
    while IFS= read -r table_name; do
      echo '[]' > "grouped_files/${table_name}.schema.json"
    done < grouped_files/table_names.txt
    
    # Create list of schema paths
    while IFS= read -r table_name; do
      echo "$(pwd)/grouped_files/${table_name}.schema.json"
    done < grouped_files/table_names.txt > grouped_files/schema_paths.txt
  >>>
  
  runtime {
    docker: docker
    memory: "4 GB"
    disks: "local-disk 50 HDD"
    preemptible: 3
    cpu: 2
  }
  
  output {
    Array[String] table_names = read_lines("grouped_files/table_names.txt")
    Array[File] file_fofns = read_lines("grouped_files/fofn_paths.txt")
    Array[File] schema_paths = read_lines("grouped_files/schema_paths.txt")
    File all_files_list = "all_files.txt"
    File stats_json = "grouped_files/stats.json"
  }
}

task LoadParquetFilesToBQ {
  input {
    String project_id
    String dataset_name
    String table_name
    File files_to_load
    File schema_path
    Int batch_size
    String? billing_project_id
    String docker
  }
  
  command <<<
    set -euo pipefail
    
    python3 /app/load_parquet_to_bq.py \
      --project-id ~{project_id} \
      --dataset-name ~{dataset_name} \
      --table-name ~{table_name} \
      --files-fofn ~{files_to_load} \
      --schema-path ~{schema_path} \
      --pending-jobs-path pending_jobs.json \
      --batch-size ~{batch_size} \
      --output-stats stats.json
  >>>
  
  runtime {
    docker: docker
    memory: "4 GB"
    disks: "local-disk 20 HDD"
    preemptible: 5
    maxRetries: 3
    cpu: 1
  }
  
  output {
    String completion_status = read_string("stats.json")
    File stats_json = "stats.json"
  }
}

task VerifyAllLoaded {
  input {
    String project_id
    String dataset_name
    File gcs_files_list
    Array[String] load_outputs
    String docker
  }
  
  meta {
    volatile: true
  }
  
  command <<<
    set -euo pipefail
    
    mkdir -p verification_output
    
    python3 /app/verify_all_loaded.py \
      --project-id ~{project_id} \
      --dataset-name ~{dataset_name} \
      --gcs-files-list ~{gcs_files_list} \
      --output-dir verification_output
  >>>
  
  runtime {
    docker: docker
    memory: "4 GB"
    disks: "local-disk 20 HDD"
    cpu: 1
  }
  
  output {
    File results_json = "verification_output/verification_results.json"
    Boolean all_loaded = read_json(results_json)["all_loaded"]
    Int total_files = read_json(results_json)["total_files"]
    Int loaded_files = read_json(results_json)["loaded_files"]
    Int missing_files = read_json(results_json)["missing_files"]
    File? missing_files_list = "verification_output/missing_files.txt"
  }
}
