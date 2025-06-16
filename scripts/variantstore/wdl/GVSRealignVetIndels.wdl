version 1.0

import "GvsUtils.wdl" as Utils

workflow GVSRealignVetIndels {
  input {
    String dataset_name
    String project_id
    Int vet_table_number
    File? gatk_override
    
    # Reference genome parameters
    String reference_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    
    # VET table configuration
    Int samples_per_vet_table = 4000
    Int? max_records_per_sample = 1000000
    
    # Optional runtime parameters
    String? gatk_docker
    String? git_branch_or_tag
    String? git_hash
  }

  # Get tool versions if not provided
  if (!defined(gatk_docker) || !defined(git_hash)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  # Input validation notes:
  # - vet_table_number should be >= 1
  # - samples_per_vet_table should be >= 1  
  # - reference_fasta should be a valid FASTA file path

  # Compute VET table name with zero-padded 3-digit format
  String vet_table_name = if (vet_table_number < 10) then "vet_00" + vet_table_number else if (vet_table_number < 100) then "vet_0" + vet_table_number else "vet_" + vet_table_number
  
  # Compute sample ID range for this VET table
  # Formula: start = (vet_table_number - 1) * samples_per_vet_table + 1
  # Range covers samples_per_vet_table samples: [start, start + samples_per_vet_table - 1] inclusive
  Int sample_id_start = (vet_table_number - 1) * samples_per_vet_table + 1
  Int sample_id_end = sample_id_start + samples_per_vet_table - 1
  
  # Create array of sample IDs to scatter over
  Array[Int] sample_ids = range(sample_id_end - sample_id_start + 1)
  Array[Int] actual_sample_ids = sample_ids # Add sample_id_start to each element
  
  # Scatter over each sample ID in the VET table range
  scatter (i in range(length(sample_ids))) {
    Int current_sample_id = sample_id_start + i
    
    call RealignSingleSample {
      input:
        dataset_name = dataset_name,
        project_id = project_id,
        vet_table_name = vet_table_name,
        sample_id = current_sample_id,
        reference_fasta = reference_fasta,
        max_records = select_first([max_records_per_sample, 1000000]),
        gatk_docker = effective_gatk_docker,
        gatk_override = gatk_override
    }
  }

  output {
    # Collect all CSV outputs from the scattered tasks
    Array[File] realigned_variants_csvs = RealignSingleSample.realigned_variants_csv
    Array[File] position_bucket_histograms = RealignSingleSample.position_bucket_histogram_csv
    Array[File] shifted_indels_per_sample_csvs = RealignSingleSample.shifted_indels_per_sample_csv
    Array[File] indel_realignment_summaries = RealignSingleSample.indel_realignment_summary_csv
    Array[File] monitoring_logs = RealignSingleSample.monitoring_log
    String recorded_git_hash = effective_git_hash
    Boolean done = true
  }
}

################################################################################

task RealignSingleSample {
  input {
    String dataset_name
    String project_id
    String vet_table_name
    Int sample_id
    String reference_fasta
    Int max_records
    
    # Runtime parameters
    String gatk_docker
    File? gatk_override
    
    # Resource parameters
    Int? memory_gb = 8
    Int? disk_gb = 500
    Int? cpu_count = 3
    Int? preemptible_tries = 3
  }

  meta {
    description: "Realign indels for a single sample using VetIndelRealigner tool"
  }

  parameter_meta {
    dataset_name: "Name of the dataset containing the VET table"
    project_id: "Google Cloud project ID containing the BigQuery dataset"
    vet_table_name: "Name of the VET table (e.g., vet_007)"
    sample_id: "Sample ID to process for indel realignment"
    reference_fasta: "Path to reference genome FASTA file"
    max_records: "Maximum number of records to process for this sample"
  }

  # Construct fully qualified VET table name
  String fq_vet_table = "~{project_id}.~{dataset_name}.~{vet_table_name}"
  
  # Sample ID filter string
  String sample_filter = "sample_id IN (~{sample_id})"
  
  # Output file prefix
  String output_prefix = "sample_~{sample_id}_~{vet_table_name}"
  
  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    # Prepend date, time and pwd to xtrace log entries
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # Start monitoring in background
    bash ~{monitoring_script} > monitoring.log &

    # Set GATK JAR if override provided
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    # Display system information
    df -h
    free -h
    
    echo "Processing sample ~{sample_id} from VET table ~{fq_vet_table}"
    echo "Using sample filter: ~{sample_filter}"
    
    # Run VetIndelRealigner
    gatk --java-options "-Xmx~{memory_gb - 1}g" VetIndelRealigner \
      --project-id ~{project_id} \
      --vet-table ~{fq_vet_table} \
      -R ~{reference_fasta} \
      --sample-id-filter "~{sample_filter}" \
      --log-cadence 50 \
      --max-records ~{max_records}
    
    # Rename output files to include sample and table information
    if [ -f "realigned_variants.csv" ]; then
      mv realigned_variants.csv ~{output_prefix}_realigned_variants.csv
    fi
    
    if [ -f "position_bucket_histogram.csv" ]; then
      mv position_bucket_histogram.csv ~{output_prefix}_position_bucket_histogram.csv
    fi
    
    if [ -f "shifted_indels_per_sample.csv" ]; then
      mv shifted_indels_per_sample.csv ~{output_prefix}_shifted_indels_per_sample.csv
    fi
    
    if [ -f "indel_realignment_summary.csv" ]; then
      mv indel_realignment_summary.csv ~{output_prefix}_indel_realignment_summary.csv
    fi
    
    echo "VetIndelRealigner completed successfully for sample ~{sample_id}"
  >>>

  runtime {
    docker: gatk_docker
    memory: "~{memory_gb} GB"
    disks: "local-disk ~{disk_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: preemptible_tries
    maxRetries: 1
    cpu: cpu_count
    noAddress: true
  }

  output {
    File realigned_variants_csv = "~{output_prefix}_realigned_variants.csv"
    File position_bucket_histogram_csv = "~{output_prefix}_position_bucket_histogram.csv"
    File shifted_indels_per_sample_csv = "~{output_prefix}_shifted_indels_per_sample.csv"
    File indel_realignment_summary_csv = "~{output_prefix}_indel_realignment_summary.csv"
    File monitoring_log = "monitoring.log"
  }
}
