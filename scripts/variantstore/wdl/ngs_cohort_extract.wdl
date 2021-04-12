version 1.0

workflow WgsCohortExtract {

  input {
    File participant_ids
    String query_project
    String wgs_dataset
    String wgs_extraction_cohorts_dataset
    String wgs_extraction_destination_dataset
    String wgs_extraction_temp_tables_dataset
    String extraction_uuid
    String output_gcs_dir
    
    # Extract parameters
    File wgs_intervals
    Int scatter_count
       
    File reference
    File reference_index
    File reference_dict
        
    String output_file_base_name
    
    String? fq_filter_set_table
    String? filter_set_name
    File? gatk_override
  }

  call CreateCohortSampleTable {
    input:
      participant_ids = participant_ids,
      wgs_dataset = wgs_dataset,
      wgs_extraction_cohorts_dataset = wgs_extraction_cohorts_dataset,
      query_project = query_project,
      table_name = extraction_uuid
  }

  call CreateCohortExtractTable {
    input:
      cohort_uuid = extraction_uuid,
      fq_cohort_sample_table = CreateCohortSampleTable.fq_cohort_sample_table,
      query_project = query_project,
      wgs_dataset = wgs_dataset,
      wgs_extraction_cohorts_dataset = wgs_extraction_cohorts_dataset,
      wgs_extraction_destination_dataset = wgs_extraction_destination_dataset,
      wgs_extraction_temp_tables_dataset = wgs_extraction_temp_tables_dataset
  }
  
  call SplitIntervals {
    input:
        intervals = wgs_intervals,
        ref_fasta = reference,
        ref_fai = reference_index,
        ref_dict = reference_dict,
        scatter_count = scatter_count
  }
    
  scatter(i in range(scatter_count) ) {
    call ExtractTask {
      input:
        gatk_override            = gatk_override,
        reference                = reference,
        reference_index          = reference_index,
        reference_dict           = reference_dict,
        fq_sample_table          = CreateCohortSampleTable.fq_cohort_sample_table,
        intervals                = SplitIntervals.interval_files[i],
        fq_cohort_extract_table  = CreateCohortExtractTable.fq_cohort_extract_table,
        read_project_id          = query_project,
        fq_filter_set_table      = fq_filter_set_table,
        filter_set_name          = filter_set_name,
        output_file              = "${output_file_base_name}_${i}.vcf.gz",
        output_gcs_dir           = output_gcs_dir
    }
  }

  output {
    String fq_cohort_extract_table = CreateCohortExtractTable.fq_cohort_extract_table
  }

}

task CreateCohortSampleTable {

  input {
    File participant_ids
    String wgs_dataset
    String wgs_extraction_cohorts_dataset
    String query_project
    String table_name
  }
  
  command <<<
    echo "SELECT
      sample_id,
      sample_name
    FROM
      \`~{wgs_dataset}.metadata\`
    WHERE
      sample_name IN " > create_cohort.sql

    PARTICIPANT_IDS=$(cat ~{participant_ids} | awk '{print "\""$0"\""}' | paste -sd ",")
    echo "($PARTICIPANT_IDS)" >> create_cohort.sql
    
    DESTINATION_DATASET=$(echo ~{wgs_extraction_cohorts_dataset} | tr '.' ':')

    bq query \
      --project_id ~{query_project} \
      --destination_table ${DESTINATION_DATASET}.~{table_name} \
      --use_legacy_sql=false \
      --max_rows=10000000 \
      --allow_large_results < create_cohort.sql
    
    echo "~{wgs_extraction_cohorts_dataset}.~{table_name}" > fq_cohort_sample_table.txt
  >>>
  
  output {
    String fq_cohort_sample_table = read_string("fq_cohort_sample_table.txt")
  }

  runtime {
    memory: "3.75 GB"
    bootDiskSizeGb: "15"
    disks: "local-disk 10 HDD"
    preemptible: 3
    docker: "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }
} 


task CreateCohortExtractTable {

  input {
    String cohort_uuid
    String fq_cohort_sample_table
    String wgs_dataset
    String wgs_extraction_cohorts_dataset
    String wgs_extraction_destination_dataset
    String wgs_extraction_temp_tables_dataset
    String query_project
  }
  
  command <<<
    set -e
    
    echo "Exporting to ~{wgs_extraction_destination_dataset}.~{cohort_uuid}"

    python /app/ngs_cohort_extract.py \
      --fq_petvet_dataset ~{wgs_dataset} \
      --fq_temp_table_dataset ~{wgs_extraction_temp_tables_dataset} \
      --fq_destination_dataset ~{wgs_extraction_destination_dataset} \
      --destination_table ~{cohort_uuid} \
      --fq_cohort_sample_names ~{fq_cohort_sample_table} \
      --min_variant_samples 0 \
      --query_project ~{query_project} \
      --fq_sample_mapping_table ~{wgs_dataset}.metadata
    
    echo ~{wgs_extraction_destination_dataset}.~{cohort_uuid} > fq_cohort_extract_table.txt
  >>>
  
  output {
    String fq_cohort_extract_table = read_string("fq_cohort_extract_table.txt")
  }

  runtime {
    memory: "3.75 GB"
    bootDiskSizeGb: "15"
    disks: "local-disk 10 HDD"
    preemptible: 3
    docker: "gcr.io/all-of-us-workbench-test/variantstore-extract-prep:2"
  }
} 

task ExtractTask {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    input {
        # ------------------------------------------------
        # Input args:
        File reference
        File reference_index
        File reference_dict
    
        String fq_sample_table

        File intervals

        String fq_cohort_extract_table
        String read_project_id
        String output_file
        String output_gcs_dir
        
        String? fq_filter_set_table
        String? filter_set_name
        
        # Runtime Options:
        File? gatk_override
        
        Int? local_sort_max_records_in_ram = 5000000
    }
    
    String outdir = sub(output_gcs_dir, "/$", "")


    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        #bq query --project_id ~{read_project_id} --max_rows=10000000 --format=csv --nouse_legacy_sql \
        #  'select sample_id, sample_name FROM `~{fq_sample_table}`' | sed -e 1d \
        #  > sample_map.csv

        gatk --java-options "-Xmx9g" \
            ExtractCohort \
                --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
                -R "~{reference}" \
                -O "~{output_file}" \
                --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
                --sample-table ~{fq_sample_table} \
                --cohort-extract-table ~{fq_cohort_extract_table} \
                -L ~{intervals} \
                --project-id ~{read_project_id} \
                ~{"--variant-filter-table " + fq_filter_set_table} \
                ~{"--filter-set-name " + filter_set_name}
        
        gsutil cp ~{output_file} ~{output_gcs_dir}/
        gsutil cp ~{output_file}.tbi ~{output_gcs_dir}/
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-gatk/gatk:4.1.8.0"
        memory: "10 GB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 15
        preemptible: 3
        cpu: 2
    }

    # ------------------------------------------------
    # Outputs:
    output {
        String sample_csv = "sample_map.csv"
        File output_vcf = "~{output_file}"
        File output_vcf_index = "~{output_file}.tbi"
    }
 }

 task SplitIntervals {
     input {
       File? intervals
       File ref_fasta
       File ref_fai
       File ref_dict
       Int scatter_count
       String? split_intervals_extra_args

       File? gatk_override
     }

     command {
         set -e
         export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

         mkdir interval-files
         gatk --java-options "-Xmx5g" SplitIntervals \
             -R ~{ref_fasta} \
             ~{"-L " + intervals} \
             -scatter ~{scatter_count} \
             -O interval-files \
             ~{split_intervals_extra_args}
         cp interval-files/*.interval_list .
     }

     runtime {
         docker: "us.gcr.io/broad-gatk/gatk:4.1.8.0"
         bootDiskSizeGb: 15
         memory: "6 GB"
         disks: "local-disk 10 HDD"
         preemptible: 3
         cpu: 1
     }

     output {
         Array[File] interval_files = glob("*.interval_list")
     }
 }

