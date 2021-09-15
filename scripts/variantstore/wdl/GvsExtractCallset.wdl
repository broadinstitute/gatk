version 1.0

workflow GvsExtractCallset {
   input {
        String data_project
        String default_dataset

        File wgs_intervals
        Int scatter_count

        File reference
        File reference_index
        File reference_dict

        String fq_cohort_extract_table_prefix
        String query_project = data_project

        Boolean do_not_filter_override = false
        String? filter_set_name
        Boolean? vqslod_filter_by_site
        String fq_filter_set_info_table = "~{data_project}.~{default_dataset}.filter_set_info"
        String fq_filter_set_site_table = "~{data_project}.~{default_dataset}.filter_set_sites"
        String fq_filter_set_tranches_table = "~{data_project}.~{default_dataset}.filter_set_tranches"

        # if these are unset, default sensitivity levels will be used
        Float? snps_truth_sensitivity_filter_level_override
        Float? indels_truth_sensitivity_filter_level_override

        File? excluded_intervals
        Boolean? emit_pls = false

        String? service_account_json_path

        String output_file_base_name
        String? output_gcs_dir
        File? gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/ah_var_store_20210914/gatk-package-4.2.0.0-406-ga9206a2-SNAPSHOT-local.jar"
        Int local_disk_for_extract = 150
    }

    String fq_samples_to_extract_table = "~{fq_cohort_extract_table_prefix}__SAMPLES"
    String fq_cohort_extract_table  = "~{fq_cohort_extract_table_prefix}__DATA"

    call SplitIntervals {
      input:
          intervals = wgs_intervals,
          ref_fasta = reference,
          ref_fai = reference_index,
          ref_dict = reference_dict,
          scatter_count = scatter_count,
          output_gcs_dir = output_gcs_dir,
          service_account_json_path = service_account_json_path
    }

    call GetBQTableLastModifiedDatetime as fq_cohort_extract_table_datetime {
        input:
            query_project = query_project,
            fq_table = fq_cohort_extract_table,
            service_account_json_path = service_account_json_path
    }

    call GetBQTableLastModifiedDatetime as fq_samples_to_extract_table_datetime {
        input:
            query_project = query_project,
            fq_table = fq_samples_to_extract_table,
            service_account_json_path = service_account_json_path
    }

    scatter(i in range(scatter_count) ) {
        call ExtractTask {
            input:
                gatk_override                   = gatk_override,
                reference                       = reference,
                reference_index                 = reference_index,
                reference_dict                  = reference_dict,
                fq_samples_to_extract_table     = fq_samples_to_extract_table,
                interval_index                  = i,
                intervals                       = SplitIntervals.interval_files[i],
                fq_cohort_extract_table         = fq_cohort_extract_table,
                read_project_id                 = query_project,
                do_not_filter_override          = do_not_filter_override,
                fq_filter_set_info_table        = fq_filter_set_info_table,
                fq_filter_set_site_table        = fq_filter_set_site_table,
                fq_filter_set_tranches_table    = fq_filter_set_tranches_table,
                filter_set_name                 = filter_set_name,
                vqslod_filter_by_site           = vqslod_filter_by_site,
                snps_truth_sensitivity_filter_level = snps_truth_sensitivity_filter_level_override,
                indels_truth_sensitivity_filter_level = indels_truth_sensitivity_filter_level_override,
                excluded_intervals              = excluded_intervals,
                emit_pls                        = emit_pls,
                service_account_json_path       = service_account_json_path,
                output_file                     = "${output_file_base_name}_${i}.vcf.gz",
                output_gcs_dir                  = output_gcs_dir,
                local_disk                      = local_disk_for_extract,
                last_modified_timestamps        = [fq_samples_to_extract_table_datetime.last_modified_timestamp, fq_cohort_extract_table_datetime.last_modified_timestamp]
        }
    }

    call SumBytes {
      input:
        file_sizes_bytes = flatten([ExtractTask.output_vcf_bytes, ExtractTask.output_vcf_index_bytes])
    }

    call CreateManifest {
      input:
        manifest_lines = ExtractTask.manifest,
        output_gcs_dir = output_gcs_dir
    }

    output {
      Array[File] output_vcfs = ExtractTask.output_vcf
      Array[File] output_vcf_indexes = ExtractTask.output_vcf_index
      Float total_vcfs_size_mb = SumBytes.total_mb
      File manifest = CreateManifest.manifest
    }
}

################################################################################
task ExtractTask {
    input {
        # ------------------------------------------------
        # Input args:
        File reference
        File reference_index
        File reference_dict

        String fq_samples_to_extract_table

        Int interval_index
        File intervals

        String fq_cohort_extract_table
        String read_project_id
        String output_file
        String? output_gcs_dir

        Boolean do_not_filter_override
        String fq_filter_set_info_table
        String fq_filter_set_site_table
        String fq_filter_set_tranches_table
        String? filter_set_name
        Boolean? vqslod_filter_by_site
        Float? snps_truth_sensitivity_filter_level
        Float? indels_truth_sensitivity_filter_level

        File? excluded_intervals
        Boolean? emit_pls

        # Runtime Options:
        String? service_account_json_path
        File? gatk_override

        Int? local_sort_max_records_in_ram = 10000000
        Int local_disk

        # for call-caching -- check if DB tables haven't been updated since the last run
        Array[String] last_modified_timestamps
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        export GATK_LOCAL_JAR="~{default="/root/gatk.jar" gatk_override}"

        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
            gcloud config set project ~{read_project_id}
        fi

        df -h

        if [ ~{do_not_filter_override} = 'true' ]; then
            FILTERING_ARGS=''
        else
            FILTERING_ARGS='--filter-set-info-table ~{fq_filter_set_info_table}
                --filter-set-site-table ~{fq_filter_set_site_table}
                --tranches-table ~{fq_filter_set_tranches_table}
                --filter-set-name ~{filter_set_name}
                ~{true='--vqslod-filter-by-site' false='' vqslod_filter_by_site}
                ~{"--snps-truth-sensitivity-filter-level " + snps_truth_sensitivity_filter_level}
                ~{"--indels-truth-sensitivity-filter-level " + indels_truth_sensitivity_filter_level}'
        fi

        gatk --java-options "-Xmx9g" \
            ExtractCohort \
                --mode GENOMES --ref-version 38 \
                -R ~{reference} \
                -O ~{output_file} \
                --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
                --sample-table ~{fq_samples_to_extract_table} \
                --cohort-extract-table ~{fq_cohort_extract_table} \
                -L ~{intervals} \
                ~{"-XL " + excluded_intervals} \
                --project-id ~{read_project_id} \
                ~{true='--emit-pls' false='' emit_pls} \
                ${FILTERING_ARGS}

        OUTPUT_FILE_BYTES="$(du -b ~{output_file} | cut -f1)"
        echo ${OUTPUT_FILE_BYTES} > vcf_bytes.txt

        OUTPUT_FILE_INDEX_BYTES="$(du -b ~{output_file}.tbi | cut -f1)"
        echo ${OUTPUT_FILE_INDEX_BYTES} > vcf_index_bytes.txt

        if [ -n "${OUTPUT_GCS_DIR}" ]; then
          gsutil cp ~{output_file} ${OUTPUT_GCS_DIR}/
          gsutil cp ~{output_file}.tbi ${OUTPUT_GCS_DIR}/
          OUTPUT_FILE_DEST="${OUTPUT_GCS_DIR}/~{output_file}"
          OUTPUT_FILE_INDEX_DEST="${OUTPUT_GCS_DIR}/~{output_file}.tbi"
        else
          OUTPUT_FILE_DEST="~{output_file}"
          OUTPUT_FILE_INDEX_DEST="~{output_file}.tbi"
        fi

        # Parent Task will collect manifest lines and create a joined file
        # Currently, the schema is `[interval_number], [output_file_location], [output_file_size_bytes], [output_file_index_location], [output_file_size_bytes]`
        echo ~{interval_index},${OUTPUT_FILE_DEST},${OUTPUT_FILE_BYTES},${OUTPUT_FILE_INDEX_DEST},${OUTPUT_FILE_INDEX_BYTES} >> manifest.txt
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "12 GB"
        disks: "local-disk ~{local_disk} HDD"
        bootDiskSizeGb: 15
        preemptible: 2
        maxRetries: 3
        cpu: 2
    }

    # files sizes are floats instead of ints because they can be larger
    output {
        File output_vcf = "~{output_file}"
        Float output_vcf_bytes = read_float("vcf_bytes.txt")
        File output_vcf_index = "~{output_file}.tbi"
        Float output_vcf_index_bytes = read_float("vcf_index_bytes.txt")
        String manifest = read_string("manifest.txt")
    }
 }

 task SplitIntervals {
    input {
        File intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count
        String? split_intervals_extra_args
        String? output_gcs_dir

        File? gatk_override
        String? service_account_json_path
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    parameter_meta {
        intervals: {
            localization_optional: true
        }
        ref_fasta: {
            localization_optional: true
        }
        ref_fai: {
            localization_optional: true
        }
        ref_dict: {
            localization_optional: true
        }
     }

     command {
         set -e
         export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

         mkdir interval-files
         gatk --java-options "-Xmx2g" SplitIntervals \
             --dont-mix-contigs \
             -R ~{ref_fasta} \
             ~{"-L " + intervals} \
             -scatter ~{scatter_count} \
             -O interval-files \
             ~{split_intervals_extra_args}
         cp interval-files/*.interval_list .

         # Drop trailing slash if one exists
         OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

         if [ -n "$OUTPUT_GCS_DIR" ]; then
             if [ ~{has_service_account_file} = 'true' ]; then
                 gsutil cp ~{service_account_json_path} local.service_account.json
                 gcloud auth activate-service-account --key-file=local.service_account.json
             fi
             gsutil -m cp *.interval_list $OUTPUT_GCS_DIR/
         fi
     }

     runtime {
         docker: "us.gcr.io/broad-gatk/gatk:4.2.0.0"
         bootDiskSizeGb: 15
         memory: "3 GB"
         disks: "local-disk 10 HDD"
         preemptible: 3
         cpu: 1
     }

     output {
         Array[File] interval_files = glob("*.interval_list")
     }
 }

task GetBQTableLastModifiedDatetime {
    # because this is being used to determine if the data has changed, never use call cache
    meta {
        volatile: true
    }

    input {
        String query_project
        String fq_table
        String? service_account_json_path
    }

    String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

    # ------------------------------------------------
    # try to get the last modified date for the table in question; fail if something comes back from BigQuwey
    # that isn't in the right format (e.g. an error)
    command <<<
        set -e
        if [ ~{has_service_account_file} = 'true' ]; then
            gsutil cp ~{service_account_json_path} local.service_account.json
            gcloud auth activate-service-account --key-file=local.service_account.json
        fi

        echo "project_id = ~{query_project}" > ~/.bigqueryrc

        # bq needs the project name to be separate by a colon
        DATASET_TABLE_COLON=$(echo ~{fq_table} | sed 's/\./:/')

        LASTMODIFIED=$(bq --location=US --project_id=~{query_project} --format=json show ${DATASET_TABLE_COLON} | python3 -c "import sys, json; print(json.load(sys.stdin)['lastModifiedTime']);")
        if [[ $LASTMODIFIED =~ ^[0-9]+$ ]]; then
            echo $LASTMODIFIED
        else
            exit 1
        fi
    >>>

    output {
        String last_modified_timestamp = read_string(stdout())
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        cpu: 1
    }
}

task SumBytes {

  input {
    Array[Float] file_sizes_bytes
  }

  command <<<
    set -e
    echo "~{sep=" " file_sizes_bytes}" | tr " " "\n" | python -c "
    import sys;
    total_bytes = sum(float(i.strip()) for i in sys.stdin);
    total_mb = total_bytes/10**6;
    print(total_mb);"
  >>>

  output {
    Float total_mb = read_float(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task CreateManifest {

    input {
        Array[String] manifest_lines
        String? output_gcs_dir
    }

    command <<<
        set -e
        MANIFEST_LINES_TXT=~{write_lines(manifest_lines)}
        echo "vcf_file_location, vcf_file_bytes, vcf_index_location, vcf_index_bytes" >> manifest.txt
        sort -n ${MANIFEST_LINES_TXT} | cut -d',' -f 2- >> manifest.txt

        # Drop trailing slash if one exists
        OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

        if [ -n "${OUTPUT_GCS_DIR}" ]; then
          gsutil cp manifest.txt ${OUTPUT_GCS_DIR}/
        fi
    >>>

    output {
        File manifest = "manifest.txt"
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        cpu: 1
    }
}
