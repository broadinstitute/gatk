version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsBenchmarkExtractTask {
    input {
        String? git_branch_or_tag
        String data_project
        String dataset_name

        File wgs_intervals

        File reference
        File reference_index
        File reference_dict

        # For reblocking v1, the default is "SIXTY" instead of "FORTY"
        String? drop_state = "FORTY"

        # NOTE: this is just the cohort table prefix, not including project or dataset qualifiers
        # without a default value, ranges users are forced to specify a value even though it is meaningless
        String extract_table_prefix = ""
        String query_project = data_project
        String fq_ranges_dataset = "~{data_project}.~{dataset_name}"

        Boolean do_not_filter_override = false
        String? filter_set_name
        Boolean? vqslod_filter_by_site
        String fq_filter_set_info_table = "~{data_project}.~{dataset_name}.filter_set_info"
        String fq_filter_set_site_table = "~{data_project}.~{dataset_name}.filter_set_sites"
        String fq_filter_set_tranches_table = "~{data_project}.~{dataset_name}.filter_set_tranches"

        # if these are unset, default sensitivity levels will be used
        Float? snps_truth_sensitivity_filter_level_override
        Float? indels_truth_sensitivity_filter_level_override

        File? excluded_intervals
        Boolean? emit_pls = false

        Int? extract_cpu_override = 2
        String? extract_memory_override = "12 GB"

        Int? extract_preemptible_override
        Int? extract_maxretries_override
        Int? split_intervals_disk_size_override

        String mode = "RANGES"

        String output_file_base_name
        String? output_gcs_dir
        File? gatk_override
        Int local_disk_for_extract = 150

        String fq_samples_to_extract_table = "~{data_project}.~{dataset_name}.~{extract_table_prefix}__SAMPLES"
        String fq_cohort_extract_table  = "~{data_project}.~{dataset_name}.~{extract_table_prefix}__DATA"

        String? gatk_docker
    }

    # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
    # no calling WDLs that might supply `git_hash`).
    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_cloud_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])

    call ExtractTask {
        input:
            gatk_override                   = gatk_override,
            reference                       = reference,
            reference_index                 = reference_index,
            reference_dict                  = reference_dict,
            fq_samples_to_extract_table     = fq_samples_to_extract_table,
            intervals                       = wgs_intervals,
            fq_cohort_extract_table         = fq_cohort_extract_table,
            read_project_id                 = query_project,
            do_not_filter_override          = do_not_filter_override,
            fq_ranges_dataset               = fq_ranges_dataset,
            fq_filter_set_info_table        = fq_filter_set_info_table,
            fq_filter_set_site_table        = fq_filter_set_site_table,
            fq_filter_set_tranches_table    = fq_filter_set_tranches_table,
            filter_set_name                 = filter_set_name,
            vqslod_filter_by_site           = vqslod_filter_by_site,
            snps_truth_sensitivity_filter_level = snps_truth_sensitivity_filter_level_override,
            indels_truth_sensitivity_filter_level = indels_truth_sensitivity_filter_level_override,
            excluded_intervals              = excluded_intervals,
            emit_pls                        = emit_pls,
            drop_state                      = drop_state,
            output_file                     = "${output_file_base_name}.vcf.gz",
            local_disk                      = local_disk_for_extract,
            extract_preemptible_override    = extract_preemptible_override,
            extract_cpu_override            = extract_cpu_override,
            extract_memory_override         = extract_memory_override,
            extract_maxretries_override     = extract_maxretries_override,
            gatk_docker                     = effective_cloud_gatk_docker,
    }

    output {
        String recorded_git_hash = GetToolVersions.git_hash
    }
}

################################################################################
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

        String fq_samples_to_extract_table

        File intervals
        String? drop_state

        String fq_cohort_extract_table
        String read_project_id
        String output_file

        Boolean do_not_filter_override
        String fq_ranges_dataset
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
        File? gatk_override
        String gatk_docker
        Int? extract_preemptible_override
        Int? extract_maxretries_override
        String? extract_memory_override
        Int? extract_cpu_override

        Int? local_sort_max_records_in_ram = 10000000
        Int local_disk

    }

    # ------------------------------------------------
    # Run our command:
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        export GATK_LOCAL_JAR="~{default="/root/gatk.jar" gatk_override}"

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
            ExtractCohortToVcf \
                --mode RANGES --vet-ranges-fq-dataset ~{fq_ranges_dataset} \
                --ref-version 38 \
                -R ~{reference} \
                -O ~{output_file} \
                --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
                --sample-table ~{fq_samples_to_extract_table} \
                ~{"--inferred-reference-state " + drop_state} \
                -L ~{intervals} \
                ~{"-XL " + excluded_intervals} \
                --project-id ~{read_project_id} \
                ~{true='--emit-pls' false='' emit_pls} \
                ${FILTERING_ARGS}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: gatk_docker
        memory: select_first([extract_memory_override, "12 GB"])
        disks: "local-disk ~{local_disk} HDD"
        bootDiskSizeGb: 15
        preemptible: select_first([extract_preemptible_override, "2"])
        maxRetries: select_first([extract_maxretries_override, "3"])
        cpu: select_first([extract_cpu_override, "2"])
    }

    output {
        File output_vcf = "~{output_file}"
        File output_vcf_index = "~{output_file}.tbi"
    }
}

task SumBytes {

  input {
    Array[Float] file_sizes_bytes
    String cloud_sdk_docker
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

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
    docker: cloud_sdk_docker
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
        String cloud_sdk_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        MANIFEST_LINES_TXT=~{write_lines(manifest_lines)}
        echo "vcf_file_location, vcf_file_bytes, vcf_index_location, vcf_index_bytes" >> manifest.txt
        sort -n ${MANIFEST_LINES_TXT} | cut -d',' -f 2- >> manifest.txt

        # Drop trailing slash if one exists
        OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

        if [ -n "$OUTPUT_GCS_DIR" ]; then
          gsutil cp manifest.txt ${OUTPUT_GCS_DIR}/
        fi
    >>>

    output {
        File manifest = "manifest.txt"
    }

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        cpu: 1
    }
}
