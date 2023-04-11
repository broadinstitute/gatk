version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsExtractCallset {
    input {
        Boolean go = true
        # The name of the bigquery dataset containing the GVS data we are extracting
        String dataset_name
        # The project id for the specified bigquery dataset
        String project_id
        # The name of the callset we created using GvsPrepareRangesCallset to extract
        String call_set_identifier

        # Chromosome code corresponding to the set of contig names used by plink that matches those used by the
        # reference (see plink chromosome codes here: https://www.cog-genomics.org/plink/2.0/data#irreg_output)
        # ExtractCohortToPgen currently only supports codes "chrM" and "MT"
        String pgen_chromosome_code = "chrM"
        # Max number of alt alleles a site can have. If a site exceeds this number, it will not be written (max 254)
        Int max_alt_alleles = 254
        # If true, does not throw an exception for samples@sites with unsupported ploidy (codes it as missing instead)
        Boolean lenient_ploidy_validation = false

        # The project id for the bigquery dataset containing the cohort we created using GvsPrepareRangesCallset
        String cohort_project_id = project_id
        # The name of the bigquery dataset containing the cohort we created using GvsPrepareRangesCallset
        String cohort_dataset_name = dataset_name
        Boolean do_not_filter_override = false
        Boolean control_samples = false
        String extract_table_prefix
        String filter_set_name
        String query_project = project_id
        # This is optional now since the workflow will choose an appropriate value below if this is unspecified.
        Int? scatter_count
        Int? memory_override
        Int? disk_override
        # If true, the output files will be named with a zero-padded interval number prefix
        # If false, the output files will be named with a non-zero-padded interval number suffix
        Boolean zero_pad_output_pgen_filenames = true

        # set to "NONE" if all the reference data was loaded into GVS in GvsImportGenomes
        String drop_state = "NONE"

        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
        Boolean use_interval_weights = true
        File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"

        File? gatk_override
        String? extract_docker_override
        String gatk_docker

        String output_file_base_name = filter_set_name

        Int? extract_maxretries_override
        Int? extract_preemptible_override
        String? output_gcs_dir
        String? split_intervals_extra_args
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        Float x_bed_weight_scaling = 4
        Float y_bed_weight_scaling = 4
        Boolean write_cost_to_db = true
    }

    File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
    File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
    File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

    String fq_gvs_dataset = "~{project_id}.~{dataset_name}"
    String fq_cohort_dataset = "~{cohort_project_id}.~{cohort_dataset_name}"

    String full_extract_prefix = if (control_samples) then "~{extract_table_prefix}_controls" else extract_table_prefix

    String fq_filter_set_info_table = "~{fq_gvs_dataset}.filter_set_info"
    String fq_filter_set_site_table = "~{fq_gvs_dataset}.filter_set_sites"
    String fq_filter_set_tranches_table = "~{fq_gvs_dataset}.filter_set_tranches"
    String fq_sample_table = "~{fq_gvs_dataset}.sample_info"
    String fq_cohort_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__DATA"
    String fq_ranges_cohort_ref_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__REF_DATA"
    String fq_ranges_cohort_vet_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__VET_DATA"

    String fq_samples_to_extract_table = "~{fq_cohort_dataset}.~{full_extract_prefix}__SAMPLES"
    Array[String] tables_patterns_for_datetime_check = ["~{full_extract_prefix}__%"]

    Boolean emit_pls = false
    Boolean emit_ads = true

    String intervals_file_extension = if (zero_pad_output_pgen_filenames) then '-~{output_file_base_name}.interval_list' else '-scattered.interval_list'

    call Utils.ScaleXYBedValues {
        input:
            interval_weights_bed = interval_weights_bed,
            x_bed_weight_scaling = x_bed_weight_scaling,
            y_bed_weight_scaling = y_bed_weight_scaling
    }

    call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
        input:
            project_id = project_id,
            fq_table = fq_sample_table
    }

    call Utils.GetNumSamplesLoaded {
        input:
            fq_sample_table = fq_samples_to_extract_table,
            project_id = project_id,
            control_samples = control_samples,
            sample_table_timestamp = SamplesTableDatetimeCheck.last_modified_timestamp
    }

    Int effective_scatter_count = if defined(scatter_count) then select_first([scatter_count])
                                  else if GetNumSamplesLoaded.num_samples < 100 then 100 # Quickstart
                                       else if GetNumSamplesLoaded.num_samples < 1000 then 500
                                            else if GetNumSamplesLoaded.num_samples < 5000 then 1000
                                                 else if GetNumSamplesLoaded.num_samples < 20000 then 2000 # Stroke Anderson
                                                      else if GetNumSamplesLoaded.num_samples < 50000 then 10000
                                                           else if GetNumSamplesLoaded.num_samples < 100000 then 20000 # Charlie
                                                                else 34000

    # WDL 1.0 trick to set a variable ('none') to be undefined.
    if (false) {
        File? none = ""
    }

    call Utils.SplitIntervalsTarred {
        input:
            intervals = interval_list,
            ref_fasta = reference,
            ref_fai = reference_index,
            ref_dict = reference_dict,
            interval_weights_bed = if (use_interval_weights) then ScaleXYBedValues.xy_scaled_bed else none,
            intervals_file_extension = intervals_file_extension,
            scatter_count = effective_scatter_count,
            output_gcs_dir = output_gcs_dir,
            split_intervals_extra_args = split_intervals_extra_args,
            split_intervals_disk_size_override = split_intervals_disk_size_override,
            split_intervals_mem_override = split_intervals_mem_override,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    call Utils.GetBQTableLastModifiedDatetime as FilterSetInfoTimestamp {
        input:
            project_id = project_id,
            fq_table = "~{fq_filter_set_info_table}"
    }

    if ( !do_not_filter_override ) {
        call Utils.ValidateFilterSetName {
            input:
                project_id = query_project,
                fq_filter_set_info_table = "~{fq_filter_set_info_table}",
                filter_set_name = filter_set_name,
                filter_set_info_timestamp = FilterSetInfoTimestamp.last_modified_timestamp
        }

        call Utils.IsVQSRLite {
            input:
                project_id = query_project,
                fq_filter_set_info_table = "~{fq_filter_set_info_table}",
                filter_set_name = filter_set_name
        }
    }

    # If we're not using the VQSR filters, set it to Lite (really shouldn't matter one way or the other)
    # Otherwise use the auto-derived flag.
    Boolean use_VQSR_lite = select_first([IsVQSRLite.is_vqsr_lite, true])

    call Utils.GetBQTablesMaxLastModifiedTimestamp {
        input:
            query_project = query_project,
            data_project = project_id,
            dataset_name = dataset_name,
            table_patterns = tables_patterns_for_datetime_check
    }

    scatter(i in range(length(SplitIntervalsTarred.interval_filenames))) {
        String interval_filename = SplitIntervalsTarred.interval_filenames[i]
        String pgen_basename = if (zero_pad_output_pgen_filenames) then sub(interval_filename, ".interval_list", "") else "~{output_file_base_name}_${i}"
        call ExtractTask {
            input:
                go                                 = select_first([ValidateFilterSetName.done, true]),
                dataset_name                       = dataset_name,
                call_set_identifier                = call_set_identifier,
                pgen_chromosome_code               = pgen_chromosome_code,
                max_alt_alleles                    = max_alt_alleles,
                lenient_ploidy_validation          = lenient_ploidy_validation,
                use_VQSR_lite                      = use_VQSR_lite,
                gatk_override                      = gatk_override,
                reference                          = reference,
                reference_index                    = reference_index,
                reference_dict                     = reference_dict,
                fq_samples_to_extract_table        = fq_samples_to_extract_table,
                interval_index                     = i,
                interval_files_tar                 = SplitIntervalsTarred.interval_files_tar,
                interval_filename                  = interval_filename,
                fq_cohort_extract_table            = fq_cohort_extract_table,
                fq_ranges_cohort_ref_extract_table = fq_ranges_cohort_ref_extract_table,
                fq_ranges_cohort_vet_extract_table = fq_ranges_cohort_vet_extract_table,
                read_project_id                    = query_project,
                do_not_filter_override             = do_not_filter_override,
                fq_filter_set_info_table           = fq_filter_set_info_table,
                fq_filter_set_site_table           = fq_filter_set_site_table,
                fq_filter_set_tranches_table       = if (use_VQSR_lite) then none else fq_filter_set_tranches_table,
                filter_set_name                    = filter_set_name,
                drop_state                         = drop_state,
                output_pgen_basename               = pgen_basename,
                zero_pad_output_pgen_filenames     = zero_pad_output_pgen_filenames,
                output_gcs_dir                     = output_gcs_dir,
                max_last_modified_timestamp        = GetBQTablesMaxLastModifiedTimestamp.max_last_modified_timestamp,
                extract_preemptible_override       = extract_preemptible_override,
                extract_maxretries_override        = extract_maxretries_override,
                disk_override                      = disk_override,
                memory_override                    = memory_override,
                emit_pls                           = emit_pls,
                emit_ads                           = emit_ads,
                write_cost_to_db                   = write_cost_to_db,
                docker_override                    = select_first([extract_docker_override, gatk_docker])
        }
    }

    call SumBytes {
        input:
            file_sizes_bytes = flatten([ExtractTask.output_pgen_bytes, ExtractTask.output_pvar_bytes, ExtractTask.output_psam_bytes])
    }

    call CreateManifest {
        input:
            manifest_lines = ExtractTask.manifest,
            output_gcs_dir = output_gcs_dir
    }

    if (control_samples == false) {

        call Utils.GetBQTableLastModifiedDatetime {
            input:
                project_id = query_project,
                fq_table = fq_samples_to_extract_table
        }

        call GenerateSampleListFile {
            input:
                fq_samples_to_extract_table = fq_samples_to_extract_table,
                samples_to_extract_table_timestamp = GetBQTableLastModifiedDatetime.last_modified_timestamp,
                output_gcs_dir = output_gcs_dir,
                query_project = query_project
        }
    }

    output {
        Array[File] output_pgens = ExtractTask.output_pgen
        Array[File] output_pvars = ExtractTask.output_pvar
        Array[File] output_psams = ExtractTask.output_psam
        File output_pgen_interval_files = SplitIntervalsTarred.interval_files_tar
        Array[String] output_pgen_interval_filenames = SplitIntervalsTarred.interval_filenames
        Float total_pgens_size_mb = SumBytes.total_mb
        File manifest = CreateManifest.manifest
        File? sample_name_list = GenerateSampleListFile.sample_name_list
        Boolean done = true
    }
}

task ExtractTask {
    input {
        Boolean go

        String dataset_name
        String call_set_identifier

        # Chromosome code corresponding to the set of contig names used by plink that matches those used by the
        # reference (see plink chromosome codes here: https://www.cog-genomics.org/plink/2.0/data#irreg_output)
        # ExtractCohortToPgen currently only supports codes "chrM" and "MT"
        String pgen_chromosome_code
        # Max number of alt alleles a site can have. If a site exceeds this number, it will not be written (max 254)
        Int? max_alt_alleles
        # If true, does not throw an exception for samples@sites with unsupported ploidy (codes it as missing instead)
        Boolean? lenient_ploidy_validation

        Boolean use_VQSR_lite

        File reference
        File reference_index
        File reference_dict

        String fq_samples_to_extract_table

        Int interval_index
        File interval_files_tar
        String interval_filename
        String drop_state

        String fq_cohort_extract_table
        String fq_ranges_cohort_ref_extract_table
        String fq_ranges_cohort_vet_extract_table
        String read_project_id
        String output_pgen_basename
        Boolean zero_pad_output_pgen_filenames
        String? output_gcs_dir

        String cost_observability_tablename = "cost_observability"

        Boolean emit_pls
        Boolean emit_ads

        Boolean do_not_filter_override
        String fq_filter_set_info_table
        String fq_filter_set_site_table
        String? fq_filter_set_tranches_table
        String? filter_set_name
        Boolean write_cost_to_db

        # Runtime Options:
        File? gatk_override
        String? docker_override
        Int? extract_preemptible_override
        Int? extract_maxretries_override
        Int? disk_override
        Int? memory_override

        Int? local_sort_max_records_in_ram = 10000000

        # for call-caching -- check if DB tables haven't been updated since the last run
        String max_last_modified_timestamp
    }
    meta {
        # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    String cost_observability_line = if (write_cost_to_db == true) then "--cost-observability-tablename ~{cost_observability_tablename}" else ""

    String inferred_reference_state = if (drop_state == "NONE") then "ZERO" else drop_state

    command <<<
        set -eux

        bash ~{monitoring_script} > monitoring.log &

        export GATK_LOCAL_JAR="~{default="/root/gatk.jar" gatk_override}"

        if [ ~{do_not_filter_override} = true ]; then
            FILTERING_ARGS=''
        elif [ ~{use_VQSR_lite} = false ]; then
            FILTERING_ARGS='--filter-set-info-table ~{fq_filter_set_info_table}
            --filter-set-site-table ~{fq_filter_set_site_table}
            --tranches-table ~{fq_filter_set_tranches_table}
            --filter-set-name ~{filter_set_name}'
        else
            FILTERING_ARGS='--filter-set-info-table ~{fq_filter_set_info_table}
            --filter-set-site-table ~{fq_filter_set_site_table}
            --filter-set-name ~{filter_set_name}'
        fi

        touch writer.log

        # Extract the intervals file from the intervals tarball
        tar -xf ~{interval_files_tar} ~{interval_filename}
        
        # Calculate the memory size we'll use for extraction as 3/4 of the total memory
        JAVA_MEM=$(echo "scale=0; 0.75 * ${MEM_SIZE}" | bc)
        JAVA_MEM=$(printf "%.0f" ${JAVA_MEM})

        gatk --java-options "-Xmx${JAVA_MEM}g" \
        ExtractCohortToPgen \
        --vet-ranges-extract-fq-table ~{fq_ranges_cohort_vet_extract_table} \
        --ref-ranges-extract-fq-table ~{fq_ranges_cohort_ref_extract_table} \
        --ref-version 38 \
        -R ~{reference} \
        -O ~{output_pgen_basename}.pgen \
        --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
        --sample-table ~{fq_samples_to_extract_table} \
        ~{"--inferred-reference-state " + inferred_reference_state} \
        -L ~{interval_filename} \
        --project-id ~{read_project_id} \
        ~{true='--emit-pls' false='' emit_pls} \
        ~{true='--emit-ads' false='' emit_ads} \
        ~{true='' false='--use-vqsr-classic-scoring' use_VQSR_lite} \
        ${FILTERING_ARGS} \
        --dataset-id ~{dataset_name} \
        --call-set-identifier ~{call_set_identifier} \
        --wdl-step GvsExtractCallset \
        --wdl-call ExtractTask \
        --shard-identifier ~{interval_filename} \
        ~{cost_observability_line} \
        --writer-log-file writer.log \
        --pgen-chromosome-code ~{pgen_chromosome_code} \
        --max-alt-alleles ~{max_alt_alleles} \
        ~{true='--lenient-ploidy-validation' false='' lenient_ploidy_validation} \
        --allow-empty-pgen


        # Drop trailing slash if one exists
        OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

        OUTPUT_FILE_BYTES="$(du -b ~{output_pgen_basename}.pgen | cut -f1)"
        echo ${OUTPUT_FILE_BYTES} > pgen_bytes.txt

        OUTPUT_FILE_PVAR_BYTES="$(du -b ~{output_pgen_basename}.pvar.zst | cut -f1)"
        echo ${OUTPUT_FILE_PVAR_BYTES} > pvar_bytes.txt

        OUTPUT_FILE_PSAM_BYTES="$(du -b ~{output_pgen_basename}.psam | cut -f1)"
        echo ${OUTPUT_FILE_PSAM_BYTES} > psam_bytes.txt

        if [ -n "${OUTPUT_GCS_DIR}" ]; then
            gsutil cp ~{output_pgen_basename}.pgen ${OUTPUT_GCS_DIR}/
            gsutil cp ~{output_pgen_basename}.pvar.zst ${OUTPUT_GCS_DIR}/
            gsutil cp ~{output_pgen_basename}.psam ${OUTPUT_GCS_DIR}/
            OUTPUT_FILE_DEST="${OUTPUT_GCS_DIR}/~{output_pgen_basename}.pgen"
            OUTPUT_FILE_PVAR_DEST="${OUTPUT_GCS_DIR}/~{output_pgen_basename}.pvar.zst"
            OUTPUT_FILE_PSAM_DEST="${OUTPUT_GCS_DIR}/~{output_pgen_basename}.psam"
        else
            OUTPUT_FILE_DEST="~{output_pgen_basename}.pgen"
            OUTPUT_FILE_PVAR_DEST="~{output_pgen_basename}.pvar.zst"
            OUTPUT_FILE_PSAM_DEST="~{output_pgen_basename}.psam"
        fi

        # Parent Task will collect manifest lines and create a joined file
        # Currently, the schema is `[interval_number], [output_file_location], [output_file_size_bytes], [output_file_pvar_location], [output_file_pvar_size_bytes], [output_file_psam_location], [output_file_psam_size_bytes]`
        echo ~{interval_index},${OUTPUT_FILE_DEST},${OUTPUT_FILE_BYTES},${OUTPUT_FILE_PVAR_DEST},${OUTPUT_FILE_PVAR_BYTES},${OUTPUT_FILE_PSAM_DEST},${OUTPUT_FILE_PSAM_BYTES} >> manifest.txt
    >>>
    runtime {
        docker: select_first([docker_override, "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds:klydon-fc3a5d17fbab35f51e58e740452ead60b908d6b7-4.4.0.0-78-gfc3a5d17f"])
        memory: select_first([memory_override, 12]) + " GB"
        disks: "local-disk " + select_first([disk_override, 150]) + " HDD"
        bootDiskSizeGb: 15
        preemptible: select_first([extract_preemptible_override, "2"])
        maxRetries: select_first([extract_maxretries_override, "3"])
        cpu: 2
    }

    # files sizes are floats instead of ints because they can be larger
    output {
        File output_pgen = "~{output_pgen_basename}.pgen"
        Float output_pgen_bytes = read_float("pgen_bytes.txt")
        File output_pvar = "~{output_pgen_basename}.pvar.zst"
        Float output_pvar_bytes = read_float("pvar_bytes.txt")
        File output_psam = "~{output_pgen_basename}.psam"
        Float output_psam_bytes = read_float("psam_bytes.txt")
        String manifest = read_string("manifest.txt")
        File writer_log = "writer.log"
        File monitoring_log = "monitoring.log"
    }
}

task SumBytes {
    input {
        Array[Float] file_sizes_bytes
    }
    meta {
        # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
    }

    command <<<
        set -e
        echo "~{sep=" " file_sizes_bytes}" | tr " " "\n" | python3 -c "
        import sys;
        total_bytes = sum(float(i.strip()) for i in sys.stdin);
        total_mb = total_bytes/10**6;
        print(total_mb);"
    >>>
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
        memory: "3 GB"
        disks: "local-disk 500 HDD"
        preemptible: 3
        cpu: 1
    }

    output {
        Float total_mb = read_float(stdout())
    }
}

task CreateManifest {
    input {
        Array[String] manifest_lines
        String? output_gcs_dir
    }
    meta {
        # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
    }

    command <<<
        set -e
        MANIFEST_LINES_TXT=~{write_lines(manifest_lines)}
        echo "pgen_file_location, pgen_file_bytes, pvar_file_location, pvar_file_bytes, psam_file_location, psam_file_bytes" >> manifest.txt
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
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
        memory: "3 GB"
        disks: "local-disk 500 HDD"
        preemptible: 3
        cpu: 1
    }
}

task GenerateSampleListFile {
    input {
        String fq_samples_to_extract_table
        String samples_to_extract_table_timestamp
        String query_project

        String? output_gcs_dir
    }
    meta {
        # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
    }

    # add labels for DSP Cloud Cost Control Labeling and Reporting
    String bq_labels = "--label service:gvs --label team:variants --label managedby:extract_callset"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Drop trailing slash if one exists
        OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

        echo "project_id = ~{query_project}" > ~/.bigqueryrc

        bq --apilog=false --project_id=~{query_project} --format=csv query --use_legacy_sql=false ~{bq_labels} 'SELECT sample_name FROM `~{fq_samples_to_extract_table}`' | sed 1d > sample-name-list.txt

        if [ -n "$OUTPUT_GCS_DIR" ]; then
        gsutil cp sample-name-list.txt ${OUTPUT_GCS_DIR}/
        fi
    >>>
    output {
        File sample_name_list = "sample-name-list.txt"
    }

    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:426.0.0-alpine"
        memory: "3 GB"
        disks: "local-disk 500 HDD"
        preemptible: 3
        cpu: 1
    }
}