version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsImportGenomes {

  input {
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go = true
    String? git_branch_or_tag
    String? git_hash
    String dataset_name
    String project_id

    Int num_samples
    File external_sample_names
    File input_vcfs
    File input_vcf_indexes

    Boolean skip_loading_vqsr_fields = false
    Boolean use_compressed_references = false
    # Turn Parquet lifecycle configuration off by default as pet service accounts don't seem to automatically get the
    # required permissions on the workspace bucket for this to work.
    Boolean configure_parquet_lifecycle = false

    # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
    String drop_state = "NONE"
    # beta customers will almost always have a naive GCP account, and as such will not be able to cross over their quotas
    # without Google shutting import down by throwing them API errors.  For them, we limit our scattering.
    Boolean is_rate_limited_beta_customer = false
    # This was determined to be the point at which we come close to but don't cross over the "AppendRows throughput per
    # project for small regions per minute per region" default quota of ~19G.  Uses up to ~90% of the quota at peaks
    # without going over
    Int beta_customer_max_scatter = 200

    String reference_name = "hg38"
    File? interval_list

    Int? load_data_scatter_width
    Int? load_data_preemptible_override
    Int? load_data_maxretries_override
    # At least one of these "load" inputs must be true
    Boolean load_vet_and_ref_ranges = true
    Boolean load_vcf_headers = false
    String? basic_docker
    String? cloud_sdk_docker
    String? variants_docker
    String? gatk_docker
    File? load_data_gatk_override
    String? billing_project_id

    Boolean use_parquet_ingest = true
    # Dump these Parquet files to a bucket.
    String? parquet_output_gcs_dir

    # Delete parquet files from GCS after successfully loading them into BigQuery
    Boolean delete_parquet_files_after_loading = true
    Boolean use_alternate_parquet_delete_strategy = false

    Boolean is_wgs = true
  }

  Int max_auto_scatter_width = if is_wgs then 25000 else 100000
  String genome_type = if is_wgs then "WGS" else "exome"

  # Broad users enjoy higher quotas and can scatter more widely than beta users before BigQuery smacks them
  # We don't expect this to be changed at runtime, so we can keep this as a constant defined in here
  Int broad_user_max_scatter = 1000

  # figure out max scatter depending on whether they're a Broad internal user or a beta customer.
  Int max_scatter_for_user =  if is_rate_limited_beta_customer then beta_customer_max_scatter
                              else broad_user_max_scatter

  if (!defined(git_hash) || !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(variants_docker) || !defined(gatk_docker)) {
    call Utils.GetToolVersions {
      input:
        git_branch_or_tag = git_branch_or_tag,
    }
  }

  String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])

  if (use_parquet_ingest && !defined(parquet_output_gcs_dir)) {
    call Utils.TerminateWorkflow as MustDefineOutputDirForParquetIngest {
      input:
        message = "use_parquet_ingest set to true but parquet_output_gcs_dir not defined",
        basic_docker = effective_basic_docker,
    }
  }

  if (load_vcf_headers && use_parquet_ingest) {
    call Utils.TerminateWorkflow as CannotLoadHeadersWithParquetIngest {
      input:
      message = "The combination of Parquet ingest and VCF header loading is not currently supported.",
      basic_docker = effective_basic_docker,
    }
  }

  call Utils.GetReference {
    input:
      reference_name = reference_name,
      basic_docker = effective_basic_docker,
  }

  File effective_interval_list = select_first([interval_list, GetReference.reference.wgs_calling_interval_list])

  if (!load_vcf_headers && !load_vet_and_ref_ranges) {
    call Utils.TerminateWorkflow as MustLoadAtLeastOneThing {
      input:
        message = "GvsImportGenomes called with both load_vcf_headers and load_vet_and_ref_ranges set to false",
        basic_docker = effective_basic_docker,
    }
  }

  if ((num_samples > max_auto_scatter_width) && !(defined(load_data_scatter_width))) {
    call Utils.TerminateWorkflow as DieDueToTooManySamplesWithoutExplicitLoadDataScatterWidth {
      input:
        message = "Importing " + num_samples + " samples but 'load_data_scatter_width' is not explicitly specified; the limit for automatic scatter width selection is " + max_auto_scatter_width + " for " + genome_type + " samples.",
        basic_docker = effective_basic_docker,
    }
  }

  # At least 1, per limits above not more than 20.
  # But if it's a beta customer, use the number computed above
  Int effective_load_data_batch_size = if (defined(load_data_scatter_width)) then num_samples / select_first([load_data_scatter_width])
                                       else if num_samples < max_scatter_for_user then 1
                                         else if is_wgs then num_samples / max_scatter_for_user
                                           else if num_samples < 5001 then 20
                                             else if num_samples < 20001 then 100
                                               else if num_samples < 50001 then 500
                                                 else 1000

  # Both preemptible and maxretries should be scaled up alongside import batch size since the likelihood of preemptions
  # and retryable random BQ import errors increases with import batch size / job run time.

  # At least 3, per limits above not more than 5.
  Int effective_load_data_preemptible = if (defined(load_data_preemptible_override)) then select_first([load_data_preemptible_override])
                                        else if effective_load_data_batch_size < 12 then 3
                                          else effective_load_data_batch_size / 4

  Int effective_load_data_maxretries = select_first([load_data_maxretries_override, 5])

  call CreateSampleDataViews {
    input:
      project_id = project_id,
      dataset_name = dataset_name,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call GetUningestedSampleIds {
    input:
      go = CreateSampleDataViews.done,
      dataset_name = dataset_name,
      project_id = project_id,
      external_sample_names = external_sample_names,
      num_samples = num_samples,
      table_name = "sample_info",
      load_vet_and_ref_ranges = load_vet_and_ref_ranges,
      load_vcf_headers = load_vcf_headers,
      cloud_sdk_docker = effective_cloud_sdk_docker,
  }

  call CurateInputLists {
    input:
      input_vcf_index_list = input_vcf_indexes,
      input_vcf_list = input_vcfs,
      input_sample_name_list = external_sample_names,
      input_samples_to_be_loaded_map = GetUningestedSampleIds.sample_map,
      variants_docker = effective_variants_docker,
  }

  call CreateFOFNs {
    input:
      batch_size = effective_load_data_batch_size,
      input_vcf_index_list = CurateInputLists.input_vcf_indexes,
      input_vcf_list = CurateInputLists.input_vcfs,
      sample_name_list = CurateInputLists.sample_name_list,
      basic_docker = effective_basic_docker,
  }

  scatter (i in range(length(CreateFOFNs.vcf_sample_name_fofns))) {
    if (use_parquet_ingest) {
      call ProcessInputGVCFs as GenerateParquetFilesFromInputGVCFs {
        input:
          index = i,
          dataset_name = dataset_name,
          project_id = project_id,
          skip_loading_vqsr_fields = skip_loading_vqsr_fields,
          drop_state = drop_state,
          drop_state_includes_greater_than = false,
          input_vcf_indexes = read_lines(CreateFOFNs.vcf_batch_vcf_index_fofns[i]),
          input_vcfs = read_lines(CreateFOFNs.vcf_batch_vcf_fofns[i]),
          interval_list = effective_interval_list,
          gatk_docker = effective_gatk_docker,
          gatk_override = load_data_gatk_override,
          load_data_preemptible = effective_load_data_preemptible,
          load_data_maxretries = effective_load_data_maxretries,
          sample_names = read_lines(CreateFOFNs.vcf_sample_name_fofns[i]),
          sample_map = GetUningestedSampleIds.sample_map,
          load_vet_and_ref_ranges = load_vet_and_ref_ranges,
          load_vcf_headers = load_vcf_headers,
          billing_project_id = billing_project_id,
          use_compressed_references = use_compressed_references,
          parquet_output_gcs_dir = parquet_output_gcs_dir,
          use_parquet_ingest = true,
      }
    }
    if (!use_parquet_ingest) { # WDL 1.1 does not have an else statement
      call ProcessInputGVCFs as LoadDataViaBigQueryWriteAPI {
        input:
          index = i,
          dataset_name = dataset_name,
          project_id = project_id,
          skip_loading_vqsr_fields = skip_loading_vqsr_fields,
          drop_state = drop_state,
          drop_state_includes_greater_than = false,
          input_vcf_indexes = read_lines(CreateFOFNs.vcf_batch_vcf_index_fofns[i]),
          input_vcfs = read_lines(CreateFOFNs.vcf_batch_vcf_fofns[i]),
          interval_list = effective_interval_list,
          gatk_docker = effective_gatk_docker,
          gatk_override = load_data_gatk_override,
          load_data_preemptible = effective_load_data_preemptible,
          load_data_maxretries = effective_load_data_maxretries,
          sample_names = read_lines(CreateFOFNs.vcf_sample_name_fofns[i]),
          sample_map = GetUningestedSampleIds.sample_map,
          load_vet_and_ref_ranges = load_vet_and_ref_ranges,
          load_vcf_headers = load_vcf_headers,
          billing_project_id = billing_project_id,
          use_compressed_references = use_compressed_references,
          use_parquet_ingest = false,
      }
    }
  }

  if (load_vcf_headers) {
    call ProcessVCFHeaders {
      input:
        variants_docker = effective_variants_docker,
        go = select_all(LoadDataViaBigQueryWriteAPI.done), # add a gate for Parquet header loading here once that's implemented
        dataset_name = dataset_name,
        project_id = project_id,
    }
  }

  if (load_vet_and_ref_ranges) {
    if (use_parquet_ingest) {
      String defined_parquet_output_dir = select_first([parquet_output_gcs_dir])

      # Set up lifecycle rules for parquet directories before loading
      if (configure_parquet_lifecycle) {
        call ConfigureParquetLifecycle {
          input:
            output_gcs_dir = defined_parquet_output_dir,
            billing_project_id = billing_project_id,
            variants_docker = effective_variants_docker,
        }
      }

      # Discover and load Parquet files into BigQuery after all data has been created.
      call DiscoverParquetFiles {
        input:
          output_gcs_dir = defined_parquet_output_dir,
          project_id = project_id,
          dataset_name = dataset_name,
          regular_table_prefixes = ["sample_chromosome_ploidy"],
          superpartitioned_table_prefixes = ["vet", "ref_ranges"],
          go = flatten([
            select_all([ConfigureParquetLifecycle.done]),
            select_all(GenerateParquetFilesFromInputGVCFs.done)
          ]),
          variants_docker = effective_variants_docker,
      }

      scatter (fofn in DiscoverParquetFiles.file_fofns) {
        call LoadParquetFilesToBQ {
          input:
            project_id = project_id,
            dataset_name = dataset_name,
            fofn_file = fofn,
            batch_size = 10000,
            variants_docker = effective_variants_docker,
        }
      }

      call VerifyParquetLoading {
        input:
          project_id = project_id,
          dataset_name = dataset_name,
          gcs_files_list = DiscoverParquetFiles.all_files_list,
          go = LoadParquetFilesToBQ.done,
          variants_docker = effective_variants_docker,
      }

      if (delete_parquet_files_after_loading && VerifyParquetLoading.all_loaded) {
        call DeleteParquetFiles {
          input:
            output_gcs_dir = defined_parquet_output_dir,
            use_alternate_delete_strategy = use_alternate_parquet_delete_strategy,
            billing_project_id = billing_project_id,
            cloud_sdk_docker = effective_cloud_sdk_docker,
        }
      }
    }

    call SetIsLoadedColumn {
      input:
        # A BQ Write API-flavored invocation of `LoadData` actually loads all data into vet and ref ranges tables, but a
        # Parquet-flavored invocation of `LoadData` only generates Parquet files from input gVCFs.
        # Because the loading of Parquet data into BigQuery is handled by a chain of WDL tasks subsequent to
        # `GenerateParquetFilesFromInputGVCFs`, the `go` trigger for setting the `is_loaded` column is the `done` output
        # of the last task in that chain, `VerifyParquetLoading`. The other component of the `go` trigger is the
        # `LoadDataViaBigQueryWriteAPI.done` corresponding to the Write API flow.
        # Intentionally using select_first to pick whichever of the two mutually exclusive code paths (Parquet vs WriteAPI) ran.
        #@ except: UnnecessaryFunctionCall
        go = select_all(select_first([[VerifyParquetLoading.done], LoadDataViaBigQueryWriteAPI.done])),
        project_id = project_id,
        dataset_name = dataset_name,
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }
  }

  output {
    Boolean done = true
    Boolean used_tighter_gcp_quotas = is_rate_limited_beta_customer
    String recorded_git_hash = effective_git_hash
    # Intentionally using select_first to pick the stderr files from whichever of the two mutually exclusive code paths (Parquet vs WriteAPI) ran.
    #@ except: UnnecessaryFunctionCall
    Array[File] load_data_stderrs = select_first([select_all(GenerateParquetFilesFromInputGVCFs.stderr), select_all(LoadDataViaBigQueryWriteAPI.stderr)])
    Boolean? parquet_loading_verified = VerifyParquetLoading.all_loaded
    Int? parquet_files_loaded = VerifyParquetLoading.loaded_files
    Int? parquet_total_files = VerifyParquetLoading.total_files
  }
}

task CreateFOFNs {
  input {
    Int batch_size
    File input_vcf_index_list
    File input_vcf_list
    File sample_name_list
    String basic_docker
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    split -a 5 -l ~{batch_size} ~{input_vcf_list} batched_vcfs.
    split -a 5 -l ~{batch_size} ~{input_vcf_index_list} batched_vcf_indexes.
    split -a 5 -l ~{batch_size} ~{sample_name_list} batched_sample_names.
  >>>
  runtime {
    docker: basic_docker
    bootDiskSizeGb: 15
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Array[File] vcf_batch_vcf_fofns = glob("batched_vcfs.*")
    Array[File] vcf_batch_vcf_index_fofns = glob("batched_vcf_indexes.*")
    Array[File] vcf_sample_name_fofns = glob("batched_sample_names.*")
  }
}

# This is the task known as `LoadData` on the ah_var_store branch, but on the Parquet branches it does not load data.
# In the Parquet flow we only generate Parquet files from input gVCFs and then stage them to GCS; the actual data
# loading is performed by a suite of other, Parquet-specific downstream tasks.
task ProcessInputGVCFs {
  input {
    Int index
    String dataset_name
    String project_id
    String? billing_project_id

    Array[File] input_vcf_indexes
    Array[File] input_vcfs
    File interval_list
    File sample_map
    Array[String] sample_names

    String? drop_state
    Boolean? drop_state_includes_greater_than = false
    Boolean force_loading_from_non_allele_specific = false
    Boolean skip_loading_vqsr_fields = false
    Boolean use_compressed_references = false
    Boolean load_vet_and_ref_ranges
    Boolean load_vcf_headers

    String? parquet_output_gcs_dir

    String gatk_docker
    File? gatk_override
    Int load_data_preemptible
    Int load_data_maxretries

    Boolean use_parquet_ingest
  }

  meta {
    description: "Generate Parquet files from input gVCFs OR load data into BigQuery using the Write API, depending on the value of `use_parquet_ingest`."
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }

    input_vcf_indexes: {
      localization_optional: true
    }
  }

  Int num_samples = length(sample_names)
  String temp_table = "~{dataset_name}.sample_names_to_load_~{index}"
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:import_genomes"
  String table_name = "sample_info"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # workaround for https://github.com/broadinstitute/cromwell/issues/3647
    export TMPDIR=/tmp

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    ## check which samples still need loading by looking in the BQ database for the loaded status of these samples

    # Create temp table with the sample_names and load external sample names into temp table -- make sure it doesn't exist already
    set +o errexit
    bq --apilog=false show --project_id=~{project_id} ~{temp_table} > /dev/null
    BQ_SHOW_RC=$?
    set -o errexit

    # If there is already a table of sample names or something else is wrong, burn it down to start fresh.
    if [ $BQ_SHOW_RC -eq 0 ]; then
      bq --apilog=false rm -t -f --project_id=~{project_id} ~{temp_table}
    fi

    echo "Creating the external sample name list table ~{temp_table}"
    bq --apilog=false --project_id=~{project_id} mk ~{temp_table} "sample_name:STRING"
    NAMES_FILE=~{write_lines(sample_names)}
    bq --apilog=false load --project_id=~{project_id} ~{temp_table} $NAMES_FILE "sample_name:STRING"

    # Get the current min/max id, or 0 if there are none. Withdrawn samples still have IDs so don't filter them out.
    # bq query --max_rows check: ok one row
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} '
      SELECT IFNULL(MIN(sample_id),0) as min, IFNULL(MAX(sample_id),0) as max FROM `~{dataset_name}.~{table_name}`
        AS samples JOIN `~{temp_table}` AS temp ON samples.sample_name = temp.sample_name' > results.csv

    # Get sample map of samples that haven't been loaded yet
    if [[ "~{load_vet_and_ref_ranges}" = "true" ]]
    then

    cat > query_vet_and_ref_ranges.sql <<'FIN_VET_REF'

      SELECT sample_id, samples.sample_name FROM `~{dataset_name}.~{table_name}` AS samples JOIN
      `~{temp_table}` AS temp ON
      samples.sample_name = temp.sample_name WHERE
      samples.sample_id NOT IN (
        SELECT DISTINCT ref.sample_id FROM
          `~{project_id}.~{dataset_name}.samples_with_reference_data` ref JOIN
          `~{project_id}.~{dataset_name}.samples_with_variant_data` vet USING (sample_id) JOIN
          `~{project_id}.~{dataset_name}.sample_chromosome_ploidy` ploidy USING(sample_id)
      ) AND
      samples.withdrawn IS NULL

    FIN_VET_REF

    cat query_vet_and_ref_ranges.sql |
      bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} \
        --max_rows ~{num_samples} > variant_and_reference_data.status_bucket.csv
    fi

    if [[ "~{load_vcf_headers}" = "true" ]]
    then

    cat > query_headers.sql <<'FIN_HEADERS'

      SELECT sample_id, samples.sample_name FROM `~{dataset_name}.~{table_name}` AS samples JOIN
      `~{temp_table}` AS temp ON
      samples.sample_name = temp.sample_name WHERE
      samples.sample_id NOT IN (
        SELECT sample_id FROM `~{project_id}.~{dataset_name}.samples_with_header_data`
      ) AND
      samples.withdrawn IS NULL
    FIN_HEADERS

    cat query_headers.sql |
      bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} \
        --max_rows ~{num_samples} > header_data.status_bucket.csv
    fi

    ## delete the table that was only needed for this ingest test
    bq --apilog=false --project_id=~{project_id} rm -f=true ~{temp_table}

    # If a given sample shows up in any status bucket it should appear in the final sample map exactly once.
    # Add a header manually:
    echo "sample_id,sample_name" > sample_map.csv
    # The real header sorts to the bottom of the file, delete that.
    cat *.status_bucket.csv | sort -u | sed '$d' >> sample_map.csv

    ## now we want to create a sub list of these samples (without the ones that have already been loaded)

    python3 /gatk/scripts/variantstore/scripts/curate_input_array_files.py \
      --sample_map_to_be_loaded_file_name sample_map.csv \
      --sample_name_list_file_name $NAMES_FILE \
      --vcf_list_file_name ~{write_lines(input_vcfs)} \
      --vcf_index_list_file_name  ~{write_lines(input_vcf_indexes)}

    # translate files created by the python script into BASH arrays---but only of the samples that aren't there already
    VCFS_ARRAY=($(cat output_vcf_list_file |tr "\n" " "))
    VCF_INDEXES_ARRAY=($(cat output_vcf_index_list_file |tr "\n" " "))
    SAMPLE_NAMES_ARRAY=($(cat output_sample_name_list_file |tr "\n" " "))

    # loop over the BASH arrays (See https://stackoverflow.com/questions/6723426/looping-over-arrays-printing-both-index-and-value)
    for i in "${!VCFS_ARRAY[@]}"; do
      gs_input_vcf="${VCFS_ARRAY[$i]}"
      gs_input_vcf_index="${VCF_INDEXES_ARRAY[$i]}"
      sample_name="${SAMPLE_NAMES_ARRAY[$i]}"

      # We always do our own localization.
      # It seems possible that the Parquet / non-Parquet branches below might be coalesced.
      if [[ "~{use_parquet_ingest}" = 'true' ]]
      then
        updated_input_vcf=input_vcf_${i}_${sample_name}.vcf.gz
        gcloud storage ~{"--billing-project " + billing_project_id} cp $gs_input_vcf $updated_input_vcf
        gcloud storage ~{"--billing-project " + billing_project_id} cp $gs_input_vcf_index ${updated_input_vcf}.tbi
      else
        gcloud storage ~{"--billing-project " + billing_project_id} cp $gs_input_vcf input_vcf_$i.vcf.gz
        gcloud storage ~{"--billing-project " + billing_project_id} cp $gs_input_vcf_index input_vcf_$i.vcf.gz.tbi
        updated_input_vcf=input_vcf_$i.vcf.gz
      fi

      gatk --java-options "-Xmx2g" CreateVariantIngestFiles \
        -V ${updated_input_vcf} \
        -L ~{interval_list} \
        ~{"--ref-block-gq-to-ignore " + drop_state} \
        --ignore-above-gq-threshold ~{drop_state_includes_greater_than} \
        --force-loading-from-non-allele-specific ~{force_loading_from_non_allele_specific} \
        --project-id ~{project_id} \
        --dataset-name ~{dataset_name} \
        --output-type ~{true="PARQUET" false="BQ" use_parquet_ingest} \
        --enable-reference-ranges ~{load_vet_and_ref_ranges} \
        --enable-vet ~{load_vet_and_ref_ranges} \
        -SN ${sample_name} \
        -SNM ~{sample_map} \
        --ref-version 38 \
        --skip-loading-vqsr-fields ~{skip_loading_vqsr_fields} \
        --enable-vcf-headers ~{load_vcf_headers} \
        --use-compressed-refs ~{use_compressed_references}

      # The Parquet / non-Parquet branches here might also be coalesced.
      if [[ "~{use_parquet_ingest}" = 'true' ]]
      then
        rm $updated_input_vcf
        rm ${updated_input_vcf}.tbi

        OUTPUT_GCS_DIR=$(echo ~{parquet_output_gcs_dir} | sed 's/\/$//')
        # the file name is a little wonky, so let's just grab the file using such a star statement
        vet_parquet_file=`ls vet_*.parquet`
        ref_parquet_file=`ls ref_*.parquet`
        ploidy_parquet_file=`ls sample_chromosome_ploidy_*.parquet`

        # parse the table superpartition out of the file name
        table_number=$(echo "$vet_parquet_file" | cut -d'_' -f2)

        # copy the vet and ref parquet files to the gcs bucket in the right place
        gcloud storage ~{"--billing-project " + billing_project_id} cp $vet_parquet_file ${OUTPUT_GCS_DIR}/vet/$table_number/$vet_parquet_file
        gcloud storage ~{"--billing-project " + billing_project_id} cp $ref_parquet_file ${OUTPUT_GCS_DIR}/ref_ranges/$table_number/$ref_parquet_file
        gcloud storage ~{"--billing-project " + billing_project_id} cp $ploidy_parquet_file ${OUTPUT_GCS_DIR}/sample_chromosome_ploidy/$ploidy_parquet_file

        # cleanup after ourselves
        rm *.parquet
      else
        rm input_vcf_$i.vcf.gz
        rm input_vcf_$i.vcf.gz.tbi
      fi

    done
  >>>

  runtime {
    docker: gatk_docker
    maxRetries: load_data_maxretries
    memory: "3.75 GB"
    disks: "local-disk 50 HDD"
    preemptible: load_data_preemptible
    cpu: 1
    noAddress: true
  }
  output {
    Boolean done = true
    File stderr = stderr()
  }
}

task ProcessVCFHeaders {
  input {
    String dataset_name
    String project_id
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Array[Boolean] go
    String variants_docker
  }
  meta {
    volatile: true
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    python3 /app/process_sample_vcf_headers.py \
      --project_id=~{project_id} \
      --dataset_name ~{dataset_name}
  >>>

  runtime {
    docker: variants_docker
    disks: "local-disk 500 HDD"
  }
}

task SetIsLoadedColumn {
  input {
    String dataset_name
    String project_id

    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Array[Boolean] go
    String cloud_sdk_docker
  }
  meta {
    # Always run. This task is idempotent and depends on upstream tasks side-effecting data into BigQuery.
    volatile: true
  }

  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:import_genomes"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # Note that we tried modifying CreateVariantIngestFiles to UPDATE sample_info.is_loaded on a per-sample basis.
    # The major issue that was found is that BigQuery allows only 20 such concurrent DML statements. Considered using
    # an exponential backoff, but at the number of samples that are being loaded this would introduce significant delays
    # in workflow processing. So this method is used to set *all* of the sample_info.is_loaded flags at one time.

    # bq query --max_rows check: ok update
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} '

    UPDATE `~{project_id}.~{dataset_name}.sample_info`
    SET is_loaded = TRUE
    WHERE
      sample_id IN (
      SELECT sample_id FROM `~{project_id}.~{dataset_name}.samples_with_all_data`
    );

    '
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "4 GB"
    disks: "local-disk 500 HDD"
    cpu: 1
  }

  output {
    Boolean done = true
  }
}

task GetUningestedSampleIds {
  input {
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go
    String dataset_name
    String project_id

    File external_sample_names
    Int num_samples
    String table_name
    String cloud_sdk_docker
    # At least one of these "load" inputs must be true
    Boolean load_vet_and_ref_ranges
    Boolean load_vcf_headers
  }
  meta {
    # Do not call cache this, we want to read the database state every time.
    volatile: true
  }

  Int samples_per_table = 4000
  # add labels for DSP Cloud Cost Control Labeling and Reporting
  String bq_labels = "--label service:gvs --label team:variants --label managedby:import_genomes"
  String temp_table="~{dataset_name}.sample_names_to_load"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    echo "project_id = ~{project_id}" > ~/.bigqueryrc

    # Create temp table with the sample_names and load external sample names into temp table
    # Make this idempotent - clean up any existing temp table from failed runs
    set +o errexit
    bq --apilog=false show --project_id=~{project_id} ~{temp_table} > /dev/null
    BQ_SHOW_RC=$?
    set -o errexit

    # If temp table already exists, clean it up (idempotent behavior for retries)
    if [ $BQ_SHOW_RC -eq 0 ]; then
      echo "Temp table ~{temp_table} already exists from previous run, cleaning up"
      bq --apilog=false --project_id=~{project_id} rm -f ~{temp_table}
    fi

    echo "Creating the external sample name list table ~{temp_table}"
    bq --apilog=false --project_id=~{project_id} mk ~{temp_table} "sample_name:STRING"
    bq --apilog=false load --project_id=~{project_id} ~{temp_table} ~{external_sample_names} "sample_name:STRING"

    # Get the current min/max id, or 0 if there are none. Withdrawn samples still have IDs so don't filter them out.
    # bq query --max_rows check: ok one row
    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} '

      SELECT IFNULL(MIN(sample_id),0) as min, IFNULL(MAX(sample_id),0) as max FROM `~{dataset_name}.~{table_name}`
        AS samples JOIN `~{temp_table}` AS temp ON samples.sample_name = temp.sample_name

    ' > results.csv

    # prep for being able to return min table id
    min_sample_id=$(tail -1 results.csv | cut -d, -f1)
    max_sample_id=$(tail -1 results.csv | cut -d, -f2)

    # no samples have been loaded or we don't have the right external_sample_names or something else is wrong, bail
    if [ $max_sample_id -eq 0 ]; then
      echo "Max id is 0. Exiting"
      exit 1
    fi

    python3 -c "from math import ceil; print(ceil($max_sample_id/~{samples_per_table}))" > max_sample_id
    python3 -c "from math import ceil; print(ceil($min_sample_id/~{samples_per_table}))" > min_sample_id

    # Get sample map of samples that haven't been loaded yet
    # Break out individual queries into "status buckets" for all of the statuses we care about.

    if [[ "~{load_vet_and_ref_ranges}" = "true" ]]
    then

    cat > query_vet_and_ref_ranges.sql <<'FIN_VET_REF'

      SELECT sample_id, samples.sample_name FROM `~{dataset_name}.~{table_name}` AS samples JOIN
      `~{temp_table}` AS temp ON
      samples.sample_name = temp.sample_name WHERE
      samples.sample_id NOT IN (
        SELECT DISTINCT ref.sample_id FROM
          `~{project_id}.~{dataset_name}.samples_with_reference_data` ref JOIN
          `~{project_id}.~{dataset_name}.samples_with_variant_data` vet USING (sample_id) JOIN
          `~{project_id}.~{dataset_name}.sample_chromosome_ploidy` ploidy USING(sample_id)
      ) AND
      samples.withdrawn IS NULL

    FIN_VET_REF

    cat query_vet_and_ref_ranges.sql |
      bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} \
        --max_rows ~{num_samples} > variant_and_reference_data.status_bucket.csv
    fi

    if [[ "~{load_vcf_headers}" = "true" ]]
    then

    cat > query_headers.sql <<'FIN_HEADERS'

      SELECT sample_id, samples.sample_name FROM `~{dataset_name}.~{table_name}` AS samples JOIN
      `~{temp_table}` AS temp ON
      samples.sample_name = temp.sample_name WHERE
      samples.sample_id NOT IN (
        SELECT sample_id FROM `~{project_id}.~{dataset_name}.samples_with_header_data`
      ) AND
      samples.withdrawn is NULL
    FIN_HEADERS

    cat query_headers.sql |
      bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} \
        --max_rows ~{num_samples} > header_data.status_bucket.csv
    fi

    ## delete the table that was only needed for this ingest test
    bq --apilog=false --project_id=~{project_id} rm -f=true ~{temp_table}

    # If a given sample shows up in any status bucket it should appear in the final sample map exactly once.
    # Add a header manually:
    echo "sample_id,sample_name" > sample_map.csv
    # The real header sorts to the bottom of the file, delete that.
    cat *.status_bucket.csv | sort -u | sed '$d' >> sample_map.csv

    cut -d, -f1 sample_map.csv > gvs_ids.csv

    ## delete the table that was only needed for this ingest
    bq --apilog=false --project_id=~{project_id} rm -f=true ~{temp_table}
  >>>
  runtime {
    docker: cloud_sdk_docker
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 5
    cpu: 1
  }
  output {
    Int max_table_id = ceil(read_float("max_sample_id"))
    Int min_table_id = ceil(read_float("min_sample_id"))
    File sample_map = "sample_map.csv"
    File gvs_ids = "gvs_ids.csv"
    Array[File] status_buckets = glob("*.status_bucket.csv")
    Array[File] queries = glob("query_*.sql")
  }
}

task CurateInputLists {
  input {
    File input_vcf_index_list
    File input_vcf_list
    File input_samples_to_be_loaded_map
    File input_sample_name_list
    String variants_docker
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    python3 /app/curate_input_array_files.py --sample_map_to_be_loaded_file_name ~{input_samples_to_be_loaded_map} \
                                             --sample_name_list_file_name ~{input_sample_name_list} \
                                             --vcf_list_file_name ~{input_vcf_list} \
                                             --vcf_index_list_file_name  ~{input_vcf_index_list}
  >>>
  runtime {
    docker: variants_docker
    memory: "3 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    cpu: 1
  }

  output {
    File input_vcf_indexes = "output_vcf_index_list_file"
    File input_vcfs = "output_vcf_list_file"
    File sample_name_list = "output_sample_name_list_file"
  }
}

task CreateSampleDataViews {
  input {
    String project_id
    String dataset_name
    String cloud_sdk_docker
  }

  String bq_labels = "--label service:gvs --label team:variants --label managedby:import_genomes"


  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    cat > query.sql <<'FIN'

      -- Because the vet and ref_ranges tables are partitioned by sample_id, their INFORMATION_SCHEMA partition ids
      -- will be stringified sample ids. These views identify which samples have vet, reference, or header data loaded.
      --
      -- The Parquet flow is not currently writing sample status rows so we use the data in INFORMATION_SCHEMA to
      -- determine load status. Conversely, data written with the BigQuery Write API seems to result in very delayed
      -- population of INFORMATION_SCHEMA, often lagging writes by several hours, which makes reading INFORMATION_SCHEMA
      -- unreliable with the Write API. The following vet and ref ranges queries UNION DISTINCT the sample_load_status
      -- table with INFORMATION_SCHEMA to reliably detect sample data regardless of how it was loaded into GVS.
      --
      -- This code also provides for a header row existence view if headers are being loaded.

      DECLARE sample_load_status_template STRING;

      -- In the future the `sample_load_status` table may no longer be needed. Only refer to `sample_load_status` in the
      -- following existence queries if the table actually exists.
      DECLARE sample_load_status_table_exists INT64;

      SET sample_load_status_template = """

        UNION DISTINCT

        SELECT sample_id FROM
        `~{project_id}.~{dataset_name}.sample_load_status`
        WHERE status = '%s'

      """;

      SET sample_load_status_table_exists = (
        SELECT COUNT(1) FROM
        `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES`
        WHERE table_name = 'sample_load_status'
      );

      BEGIN
      DECLARE variants_load_status_clause STRING;
      DECLARE create_variant_data_view STRING;

      IF sample_load_status_table_exists > 0 THEN
        SET variants_load_status_clause = format(sample_load_status_template, 'VARIANTS_LOADED');
      ELSE
        SET variants_load_status_clause = '';
      END IF;

      SET create_variant_data_view = """

        CREATE OR REPLACE VIEW `~{project_id}.~{dataset_name}.samples_with_variant_data` AS
        (
          SELECT CAST(partition_id AS INT64) AS sample_id
          FROM `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS`
          WHERE
          NOT STARTS_WITH(partition_id, '__') AND total_logical_bytes > 0 AND REGEXP_CONTAINS(table_name, '^vet_[0-9]+$')

      """ || variants_load_status_clause || ");";

      -- debug
      SELECT create_variant_data_view;
      EXECUTE IMMEDIATE create_variant_data_view;
      END;

      BEGIN
      DECLARE references_load_status_clause STRING;
      DECLARE create_reference_data_view STRING;

      IF sample_load_status_table_exists > 0 THEN
        SET references_load_status_clause = format(sample_load_status_template, 'REFERENCES_LOADED');
      ELSE
        SET references_load_status_clause = '';
      END IF;

      SET create_reference_data_view = """

      CREATE OR REPLACE VIEW `~{project_id}.~{dataset_name}.samples_with_reference_data` AS
      (
        SELECT CAST(partition_id AS INT64) AS sample_id
        FROM `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.PARTITIONS`
        WHERE
        NOT STARTS_WITH(partition_id, '__') AND total_logical_bytes > 0 AND REGEXP_CONTAINS(table_name, '^ref_ranges_[0-9]+$')

      """ || references_load_status_clause || ");";

      -- debug
      SELECT create_reference_data_view;
      EXECUTE IMMEDIATE create_reference_data_view;
      END;

      -- The header data view is created conditionally as the header tables are created conditionally.
      BEGIN

      DECLARE header_table_exists INT64;
      DECLARE query_header_existence_clause STRING;
      DECLARE headers_load_status_clause STRING;
      DECLARE create_header_data_view STRING;
      DECLARE create_all_sample_data_view STRING;

      SET header_table_exists = (
        SELECT COUNT(1) FROM
        `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES`
        WHERE table_name = 'vcf_header_lines_scratch'
      );

      IF header_table_exists > 0 THEN
        IF sample_load_status_table_exists > 0 THEN
          SET headers_load_status_clause = format(sample_load_status_template, 'HEADERS_LOADED');
        ELSE
          SET headers_load_status_clause = '';
        END IF;

        SET create_header_data_view = """

          CREATE OR REPLACE VIEW `~{project_id}.~{dataset_name}.samples_with_header_data` AS
          (
          SELECT DISTINCT sample_id FROM `~{project_id}.~{dataset_name}.vcf_header_lines_scratch`

          """ || headers_load_status_clause || ");";

        -- debug
        SELECT create_header_data_view;
        EXECUTE IMMEDIATE create_header_data_view;

        SET query_header_existence_clause = """

        JOIN `~{project_id}.~{dataset_name}.samples_with_header_data` header USING(sample_id)

        """;
      ELSE
        SET query_header_existence_clause = '';
      END IF;

      SET create_all_sample_data_view = """

        CREATE OR REPLACE VIEW `~{project_id}.~{dataset_name}.samples_with_all_data` AS
        (
          SELECT DISTINCT ref.sample_id FROM
            `~{project_id}.~{dataset_name}.samples_with_reference_data` ref JOIN
            `~{project_id}.~{dataset_name}.samples_with_variant_data` vet USING (sample_id) JOIN
            `~{project_id}.~{dataset_name}.sample_chromosome_ploidy` ploidy USING(sample_id)
      """ || query_header_existence_clause || """
        );
      """;

      EXECUTE IMMEDIATE create_all_sample_data_view;
      END;

    FIN

    bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} < query.sql

  >>>

  runtime {
    docker: cloud_sdk_docker
    memory: "4 GB"
    disks: "local-disk 500 HDD"
    preemptible: 3
    cpu: 1
  }

  output {
    Boolean done = true
    File query = "query.sql"
  }
}

task ConfigureParquetLifecycle {
  input {
    String output_gcs_dir
    # TODO: billing_project_id is declared but not passed to load_parquet_to_bq.py - see VS-1955.
    String? billing_project_id
    String variants_docker
  }

  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    # Extract bucket name from GCS path
    BUCKET_NAME=$(echo ~{output_gcs_dir} | sed 's|gs://||' | cut -d'/' -f1)

    # Extract bucket path prefix (if any) to ensure lifecycle rules are applied to the correct subdirectories
    # For example, if output_gcs_dir is gs://my-bucket/path/to/data/, we want the prefix to be path/to/data/ to apply rules to that subdirectory rather than the whole bucket
    # First, remove gs:// prefix and trailing slash, then check if there's a path component
    TEMP_PATH=$(echo ~{output_gcs_dir} | sed 's|gs://||' | sed 's/\/$//')

    # If TEMP_PATH contains a /, extract everything after the first /, otherwise set to empty
    if [[ "$TEMP_PATH" == */* ]]; then
      BUCKET_PATH_PREFIX=$(echo "$TEMP_PATH" | cut -d'/' -f2-)
      # Append trailing slash since we have a path
      BUCKET_PATH_PREFIX="${BUCKET_PATH_PREFIX}/"
    else
      BUCKET_PATH_PREFIX=""
    fi

    echo "Configuring lifecycle for bucket: ${BUCKET_NAME}"
    echo "Path prefix: '${BUCKET_PATH_PREFIX}'"

    # Get existing lifecycle configuration if any
    set +e
    gcloud storage buckets describe gs://${BUCKET_NAME} \
      ~{"--billing-project " + billing_project_id} \
      --format="json(lifecycle_config)" > existing_lifecycle.json 2>/dev/null
    EXISTING_RC=$?
    set -e

    if [ $EXISTING_RC -ne 0 ]; then
      echo "Error encountered retrieving lifecycle rules for bucket $BUCKET_NAME - does that bucket exist?"
      exit 1;
    fi

    # Create the new lifecycle *rule* for parquet directories
    cat > new_lifecycle_rule.json << EOF
    {
      "action": {"type": "Delete"},
      "condition": {
        "age": 14,
        "matchesPrefix": ["${BUCKET_PATH_PREFIX}vet/", "${BUCKET_PATH_PREFIX}ref_ranges/", "${BUCKET_PATH_PREFIX}sample_chromosome_ploidy/"]
      }
    }
EOF

    # If here, we successfully found a lifecycle config (even if it's empty), check if it's empty or null
    if [ -s existing_lifecycle.json ] && [ "$(cat existing_lifecycle.json)" != "null" ]; then
      # Note: The gcloud command returns lifecycle_config with a key of "lifecycle_config" but the gcloud buckets update command expects the key to be "lifecycle", so we need to rename that key before merging with jq
      jq '{lifecycle: .lifecycle_config}' existing_lifecycle.json > temp.json
      mv temp.json existing_lifecycle.json
    else
      echo "No existing lifecycle configuration found (file is empty or contains the string 'null'), starting with empty lifecycle configuration"
      # Create the new lifecycle configuration with no rules (we'll add the rule further on) for parquet directories
      cat > existing_lifecycle.json << EOF
      {
        "lifecycle": {
          "rule": [
          ]
        }
      }
EOF
    fi

    # Now use jq to merge the new lifecycle rule with the existing lifecycle configuration,
    # but only add it if there isn't already a rule with the same condition.matchesPrefix values
    jq --slurpfile new_rule new_lifecycle_rule.json '
      . as $cfg
      | $new_rule[0] as $nr
      | ($cfg.lifecycle.rule // []) as $rules
      | if ($rules | any(.condition.matchesPrefix == $nr.condition.matchesPrefix))
        then $cfg
        else $cfg | .lifecycle.rule += [$nr]
        end
    ' existing_lifecycle.json > updated_lifecycle.json

    # Apply the updated lifecycle configuration
    gcloud storage buckets update gs://${BUCKET_NAME} \
      ~{"--billing-project " + billing_project_id} \
      --lifecycle-file=updated_lifecycle.json

    echo "✓ Lifecycle rule applied: After 14 days, it will delete files in the bucket: ${BUCKET_NAME}, with path prefixes ${BUCKET_PATH_PREFIX}vet/, ${BUCKET_PATH_PREFIX}ref_ranges/ and ${BUCKET_PATH_PREFIX}/sample_chromosome_ploidy"
  >>>

  runtime {
    docker: variants_docker
    memory: "1 GB"
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
    Array[String] regular_table_prefixes
    Array[String] superpartitioned_table_prefixes
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Array[Boolean] go
    String? billing_project_id
    String variants_docker
  }

  meta {
    volatile: true
  }

  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    # Normalize GCS path to ensure exactly one trailing slash
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    # List all objects, filter for Parquet files
    echo "Listing files in ${OUTPUT_GCS_DIR}..."
    gcloud storage ls --recursive ~{"--billing-project " + billing_project_id} \
      "${OUTPUT_GCS_DIR}/" > all_objects.txt || true

    grep '\.parquet$' all_objects.txt > all_files.txt || touch all_files.txt

    FILE_COUNT=$(wc -l < all_files.txt)
    echo "Found $FILE_COUNT Parquet files"

    # Parse and group by table
    python3 /app/parse_and_group_files.py \
      --input all_files.txt \
      --output-dir grouped_files \
      --project-id ~{project_id} \
      --dataset ~{dataset_name} \
      --superpartitioned-table-prefixes ~{sep=" " superpartitioned_table_prefixes} \
      --regular-table-prefixes ~{sep=" " regular_table_prefixes}
  >>>

  runtime {
    docker: variants_docker
    memory: "4 GB"
    disks: "local-disk 50 HDD"
    preemptible: 3
    cpu: 2
  }

  output {
    Array[File] file_fofns = glob("grouped_files/*.fofn")
    File all_files_list = "all_files.txt"
    File stats_json = "grouped_files/stats.json"
  }
}

task LoadParquetFilesToBQ {
  input {
    String project_id
    String dataset_name
    File fofn_file
    Int batch_size
    # TODO: billing_project_id is declared but not passed to load_parquet_to_bq.py - see VS-1955.
    String? billing_project_id
    String variants_docker
  }

  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail
    # Table name is extracted from FOFN filename by the Python script
    python3 /app/load_parquet_to_bq.py \
      --project-id ~{project_id} \
      --dataset-name ~{dataset_name} \
      --files-fofn ~{fofn_file} \
      --batch-size ~{batch_size} \
      --output-stats stats.json
  >>>

  runtime {
    docker: variants_docker
    memory: "4 GB"
    disks: "local-disk 20 HDD"
    preemptible: 5
    maxRetries: 3
    cpu: 1
  }

  output {
    Boolean done = true
    File stats_json = "stats.json"
  }
}

task VerifyParquetLoading {
  input {
    String project_id
    String dataset_name
    File gcs_files_list
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Array[Boolean] go
    String variants_docker
  }

  meta {
    volatile: true
  }

  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail
    mkdir -p verification_output

    python3 /app/verify_all_loaded.py \
      --project-id ~{project_id} \
      --dataset-name ~{dataset_name} \
      --gcs-files-list ~{gcs_files_list} \
      --output-dir verification_output
  >>>

  runtime {
    docker: variants_docker
    memory: "4 GB"
    disks: "local-disk 20 HDD"
    cpu: 1
  }

  output {
    # TODO: Sprocket flags read_json indexing as invalid on Union type; fix by upgrading to WDL 1.1 and using struct coercion — see VS-1957.
    File results_json = "verification_output/verification_results.json"
    Boolean all_loaded = read_json(results_json)["all_loaded"]
    Int total_files = read_json(results_json)["total_files"]
    Int loaded_files = read_json(results_json)["loaded_files"]
    Int missing_files = read_json(results_json)["missing_files"]
    File? missing_files_list = "verification_output/missing_files.txt"
    Boolean done = true
  }
}

task DeleteParquetFiles {
  input {
    String output_gcs_dir
    Boolean use_alternate_delete_strategy = false

    String? billing_project_id
    String cloud_sdk_docker
  }

  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    # Normalize GCS path by removing any trailing slash
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    if [ "~{use_alternate_delete_strategy}" = "false" ]; then
      gcloud storage rm --recursive ~{"--billing-project " + billing_project_id} "${OUTPUT_GCS_DIR}/"'**/*.parquet'
    else
      # List the contents of the vet and ref_ranges directories for subsequent deletion in the loop below
      echo "Listing directories under ${OUTPUT_GCS_DIR}/vet/ and ${OUTPUT_GCS_DIR}/ref_ranges/ ${OUTPUT_GCS_DIR}/sample_chromosome_ploidy/ for deletion..."
      gcloud storage ls ~{"--billing-project " + billing_project_id} \
        "${OUTPUT_GCS_DIR}/vet/" "${OUTPUT_GCS_DIR}/ref_ranges/" > parquet_dirs.txt
      echo "${OUTPUT_GCS_DIR}/sample_chromosome_ploidy/" >> parquet_dirs.txt

      # Iterate over all Google Cloud paths in parquet_dirs.txt and delete all objects therein
      echo "Deleting Parquet files..."
      while IFS= read -r gcs_path; do
        if [ -n "$gcs_path" ]; then
          echo "Deleting objects in: $gcs_path"
          gcloud storage rm ~{"--billing-project " + billing_project_id} "$gcs_path" --recursive
        fi
      done < parquet_dirs.txt
    fi

    echo "✓ Completed deletion of Parquet files."

  >>>
  output {
    Boolean done = true
  }

  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 100 HDD"
    preemptible: 3
    cpu: 1
  }
}
