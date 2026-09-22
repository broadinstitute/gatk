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

    # Independent post-load structural checks (VS-1989). The vet duplication screen flags any sample
    # whose vet row count is >= parquet_vet_duplication_threshold times the callset median; the
    # truncation screen flags any sample at or below median / parquet_vet_truncation_threshold. A flag
    # does not fail the load: the flagged samples' Parquet is moved to a quarantine prefix out of the
    # bulk delete's reach and the rest of the callset is deleted as normal. Set
    # parquet_allow_flagged_vet_loads to waive both screens and delete everything anyway.
    Float parquet_vet_duplication_threshold = 1.6
    # The two thresholds are equal by default but deliberately separate, because the evidence behind
    # them is not. Foxtrot calibration measured the high side only (1.6x flagged 2 vet samples of
    # 540,545); nothing has been measured below the median, and a variant-count distribution has no
    # reason to be symmetric -- its upper tail is bounded by biology while its lower tail absorbs
    # low-coverage samples and more aggressive GQ dropping. So the low side must be movable, and
    # switchable off, without disturbing the calibrated high side. Set to 0 to disable the truncation
    # screen entirely while leaving the duplication screen running.
    Float parquet_vet_truncation_threshold = 1.6
    Boolean parquet_allow_flagged_vet_loads = false
    # Whether a quarantine should also abort the workflow. Without this a flagged run is simply a green
    # run that quietly set some Parquet aside, and the only trace is a workflow output nobody reads.
    # Aborting is safe here because it happens after the files are already quarantined, so it destroys
    # nothing and reverses nothing -- it is purely a notification that the run needs a human.
    Boolean parquet_fail_on_quarantine = true
    # If set, the exact per-sample ploidy row count to validate against (e.g. 24 for WGS) instead of
    # the callset mode. Leave unset to infer the reference from the data (correct for exome/BGE/chrM).
    Int? parquet_expected_ploidy_rows_per_sample

    Boolean is_wgs = true
  }

  parameter_meta {
    # VS-1989 independent post-load structural checks; documented so these verification controls are
    # discoverable (womtool inputs, Terra, integration tests) alongside use_parquet_ingest.
    parquet_vet_duplication_threshold: "VS-1989 post-load verification: ratio-to-callset-median at or above which a vet sample's row count is flagged as a possible duplicate. Must be > 1; default 1.6, calibrated against Foxtrot."
    parquet_vet_truncation_threshold: "VS-1989 post-load verification: ratio whose reciprocal sets the low-side floor -- a vet sample at or below median/ratio is flagged as possibly truncated. Must be > 1, or 0 to disable the truncation screen; default 1.6 (i.e. 0.625x). Separate from parquet_vet_duplication_threshold because only the high side has been calibrated, so the low side can be retuned or switched off on its own."
    parquet_allow_flagged_vet_loads: "VS-1989 post-load verification: when false (default), the Parquet of any sample a vet duplication- or truncation-screen flag names is moved to a quarantine prefix instead of deleted (the load itself still succeeds, and the unflagged samples' Parquet is deleted as normal); when true the screens are waived and everything is deleted despite a flag. Family completeness and ploidy cardinality are exact checks that always gate load completeness regardless."
    parquet_fail_on_quarantine: "VS-1989 post-load verification: when true (default), a run that quarantined the Parquet of a duplication-flagged sample aborts after the quarantine completes, so the run is not silently green. Set false to leave the quarantine advisory (reported only through the parquet_quarantined_* outputs and the quarantine directory's README). Truncation-only flags never abort, because that threshold is not yet calibrated."
    parquet_expected_ploidy_rows_per_sample: "VS-1989 post-load verification: exact per-sample sample_chromosome_ploidy row count to validate against (e.g. 24 for WGS) instead of the inferred callset mode; leave unset to infer from the data (correct for exome/BGE/chrM)."
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

  # Compute effective scatter width, clamped to not exceed number of samples
  Int effective_scatter_width = if (defined(load_data_scatter_width)) then
                                  if select_first([load_data_scatter_width]) < num_samples then select_first([load_data_scatter_width]) else num_samples
                                else num_samples

  # Compute batch size using ceiling division to ensure we don't exceed the requested scatter width.
  # Using (x + y - 1) / y implementation of ceil(x/y)
  # Example: 10 samples, scatter_width 6 -> batch_size = ceil(10/6) = 2 -> actual tasks = ceil(10/2) ≈ 5 ≤ 6
  Int effective_load_data_batch_size = if (defined(load_data_scatter_width)) then
                                         (num_samples + effective_scatter_width - 1) / effective_scatter_width
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

  # Load whatever Parquet we generated into BigQuery. Headers and vet/ref/ploidy are discovered and
  # loaded together here; header data lands in the vcf_header_lines_scratch table. A caller wanting the
  # AoU-style gated phasing (load + sanity-check headers, THEN load the rest) does so by running this
  # workflow twice -- once with load_vcf_headers=true/load_vet_and_ref_ranges=false, then the reverse.
  # See the VS-1968 design doc, section 1.5. (No intra-run gate between the two phases is implemented.)
  if (use_parquet_ingest && (load_vet_and_ref_ranges || load_vcf_headers)) {
    String defined_parquet_output_dir = select_first([parquet_output_gcs_dir])

    # Table prefixes to discover/load/verify. vcf_header_lines_scratch is included only when headers
    # were generated (its table may not exist otherwise). vet/ref_ranges/ploidy prefixes are harmless
    # to include for a headers-only run -- no such files will be present to match.
    Array[String] parquet_regular_prefixes = if load_vcf_headers
      then ["sample_chromosome_ploidy", "vcf_header_lines_scratch"]
      else ["sample_chromosome_ploidy"]
    Array[String] parquet_superpartitioned_prefixes = ["vet", "ref_ranges"]

    # Parquet belonging to a sample the heuristic vet screens flag is moved under this subdirectory of
    # the output dir instead of being deleted, so it survives for inspection and re-ingest while the
    # rest of the callset's Parquet is deleted normally (VS-1989). Two tasks have to agree on the name
    # -- QuarantineFlaggedParquetFiles writes it and DiscoverParquetFiles must not re-discover it -- so
    # it is declared once here and threaded to both.
    #
    # The name is deliberately outside the prefixes ConfigureParquetLifecycle matches (vet/,
    # ref_ranges/, sample_chromosome_ploidy/, vcf_header_lines_scratch/) and outside the directory list
    # DeleteParquetFiles' alternate strategy walks. That is what exempts quarantined files from the
    # bucket's own 14-day Delete rule and from that strategy. Any lifecycle rule the operator has
    # configured on the bucket independently of this workflow is of course still theirs to reckon with.
    String parquet_quarantine_subdir = "quarantine"

    # Appended to each quarantined object's name so it no longer ends in ".parquet". The two deletion
    # strategies in DeleteParquetFiles are defeated by different halves of this: the default strategy's
    # whole-output-dir "*.parquet" glob by the rename, the alternate strategy's per-table-directory
    # deletes by the location. DeleteParquetFiles asserts this value does not itself end in ".parquet".
    String parquet_quarantine_suffix = ".quarantined"

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
        regular_table_prefixes = parquet_regular_prefixes,
        superpartitioned_table_prefixes = parquet_superpartitioned_prefixes,
        quarantine_subdir = parquet_quarantine_subdir,
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
        regular_table_prefixes = parquet_regular_prefixes,
        superpartitioned_table_prefixes = parquet_superpartitioned_prefixes,
        vet_duplication_threshold = parquet_vet_duplication_threshold,
        vet_truncation_threshold = parquet_vet_truncation_threshold,
        allow_flagged_vet_loads = parquet_allow_flagged_vet_loads,
        expected_ploidy_rows_per_sample = parquet_expected_ploidy_rows_per_sample,
        verification_diagnostics_gcs_dir = defined_parquet_output_dir + "/verification_diagnostics",
        billing_project_id = billing_project_id,
        go = LoadParquetFilesToBQ.done,
        variants_docker = effective_variants_docker,
    }

    # Move the flagged samples' Parquet out of the bulk delete's reach before that delete runs. Called
    # unconditionally, and a no-op on the usual clean run where the work list is empty, so that the
    # ordering guarantee is structural (DeleteParquetFiles consumes this task's `done`) rather than a
    # boolean anyone could get wrong. A failed quarantine followed by a bulk delete is the one unsafe
    # ordering, so this task fails loudly and takes the delete down with it.
    #
    # Not gated on delete_parquet_files_after_loading, because this workflow's own delete is not the
    # only thing that removes the Parquet: the lifecycle rule ConfigureParquetLifecycle installs deletes
    # it after 14 days regardless. Quarantining is what a flagged sample needs in either case.
    call QuarantineFlaggedParquetFiles {
      input:
        output_gcs_dir = defined_parquet_output_dir,
        quarantine_subdir = parquet_quarantine_subdir,
        quarantine_suffix = parquet_quarantine_suffix,
        quarantine_files_list = VerifyParquetLoading.quarantine_files_list,
        billing_project_id = billing_project_id,
        cloud_sdk_docker = effective_cloud_sdk_docker,
    }

    if (delete_parquet_files_after_loading && VerifyParquetLoading.safe_to_delete_parquet) {
      call DeleteParquetFiles {
        input:
          output_gcs_dir = defined_parquet_output_dir,
          quarantine_suffix = parquet_quarantine_suffix,
          use_alternate_delete_strategy = use_alternate_parquet_delete_strategy,
          billing_project_id = billing_project_id,
          go = QuarantineFlaggedParquetFiles.done,
          cloud_sdk_docker = effective_cloud_sdk_docker,
      }
    }

    # Make a quarantine visible. Everything upstream of here succeeded -- the samples are loaded and
    # their Parquet is safe -- so without this the run is green and the only trace is a workflow output
    # nobody reads. Aborting here is therefore a notification and nothing more: it destroys nothing,
    # rolls nothing back, and re-running with parquet_fail_on_quarantine = false (or
    # parquet_allow_flagged_vet_loads = true, if the flag has been reviewed and dismissed) picks up from
    # a fully loaded dataset.
    #
    # Conditioned on QuarantineFlaggedParquetFiles.quarantined_files rather than on the equivalent
    # VerifyParquetLoading count so the data dependency places the abort after the move; see that
    # output's comment. DeleteParquetFiles depends on the same task, so it and this abort are siblings
    # and Cromwell may start the delete first. That race is benign -- the delete skips the quarantined
    # files by both their suffix and their location, and deleting the unflagged samples' Parquet is
    # wanted either way.
    #
    # Only duplication flags abort. A truncation-only flag still quarantines and still reports, but the
    # truncation threshold has been measured on its high side only, and an uncalibrated heuristic should
    # not be able to fail a completed 500k-sample ingest.
    if (parquet_fail_on_quarantine && QuarantineFlaggedParquetFiles.quarantined_files > 0 && VerifyParquetLoading.quarantined_duplication_samples > 0) {
      call Utils.TerminateWorkflow as ParquetWasQuarantined {
        input:
          message = "Parquet ingest completed and all data loaded, but the VS-1989 vet duplication screen flagged " +
                    VerifyParquetLoading.quarantined_duplication_samples + " sample(s). Their Parquet (" +
                    QuarantineFlaggedParquetFiles.quarantined_files + " file(s), including any truncation-flagged samples) " +
                    "has been moved to " + defined_parquet_output_dir + "/" + parquet_quarantine_subdir +
                    "/ and is exempt from both the bulk delete and the bucket's 14-day lifecycle rule; see the README.txt " +
                    "there. Review those samples, then re-run with parquet_fail_on_quarantine = false to finish, or with " +
                    "parquet_allow_flagged_vet_loads = true to waive the screens. THE DATA IS LOADED -- this failure is a " +
                    "notification, not an incomplete ingest.",
          basic_docker = effective_basic_docker,
      }
    }
  }

  # Merge the loaded header scratch data into vcf_header_lines / sample_vcf_header. Gate on whichever
  # load path ran: the Parquet header load (VerifyParquetLoading) or the BQ Write API.
  # Declared before SetIsLoadedColumn because that task may depend on this one's `done` (see below).
  if (load_vcf_headers) {
    call ProcessVCFHeaders {
      input:
        variants_docker = effective_variants_docker,
        go = flatten([
          select_all([VerifyParquetLoading.done]),
          select_all(LoadDataViaBigQueryWriteAPI.done)
        ]),
        dataset_name = dataset_name,
        project_id = project_id,
    }
  }

  if (load_vet_and_ref_ranges) {
    call SetIsLoadedColumn {
      input:
        # A BQ Write API-flavored invocation of `LoadData` actually loads all data into vet and ref ranges tables, but a
        # Parquet-flavored invocation of `LoadData` only generates Parquet files from input gVCFs.
        # Because the loading of Parquet data into BigQuery is handled by a chain of WDL tasks subsequent to
        # `GenerateParquetFilesFromInputGVCFs`, the `go` trigger for setting the `is_loaded` column is the `done` output
        # of the last task in that chain, `VerifyParquetLoading`. The other component of the `go` trigger is the
        # `LoadDataViaBigQueryWriteAPI.done` corresponding to the Write API flow.
        # Intentionally using select_first to pick whichever of the two mutually exclusive code paths (Parquet vs WriteAPI) ran.
        #
        # Also gate on ProcessVCFHeaders.done: this task computes is_loaded from `samples_with_all_data`, which
        # JOINs `samples_with_header_data` (backed by `sample_vcf_header` on the Parquet path). That table is
        # populated by ProcessVCFHeaders, and the Parquet path writes no HEADERS_LOADED status to rescue the
        # view. When a single run sets both load_vcf_headers and load_vet_and_ref_ranges, without this edge the
        # two tasks race and is_loaded is silently left FALSE. select_all makes it a no-op for data-only runs
        # (headers loaded in a prior run, so sample_vcf_header is already populated).
        #@ except: UnnecessaryFunctionCall
        go = flatten([
          select_all(select_first([[VerifyParquetLoading.done], LoadDataViaBigQueryWriteAPI.done])),
          select_all([ProcessVCFHeaders.done])
        ]),
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
    Boolean? parquet_safe_to_delete = VerifyParquetLoading.safe_to_delete_parquet
    Int? parquet_files_loaded = VerifyParquetLoading.loaded_files
    Int? parquet_total_files = VerifyParquetLoading.total_files
    # Independent structural-check observability (VS-1989). parquet_loading_verified (all_loaded) is the
    # factual "load complete?" verdict; parquet_safe_to_delete (safe_to_delete_parquet) is the bulk
    # deletion gate. A vet screen flag no longer separates them: the flagged samples' Parquet is moved
    # under the quarantine subdirectory and the rest of the callset is deleted as usual, so a flagged
    # run is both all_loaded and safe_to_delete. They differ only in the one case where a flagged
    # sample had no Parquet path to move, which leaves the delete blocked because running it would
    # destroy exactly the files the screens asked to keep.
    #
    # The screen flags and the quarantine counts are surfaced because a flagged run is now a SUCCESSFUL
    # run that quietly set some Parquet aside -- without these outputs nothing above this workflow would
    # say so. parquet_quarantined_samples being non-zero on a green run is the signal to go and look --
    # green rather than aborted only when parquet_fail_on_quarantine was turned off, or when every
    # flagged sample was flagged by the truncation screen alone.
    #
    # The family-completeness / ploidy-cardinality components are exact checks that fail the fail-loud
    # VerifyParquetLoading task, whose outputs Cromwell then never delocalizes -- so they could only
    # ever be read as true and are not published; the full verdict is written to
    # verification_results.json, copied to a durable diagnostics path on failure.
    Boolean? parquet_vet_duplication_flagged = VerifyParquetLoading.vet_duplication_flagged
    Boolean? parquet_vet_truncation_flagged = VerifyParquetLoading.vet_truncation_flagged
    Int? parquet_quarantined_samples = VerifyParquetLoading.quarantined_samples
    Int? parquet_quarantined_duplication_samples = VerifyParquetLoading.quarantined_duplication_samples
    Int? parquet_quarantined_files = VerifyParquetLoading.quarantined_files
    File? parquet_quarantine_files_list = VerifyParquetLoading.quarantine_files_list
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

        if [[ "~{load_vet_and_ref_ranges}" = 'true' ]]
        then
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
        fi

        if [[ "~{load_vcf_headers}" = 'true' ]]
        then
          # Header parquet is named vcf_header_lines_scratch_<sampleId>.parquet so DiscoverParquetFiles
          # groups it under the vcf_header_lines_scratch table.
          header_parquet_file=`ls vcf_header_lines_scratch_*.parquet`
          gcloud storage ~{"--billing-project " + billing_project_id} cp $header_parquet_file ${OUTPUT_GCS_DIR}/vcf_header_lines_scratch/$header_parquet_file
        fi

        # cleanup after ourselves
        rm -f *.parquet
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

  output {
    # Exists so downstream tasks (e.g. SetIsLoadedColumn) can order themselves after the header merge.
    Boolean done = true
  }

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

      -- Probe sample_vcf_header, the table the view below actually reads (it is created alongside
      -- vcf_header_lines_scratch when load_vcf_headers is set), so the guard matches its dependency.
      SET header_table_exists = (
        SELECT COUNT(1) FROM
        `~{project_id}.~{dataset_name}.INFORMATION_SCHEMA.TABLES`
        WHERE table_name = 'sample_vcf_header'
      );

      IF header_table_exists > 0 THEN
        IF sample_load_status_table_exists > 0 THEN
          SET headers_load_status_clause = format(sample_load_status_template, 'HEADERS_LOADED');
        ELSE
          SET headers_load_status_clause = '';
        END IF;

        -- Read the DURABLE sample->header map, not the transient scratch table. `vcf_header_lines_scratch`
        -- is emptied by ProcessVCFHeaders after each merge, so a view over it goes empty once headers are
        -- loaded; combined with the Parquet path not writing HEADERS_LOADED, that left this view empty for
        -- Parquet callsets and prevented SetIsLoadedColumn from ever setting is_loaded. sample_vcf_header
        -- persists and is populated for both the BQ and Parquet paths. See the VS-1968 design doc.
        --
        -- Safe for the BQ path: the HEADERS_LOADED status clause (UNION'd below) still covers the
        -- during-ingest window before the scratch->final merge, and post-merge sample_vcf_header is
        -- populated for BQ as well -- so switching the data-presence source from scratch to
        -- sample_vcf_header does not lose any BQ samples from this view.
        SET create_header_data_view = """

          CREATE OR REPLACE VIEW `~{project_id}.~{dataset_name}.samples_with_header_data` AS
          (
          SELECT DISTINCT sample_id FROM `~{project_id}.~{dataset_name}.sample_vcf_header`

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
        "matchesPrefix": ["${BUCKET_PATH_PREFIX}vet/", "${BUCKET_PATH_PREFIX}ref_ranges/", "${BUCKET_PATH_PREFIX}sample_chromosome_ploidy/", "${BUCKET_PATH_PREFIX}vcf_header_lines_scratch/"]
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
    # Subdirectory of output_gcs_dir holding Parquet quarantined by a previous run (VS-1989). Excluded
    # from the listing below: those files are already loaded, and re-discovering them would double-count
    # their samples and hand the loader a second copy of work it has done. The quarantine rename also
    # takes them out of the ".parquet" filter on its own; this exclusion is by location as well, so the
    # property holds even for a file someone has restored in place.
    String quarantine_subdir = "quarantine"
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

    # List all objects, filter for Parquet files.
    #
    # We deliberately tolerate exactly one gcloud failure mode -- stderr containing "One or more
    # URLs matched no objects" -- and treat it as a legitimately empty directory (proceed with an
    # empty file list). Every other non-zero exit fails the task. This couples GVS idempotency to
    # gcloud's error wording; we go into it with open eyes:
    #   - It is fail-safe. If a future cloud-sdk reworded this message, the empty-directory case
    #     would stop matching and fall through to the `else`, failing the task loudly (retried, then
    #     a visible workflow failure) rather than silently loading nothing.
    #   - gcloud is version-pinned in the Variants image, so the wording can only change at a
    #     deliberate cloud-sdk bump -- a reviewed, integration-tested event -- not under a running
    #     pipeline.
    #   - The match is kept specific (not broadened) so a genuine listing error is never mis-read
    #     as "empty" -- that would be the one silent, dangerous direction.
    # A canary asserts this contract against real gcloud + GCS. build_docker.sh runs it inside every
    # freshly built Variants image (so a cloud-sdk bump that reworded this is caught at rebuild time);
    # it can also be run standalone:
    #   scripts/variantstore/scripts/test/gcs_listing_canary/run_gcs_listing_canary.sh
    # If you change the grep pattern below, change the SENTINEL in that script too.
    echo "Listing files in ${OUTPUT_GCS_DIR}..."
    set +o errexit
    gcloud storage ls --recursive ~{"--billing-project " + billing_project_id} \
      "${OUTPUT_GCS_DIR}/" > all_objects.txt 2> gcloud_ls_stderr.txt
    GCLOUD_LS_EXIT_CODE=$?
    set -o errexit

    if [[ ${GCLOUD_LS_EXIT_CODE} -ne 0 ]]; then
      if grep -q 'One or more URLs matched no objects' gcloud_ls_stderr.txt; then
        echo "No objects found under ${OUTPUT_GCS_DIR}/, proceeding with an empty file list."
        : > all_objects.txt
      else
        echo "gcloud storage ls failed unexpectedly:" >&2
        cat gcloud_ls_stderr.txt >&2
        exit "${GCLOUD_LS_EXIT_CODE}"
      fi
    fi

    grep '\.parquet$' all_objects.txt > all_parquet_objects.txt || touch all_parquet_objects.txt

    # Drop anything a previous run quarantined (VS-1989). Those files were already loaded -- they were
    # set aside for inspection, not left unloaded -- so discovering them again would add a duplicate
    # (table, sample_id) entry for every quarantined sample, inflating the file counts the verifier
    # reconciles and handing LoadParquetFilesToBQ a second copy of the same work.
    # awk index() rather than grep -v "^..." so the prefix is matched literally: a GCS bucket name may
    # contain dots, which as a regex would match any character and could over-exclude.
    awk -v pfx="${OUTPUT_GCS_DIR}/~{quarantine_subdir}/" 'index($0, pfx) != 1' \
      all_parquet_objects.txt > all_files.txt
    QUARANTINED_COUNT=$(( $(wc -l < all_parquet_objects.txt) - $(wc -l < all_files.txt) ))
    if [[ ${QUARANTINED_COUNT} -gt 0 ]]; then
      echo "Skipping ${QUARANTINED_COUNT} previously quarantined Parquet file(s) under ${OUTPUT_GCS_DIR}/~{quarantine_subdir}/"
    fi

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
    maxRetries: 3
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
    Array[String] regular_table_prefixes = ["sample_chromosome_ploidy"]
    Array[String] superpartitioned_table_prefixes = ["vet", "ref_ranges"]
    # Independent structural checks (VS-1989): the two screen ratios and whether a screen flag is
    # waived. By default (false) a flag quarantines the flagged samples' Parquet; true waives the
    # screens. vet_truncation_threshold is separate from vet_duplication_threshold -- only the high
    # side is calibrated -- and 0 disables the truncation screen alone.
    Float vet_duplication_threshold = 1.6
    Float vet_truncation_threshold = 1.6
    Boolean allow_flagged_vet_loads = false
    # Exact per-sample ploidy row count to validate against (e.g. 24 for WGS); unset infers the mode.
    Int? expected_ploidy_rows_per_sample
    # Optional durable location for the verdict JSON. This task is fail-loud -- a bad load exits non-zero
    # and aborts the workflow -- and Cromwell does not delocalize a failed task's outputs, so copying the
    # results JSON here keeps the diagnostic recoverable at a predictable path on failure. Unset -> no copy
    # (the JSON still lives only in the task execution directory).
    String? verification_diagnostics_gcs_dir
    # Billing project for the diagnostics copy, required when the target bucket is requester-pays.
    String? billing_project_id
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

    # This verification is fail-loud: verify_all_loaded.py exits non-zero when the load is incomplete, and
    # a non-zero task aborts the workflow before any irreversible Parquet deletion. Capture that exit code
    # (rather than letting errexit abort here) so we can first copy the verdict JSON somewhere durable --
    # Cromwell never delocalizes a failed task's outputs -- and then re-raise the code below.
    rc=0
    python3 /app/verify_all_loaded.py \
      --project-id ~{project_id} \
      --dataset-name ~{dataset_name} \
      --gcs-files-list ~{gcs_files_list} \
      --regular-table-prefixes ~{sep=" " regular_table_prefixes} \
      --superpartitioned-table-prefixes ~{sep=" " superpartitioned_table_prefixes} \
      --vet-duplication-threshold ~{vet_duplication_threshold} \
      --vet-truncation-threshold ~{vet_truncation_threshold} \
      ~{true="--allow-flagged-vet-loads" false="" allow_flagged_vet_loads} \
      ~{"--expected-ploidy-rows-per-sample " + expected_ploidy_rows_per_sample} \
      --output-dir verification_output || rc=$?

    # Copy the verdict JSON to a durable location if one was configured, so the diagnostic survives the
    # fail-loud abort above. The copy must only happen on failure ($rc -ne 0) so a successful retry does
    # not overwrite the failure diagnostic that this path is meant to preserve. The copy must never mask
    # the verification verdict, hence the trailing `|| true`.
    diagnostics_dir='~{default="" verification_diagnostics_gcs_dir}'
    if [[ $rc -ne 0 && -n "${diagnostics_dir}" && -f verification_output/verification_results.json ]]
    then
      gcloud storage cp ~{"--billing-project " + billing_project_id} \
        verification_output/verification_results.json "${diagnostics_dir%/}/verification_results.json" || true
    fi

    # A quarantine on an otherwise-green run needs the same durable evidence: the task succeeds, so the
    # only record of why those files were set aside would otherwise be this execution directory, which
    # is exactly what gets cleaned up. Written under a distinct name so this path -- which a successful
    # retry does re-run -- cannot overwrite a preserved failure diagnostic above.
    if [[ $rc -eq 0 && -n "${diagnostics_dir}" && -s verification_output/quarantine_files.txt ]]
    then
      gcloud storage cp ~{"--billing-project " + billing_project_id} \
        verification_output/verification_results.json \
        "${diagnostics_dir%/}/verification_results.quarantine.json" || true
      gcloud storage cp ~{"--billing-project " + billing_project_id} \
        verification_output/quarantine_files.txt "${diagnostics_dir%/}/quarantine_files.txt" || true
    fi

    exit $rc
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
    # The bulk deletion gate DeleteParquetFiles is conditioned on: all_loaded, plus every sample the
    # heuristic vet screens flagged having Parquet the quarantine step can move aside. A screen flag by
    # itself no longer withholds the delete -- the flagged samples' files are quarantined and the rest
    # of the callset is deleted -- so on a flagged-but-complete load this reads true alongside
    # all_loaded. It is still distinct from all_loaded because a flagged sample with no resolvable
    # Parquet path leaves nothing to quarantine, and deleting then would destroy the very files the
    # screens asked to keep. Meaningful precisely when the task succeeds, so -- unlike the exact-check
    # components below -- it is safe to publish as a task output.
    Boolean safe_to_delete_parquet = read_json(results_json)["safe_to_delete_parquet"]
    Int total_files = read_json(results_json)["total_files"]
    Int loaded_files = read_json(results_json)["loaded_files"]
    Int missing_files = read_json(results_json)["missing_files"]
    File? missing_files_list = "verification_output/missing_files.txt"
    # Independent structural-check observability (VS-1989), read shallowly. Both vet screens --
    # duplication and truncation -- are exposed as task outputs: they never fail all_loaded, so they
    # are meaningful precisely when the task succeeds. Since a flag no longer blocks deletion, these
    # and the quarantine counts are the only thing that says a green run set some Parquet aside.
    #
    # quarantine_files_list is the work list QuarantineFlaggedParquetFiles consumes, and is written
    # unconditionally (empty on a clean run) so that task can be unconditional rather than gated on an
    # optional output. It is derived from the uncapped screen detail, not from the capped per-sample
    # lists in results_json, so it cannot silently omit a flagged file.
    #
    # The family-completeness / ploidy-cardinality components are exact checks -- if any fails,
    # verify_all_loaded.py exits non-zero and the task fails, at which point Cromwell does not evaluate
    # outputs at all. Publishing them as task outputs could only ever yield true, so they are omitted;
    # the complete verdict lives in results_json (copied to a durable diagnostics path on failure).
    Boolean vet_duplication_flagged = read_json(results_json)["vet_duplication_flagged"]
    Boolean vet_truncation_flagged = read_json(results_json)["vet_truncation_flagged"]
    File quarantine_files_list = "verification_output/quarantine_files.txt"
    Int quarantined_files = read_json(results_json)["quarantine_files"]
    Int quarantined_samples = read_json(results_json)["quarantine_sample_count"]
    # Of those samples, the ones the duplication screen objected to. The workflow's abort condition
    # keys on this rather than on quarantined_samples because the duplication threshold is the only
    # calibrated one -- Foxtrot flagged 2 samples of 540,545 -- while the truncation threshold has been
    # measured on its high side only, and a heuristic that has never had its false-positive rate
    # measured should not be able to fail a completed 500k-sample ingest.
    Int quarantined_duplication_samples = read_json(results_json)["quarantine_duplication_sample_count"]
    Boolean done = true
  }
}

task QuarantineFlaggedParquetFiles {
  meta {
    description: "Move the Parquet of samples the heuristic vet screens flagged out of the bulk delete's reach, so it survives for inspection and re-ingest (VS-1989)."
    volatile: true
  }

  input {
    String output_gcs_dir
    String quarantine_subdir
    # Appended to each quarantined object's name so that it no longer ends in ".parquet". Two things
    # then have to go wrong at once for a quarantined file to be lost, rather than one: DeleteParquetFiles'
    # default strategy globs "*.parquet" across the whole output dir and so is defeated by the rename,
    # while its alternate strategy deletes only the four table directories and so is defeated by the
    # location. Same for the bucket's 14-day lifecycle rule, which matches those four prefixes.
    String quarantine_suffix
    # Work list from VerifyParquetLoading, one gs:// path per line, empty on a clean run.
    File quarantine_files_list
    # A flag rate this high is a miscalibrated screen, not a callset with this many bad samples, so
    # refuse rather than grind through a per-object move. Foxtrot's calibration flagged 2 samples out
    # of 540,545 on vet; the same threshold applied to ref_ranges flagged 38,629, which is the shape of
    # mistake this guards against.
    Int max_quarantine_files = 10000

    String? billing_project_id
    String cloud_sdk_docker
  }

  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    # Normalize GCS path by removing any trailing slash
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')
    QUARANTINE_DIR="${OUTPUT_GCS_DIR}/~{quarantine_subdir}"

    FILE_COUNT=$(grep -c . ~{quarantine_files_list} || true)

    # Published as a task output so the workflow can abort a green-but-quarantined run off this count.
    # Written here, before the early exit below, so the clean-run path reports 0 rather than no output.
    echo "${FILE_COUNT}" > quarantine_count.txt

    if [[ "${FILE_COUNT}" -eq 0 ]]; then
      echo "No samples were flagged by the duplication or truncation screens; nothing to quarantine."
      exit 0
    fi

    if [[ "${FILE_COUNT}" -gt ~{max_quarantine_files} ]]; then
      echo "ERROR: ${FILE_COUNT} Parquet files were flagged for quarantine, over the limit of ~{max_quarantine_files}." >&2
      echo "A flag rate this high means the screens are miscalibrated for this callset, not that this" >&2
      echo "many samples are genuinely bad. Refusing to move them object by object." >&2
      echo "Review verification_results.json, then either retune the screen that fired --" >&2
      echo "parquet_vet_duplication_threshold, or parquet_vet_truncation_threshold (0 disables it) --" >&2
      echo "or set parquet_allow_flagged_vet_loads to waive both screens for this run." >&2
      exit 1
    fi

    echo "Quarantining ${FILE_COUNT} Parquet file(s) to ${QUARANTINE_DIR}/"

    # The move preserves each file's path relative to the output dir, so vet/001/<file> lands at
    # quarantine/vet/001/<file>. That keeps which table and sample a quarantined file belongs to
    # readable from its path, and makes a collision between two families' files impossible.
    while IFS= read -r src
    do
      if [[ -z "${src}" ]]
      then
        continue
      fi
      rel="${src#"${OUTPUT_GCS_DIR}/"}"
      if [[ "${rel}" == "${src}" ]]
      then
        # Nothing stripped, so this path is not under the output dir. Refuse rather than invent a
        # destination for it: an unexpected path here means the work list and this task disagree about
        # the layout, and guessing would move a file somewhere nobody will look for it.
        echo "ERROR: ${src} is not under ${OUTPUT_GCS_DIR}/; refusing to move it." >&2
        exit 1
      fi
      gcloud storage mv ~{"--billing-project " + billing_project_id} \
        "${src}" "${QUARANTINE_DIR}/${rel}~{quarantine_suffix}"
    done < ~{quarantine_files_list}

    # Leave the restore instructions where whoever finds the quarantine will find them too. The .txt
    # extension keeps this file out of every Parquet glob in this workflow.
    {
      echo "Parquet quarantined by GvsImportGenomes because the VS-1989 duplication or truncation"
      echo "screen flagged the owning sample. These files were loaded to BigQuery; they were held back"
      echo "from deletion so the load can be inspected, not because the load failed."
      echo
      echo "Each object keeps its path relative to ${OUTPUT_GCS_DIR}/, with '~{quarantine_suffix}'"
      echo "appended to its name. To restore one for re-ingest, copy it back and drop that suffix:"
      echo
      echo "  gcloud storage cp <object> ${OUTPUT_GCS_DIR}/<relative path without the suffix>"
      echo
      echo "The suffix is what keeps these objects out of the bulk delete's '*.parquet' glob, so do not"
      echo "strip it in place."
      echo
      echo "Nothing deletes this directory: it is outside the four prefixes the bucket's 14-day"
      echo "lifecycle rule matches. Clean it up by hand once the samples have been reviewed."
      echo
      echo "Quarantined files:"
      cat ~{quarantine_files_list}
    } > quarantine_README.txt
    gcloud storage cp ~{"--billing-project " + billing_project_id} \
      quarantine_README.txt "${QUARANTINE_DIR}/README.txt"

    echo "✓ Quarantined ${FILE_COUNT} Parquet file(s) under ${QUARANTINE_DIR}/."
    echo "These are exempt from both the bulk delete and the bucket's 14-day lifecycle rule, so they"
    echo "will persist until someone removes them. Review them and clean up when done."
  >>>

  output {
    Boolean done = true
    # How many files this task actually moved. The workflow's abort condition reads this rather than
    # VerifyParquetLoading's equivalent count, because the data dependency is what sequences the abort
    # after the move: terminating off the verification output instead would let Cromwell kill the
    # workflow while the flagged files are still under the prefixes the 14-day lifecycle rule matches,
    # neither quarantined nor deleted.
    Int quarantined_files = read_int("quarantine_count.txt")
  }

  runtime {
    docker: cloud_sdk_docker
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    # Deliberately not preemptible and not retried. A partial move leaves some sources already gone, so
    # a second attempt would fail on the missing source and fail permanently anyway; failing on the
    # first attempt keeps the failure legible, and DeleteParquetFiles consumes this task's `done`, so a
    # failure here holds the delete back rather than letting it run over an incomplete quarantine.
    preemptible: 0
    cpu: 1
  }
}

task DeleteParquetFiles {
  input {
    String output_gcs_dir
    # Suffix appended to every quarantined object's name (VS-1989). Passed in so the assertion below
    # can state, and check, this task's dependency on it: the default strategy's glob deletes
    # "*.parquet" across the whole output dir, and a quarantined object survives that glob precisely
    # because its name no longer ends in ".parquet".
    String quarantine_suffix
    Boolean use_alternate_delete_strategy = false

    String? billing_project_id
    # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
    # is passed here to prevent this task from running until the upstream task has completed.
    #@ except: UnusedInput
    Boolean go
    String cloud_sdk_docker
  }

  command <<<
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o xtrace -o pipefail

    # Normalize GCS path by removing any trailing slash
    OUTPUT_GCS_DIR=$(echo ~{output_gcs_dir} | sed 's/\/$//')

    # Neither strategy below may touch quarantined Parquet (VS-1989). Two independent things keep it
    # safe, and this task depends on both:
    #
    #  * The default strategy globs "*.parquet" across the whole output dir, so it is the quarantine
    #    SUFFIX that saves those objects -- QuarantineFlaggedParquetFiles renames each one to end in
    #    the suffix below, which the glob no longer matches.
    #  * The alternate strategy deletes only the four named table directories, so it is the quarantine
    #    LOCATION that saves them -- they live in a sibling subdirectory that is not in that list.
    #
    # Assert the first of those, since it is a property of a string and so can be checked here. A
    # suffix ending in ".parquet" would make the default strategy delete the whole quarantine.
    case "~{quarantine_suffix}" in
      *.parquet)
        echo "ERROR: quarantine suffix '~{quarantine_suffix}' ends in .parquet, so the deletion glob" >&2
        echo "below would match quarantined objects and destroy them. Choose another suffix." >&2
        exit 1
        ;;
      "")
        echo "ERROR: quarantine suffix is empty, so quarantined objects keep their .parquet names and" >&2
        echo "the deletion glob below would match them. Choose a non-empty suffix." >&2
        exit 1
        ;;
    esac

    if [ "~{use_alternate_delete_strategy}" = "false" ]; then
      gcloud storage rm --recursive ~{"--billing-project " + billing_project_id} "${OUTPUT_GCS_DIR}/"'**/*.parquet'
    else
      # List the contents of the vet and ref_ranges directories for subsequent deletion in the loop below
      echo "Listing directories under ${OUTPUT_GCS_DIR}/vet/ and ${OUTPUT_GCS_DIR}/ref_ranges/ ${OUTPUT_GCS_DIR}/sample_chromosome_ploidy/ for deletion..."
      gcloud storage ls ~{"--billing-project " + billing_project_id} \
        "${OUTPUT_GCS_DIR}/vet/" "${OUTPUT_GCS_DIR}/ref_ranges/" > parquet_dirs.txt
      echo "${OUTPUT_GCS_DIR}/sample_chromosome_ploidy/" >> parquet_dirs.txt
      # Only present when headers were loaded; guard on existence so data-only runs don't fail here.
      if gcloud storage ls ~{"--billing-project " + billing_project_id} "${OUTPUT_GCS_DIR}/vcf_header_lines_scratch/" >/dev/null 2>&1; then
        echo "${OUTPUT_GCS_DIR}/vcf_header_lines_scratch/" >> parquet_dirs.txt
      fi

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
