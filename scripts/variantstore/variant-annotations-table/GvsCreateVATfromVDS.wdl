version 1.0

import "../wdl/GvsUtils.wdl" as Utils
import "GvsCreateVATFilesFromBigQuery.wdl" as GvsCreateVATFilesFromBigQuery

workflow GvsCreateVATfromVDS {
    input {
        String project_id
        String dataset_name
        File ancestry_file
        String filter_set_name
        File? sites_only_vcf
        File? sites_only_vcf_index
        File? vds_path
        File? sites_to_exclude
        String output_path

        String? basic_docker
        String? git_branch_or_tag
        String? hail_version
        File? hail_wheel
        String? vat_version
        String? workspace_project

        Boolean generate_vep_and_loftee_annotations = true
        Boolean leave_hail_cluster_running_at_end = false
        Boolean use_tiny_dataproc_cluster = false
        Boolean use_tiny_vep_annotation_load_runtime = false
        Int? merge_vcfs_disk_size_override
        Int? split_intervals_disk_size_override
        Int? split_intervals_mem_override
        Int? split_intervals_scatter_count
        Boolean use_reference_disk = true

        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? gatk_docker
        String? variants_docker
        String? variants_nirvana_docker
        String? vep_loftee_docker

        String? vep_loftee_data_table_raw
        String? vep_loftee_data_table_cooked

        String loftee_references_dir = "gs://gvs-internal/loftee/"
    }

    parameter_meta {
        project_id: {
            help: "Google project ID for the GVS BigQuery dataset"
        }
        dataset_name: {
            help: "BigQuery dataset name for GVS"
        }
        ancestry_file: {
            help: "TSV file in GCS where the first column is the research ID and the last is the derived ancestry"
        }
        filter_set_name: {
            help: "name of the filter set used to generate the callset in GVS"
        }
        sites_only_vcf: {
            help: "Optional sites-only VCF file. If defined, generation of a sites-only VCF from the VDS will be skipped. If defined, then 'vds_path' must NOT be defined."
        }
        sites_only_vcf_index: {
            help: "Index to accompany sites-only VCF file documented above."
        }
        output_path: {
            help: "GCS location (with a trailing '/') to put temporary and output files for the VAT pipeline"
        }
        vds_path: {
            help: "Optional top-level directory of the GVS VDS to be used to create the VAT. If defined, then 'sites_only_vcf' must NOT be defined"
        }
        sites_to_exclude: {
            help: "An optional file of sites to exclude from the sites-only VCF. It may become necessary to specify this if annotations for a particular position have issues that prevent Nirvana from running successfully, e.g. chr2:20447683 observed in AnVIL 3K data. The format is one bcftools-style region per line, e.g. 'chr2:20447683', no header."
        }
        generate_vep_and_loftee_annotations: {
            help: "If true (the default), run VEP + LOFTEE and load the resulting annotations into the VAT. Set to false to skip VEP + LOFTEE annotation entirely, which can significantly reduce runtime and cost for large callsets (500K+ samples). The VAT BigQuery schema is unchanged regardless of this setting; VEP/LOFTEE columns will simply be null when disabled."
        }
    }

    String region = "us-central1"
    File mane_annotation_file = "gs://gvs_quickstart_storage/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.summary.txt"

    # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
    # no calling WDLs that might supply `git_hash`).
    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_variants_nirvana_docker = select_first([variants_nirvana_docker, GetToolVersions.variants_nirvana_docker])
    String effective_vep_loftee_docker = select_first([vep_loftee_docker, GetToolVersions.vep_loftee_docker])
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])
    String effective_google_project = select_first([workspace_project, GetToolVersions.google_project])

    # If the vat version is undefined or v1 then the vat tables would be named like filter_vat, otherwise filter_vat_v2.
    String effective_vat_version = if (defined(vat_version) && select_first([vat_version]) != "v1") then "_" + select_first([vat_version]) else ""
    String effective_vat_table_name = filter_set_name + "_vat" + effective_vat_version

    Int vep_annotation_load_cpu         = if use_tiny_vep_annotation_load_runtime then 2                     else 16
    String vep_annotation_load_memory   = if use_tiny_vep_annotation_load_runtime then "7 GB"                else "16 GB"
    String vep_annotation_load_disks    = if use_tiny_vep_annotation_load_runtime then "local-disk 1000 HDD" else "local-disk 4000 SSD"
    Int vep_annotation_load_preemptible = if use_tiny_vep_annotation_load_runtime then 2                     else 0

    String output_path_without_a_trailing_slash = sub(output_path, "/$", "")
    String effective_output_path = if (output_path == output_path_without_a_trailing_slash) then output_path + "/" else output_path
    String variant_transcripts_output_path = effective_output_path + 'variant_transcripts/'
    String genes_output_path = effective_output_path + 'genes/'

    if ((defined(sites_only_vcf)) && (defined(vds_path))) {
        call Utils.TerminateWorkflow as IfSitesOnlyVcfSetDontSetCreateParameters {
            input:
                message = "Error: If 'sites_only_vcf' is set as an input, you may not set 'vds_path'",
                basic_docker = effective_basic_docker,
        }
    }

    if (!defined(sites_only_vcf) && (!defined(vds_path))) {
        call Utils.TerminateWorkflow as MustSetSitesOnlyVcfCreationParameters {
            input:
                message = "Error: If 'sites_only_vcf' is not set as an input, you MUST set 'vds_path'",
                basic_docker = effective_basic_docker,
        }
    }

    if (defined(sites_only_vcf_index) && !defined(sites_only_vcf)) {
        call Utils.TerminateWorkflow as IfSitesOnlyVcfIndexSpecifiedMustSpecifySitesOnlyVcf {
            input:
                message = "Error: If 'sites_only_vcf_index' is set as an input, you must also set 'sites_only_vcf'",
                basic_docker = effective_basic_docker,
        }
    }

    call Utils.GetHailScripts {
        input:
            variants_docker = effective_variants_docker,
    }

    call Utils.GetReference {
        input:
            basic_docker = effective_basic_docker,
    }

    String interval_list = GetReference.reference.wgs_calling_interval_list
    String reference_fasta = GetReference.reference.reference_fasta

    if (defined(sites_only_vcf) || (defined(vds_path))) {
        if (!defined(split_intervals_scatter_count)) {
            call Utils.GetBQTableLastModifiedDatetime as SampleDateTime {
                input:
                    project_id = project_id,
                    fq_table = "~{project_id}.~{dataset_name}.sample_info",
                    cloud_sdk_docker = effective_cloud_sdk_docker,
            }

            call Utils.GetNumSamplesLoaded {
                input:
                    fq_sample_table ="~{project_id}.~{dataset_name}.sample_info",
                    project_id = project_id,
                    sample_table_timestamp = SampleDateTime.last_modified_timestamp,
                    cloud_sdk_docker = effective_cloud_sdk_docker,
            }

            Int calculated_scatter_count = if (GetNumSamplesLoaded.num_samples < 11) then 10 else
                                               if (GetNumSamplesLoaded.num_samples < 250000) then 500 else
                                                   if (GetNumSamplesLoaded.num_samples < 450000) then 1000 else 2000
        }

        Int effective_scatter_count = select_first([split_intervals_scatter_count, calculated_scatter_count])

        call LoadManeDataIntoBigQuery {
            input:
                project_id = project_id,
                dataset_name = dataset_name,
                mane_table_name = "mane_annotations",
                mane_data_file = mane_annotation_file,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call MakeSubpopulationFilesAndReadSchemaFiles {
            input:
                input_ancestry_file = ancestry_file,
                variants_docker = effective_variants_docker,
        }

        if (!defined(sites_only_vcf)) {
            call GenerateSitesOnlyVcf {
                input:
                    use_tiny_dataproc_cluster = use_tiny_dataproc_cluster,
                    vds_path = select_first([vds_path]),
                    workspace_project = effective_google_project,
                    hail_version = effective_hail_version,
                    hail_wheel = hail_wheel,
                    ancestry_file_path = MakeSubpopulationFilesAndReadSchemaFiles.ancestry_file_path,
                    workspace_bucket = GetToolVersions.workspace_bucket,
                    region = region,
                    leave_cluster_running_at_end = leave_hail_cluster_running_at_end,
                    cloud_sdk_docker = effective_cloud_sdk_docker,
                    run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
                    create_vat_inputs_script = GetHailScripts.create_vat_inputs_script,
                    hail_create_vat_inputs_script = GetHailScripts.hail_create_vat_inputs_script,
            }
        }

        if (defined(sites_to_exclude)) {
            call ExcludeSitesFromSitesOnlyVcf {
                input:
                    sites_to_exclude = select_first([sites_to_exclude]),
                    input_sites_only_vcf = select_first([sites_only_vcf, GenerateSitesOnlyVcf.sites_only_vcf]),
                    variants_docker = effective_variants_docker,
            }
        }

        call Utils.CopyFile as CopySitesOnlyVcf {
            input:
                input_file = select_first([ExcludeSitesFromSitesOnlyVcf.output_sites_only_vcf, sites_only_vcf, GenerateSitesOnlyVcf.sites_only_vcf]),
                output_gcs_dir = effective_output_path + "sites_only_vcf",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        if (!defined(ExcludeSitesFromSitesOnlyVcf.output_sites_only_vcf_index) && !defined(sites_only_vcf_index)) {
            call Utils.IndexVcf {
                input:
                    input_vcf = CopySitesOnlyVcf.output_file_path,
                    gatk_docker = effective_gatk_docker,
            }
        }

        call Utils.CopyFile as CopySitesOnlyVcfIndex {
            input:
                input_file = select_first([ExcludeSitesFromSitesOnlyVcf.output_sites_only_vcf_index, sites_only_vcf_index, IndexVcf.output_vcf_index]),
                output_gcs_dir = effective_output_path + "sites_only_vcf",
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call Utils.SplitIntervals {
            input:
                intervals = interval_list,
                ref_fasta = reference_fasta,
                scatter_count = effective_scatter_count,
                output_gcs_dir = effective_output_path + "intervals",
                split_intervals_disk_size_override = split_intervals_disk_size_override,
                split_intervals_mem_override = split_intervals_mem_override,
                gatk_docker = effective_gatk_docker,
        }

        String sites_only_vcf_basename = basename(CopySitesOnlyVcf.output_file_path, ".sites-only.vcf.bgz")

        scatter(i in range(length(SplitIntervals.interval_files))) {
            String interval_file_basename = basename(SplitIntervals.interval_files[i], ".interval_list")
            String vcf_filename = interval_file_basename + "." + sites_only_vcf_basename

            call Utils.SelectVariants {
                input:
                    input_vcf = CopySitesOnlyVcf.output_file_path,
                    input_vcf_index = CopySitesOnlyVcfIndex.output_file_path,
                    interval_list = SplitIntervals.interval_files[i],
                    output_basename = vcf_filename,
                    gatk_docker = effective_gatk_docker,
            }

            call RemoveDuplicatesFromSitesOnlyVCF {
                input:
                    sites_only_vcf = SelectVariants.output_vcf,
                    ref = reference_fasta,
                    variants_docker = effective_variants_docker,
            }

            call StripCustomAnnotationsFromSitesOnlyVCF {
                input:
                    input_vcf = RemoveDuplicatesFromSitesOnlyVCF.output_vcf,
                    custom_annotations_header = MakeSubpopulationFilesAndReadSchemaFiles.custom_annotations_template_file,
                    output_vcf_name = "${vcf_filename}.unannotated.sites_only.vcf",
                    output_custom_annotations_filename = "${vcf_filename}.custom_annotations.tsv",
                    variants_docker = effective_variants_docker,
            }

            if (generate_vep_and_loftee_annotations) {
                call GenerateVepAndLofteeAnnotations {
                    input:
                        vep_loftee_docker = effective_vep_loftee_docker,
                        vep_cache = loftee_references_dir + "homo_sapiens_vep_115_GRCh38.tar.gz",
                        loftee_human_ancestor_fa_gz = loftee_references_dir + "human_ancestor.fa.gz",
                        loftee_human_ancestor_fa_gz_fai = loftee_references_dir + "human_ancestor.fa.gz.fai",
                        loftee_human_ancestor_fa_gz_gzi = loftee_references_dir + "human_ancestor.fa.gz.gzi",
                        loftee_gerp_scores = loftee_references_dir + "gerp_conservation_scores.homo_sapiens.GRCh38.bw",
                        loftee_phylo_csf_database = loftee_references_dir + "loftee.sql",
                        input_vcf = StripCustomAnnotationsFromSitesOnlyVCF.output_vcf,
                }
            }

            ## Use Nirvana to annotate the sites-only VCF and include the AC/AN/AF calculations as custom annotations
            call AnnotateVCF {
                input:
                    cromwell_root = GetToolVersions.cromwell_root,
                    input_vcf = StripCustomAnnotationsFromSitesOnlyVCF.output_vcf,
                    output_annotated_file_name = "${vcf_filename}_annotated",
                    custom_annotations_file = StripCustomAnnotationsFromSitesOnlyVCF.output_custom_annotations_file,
                    variants_nirvana_docker = effective_variants_nirvana_docker,
                    use_reference_disk = use_reference_disk,
            }

            call PrepVtAnnotationJson {
                input:
                    positions_annotation_json = AnnotateVCF.positions_annotation_json,
                    output_file_suffix = "${vcf_filename}.json.gz",
                    output_path = variant_transcripts_output_path,
                    variants_docker = effective_variants_docker,
            }

            call PrepGenesAnnotationJson {
                input:
                    genes_annotation_json = AnnotateVCF.genes_annotation_json,
                    output_file_suffix = "${vcf_filename}.json.gz",
                    output_path = genes_output_path,
                    variants_docker = effective_variants_docker,
            }
        }

        if (generate_vep_and_loftee_annotations) {
            call BigQueryLoadRawVepAndLofteeAnnotations {
                input:
                    vep_loftee_raw_output = select_all(GenerateVepAndLofteeAnnotations.output_file),
                    project_id = project_id,
                    dataset_name = dataset_name,
                    raw_data_table = select_first([vep_loftee_data_table_raw, "vep_loftee_data_table_raw"]),
                    raw_data_table_schema = MakeSubpopulationFilesAndReadSchemaFiles.vep_loftee_raw_schema_json_file,
                    variants_docker = effective_variants_docker,
                    runtime_cpu = vep_annotation_load_cpu,
                    runtime_memory = vep_annotation_load_memory,
                    runtime_disks = vep_annotation_load_disks,
                    runtime_preemptible = vep_annotation_load_preemptible,
            }

            call BigQueryCookVepAndLofteeRawAnnotations {
                input:
                    go = BigQueryLoadRawVepAndLofteeAnnotations.done,
                    project_id = project_id,
                    dataset_name = dataset_name,
                    raw_data_table = select_first([vep_loftee_data_table_raw, "vep_loftee_data_table_raw"]),
                    cooked_data_table = select_first([vep_loftee_data_table_cooked, "vep_loftee_data_table_cooked"]),
                    cooked_data_table_schema = MakeSubpopulationFilesAndReadSchemaFiles.vep_loftee_cooked_schema_json_file,
                    variants_docker = effective_variants_docker,
            }
        }

        call Utils.MergeTsvs {
            input:
                input_files = RemoveDuplicatesFromSitesOnlyVCF.track_dropped,
                output_file_name = "${sites_only_vcf_basename}.dropped.tsv",
                basic_docker = effective_basic_docker,
        }

        call Utils.MergeTsvs as MergeDroppedSynonyms {
            input:
                input_files = RemoveDuplicatesFromSitesOnlyVCF.filtered_synonyms,
                output_file_name = "${sites_only_vcf_basename}.detailed-dropped.tsv",
                basic_docker = effective_basic_docker,
        }

        call BigQueryLoadJson {
            input:
                vat_schema = MakeSubpopulationFilesAndReadSchemaFiles.vat_schema_json_file,
                variant_transcript_schema = MakeSubpopulationFilesAndReadSchemaFiles.variant_transcript_schema_json_file,
                genes_schema = MakeSubpopulationFilesAndReadSchemaFiles.genes_schema_json_file,
                mane_table_name = LoadManeDataIntoBigQuery.mane_table,
                vep_loftee_cooked_table_name = BigQueryCookVepAndLofteeRawAnnotations.cooked_table_name,
                run_vep_loftee_update = generate_vep_and_loftee_annotations,
                project_id = project_id,
                dataset_name = dataset_name,
                variant_transcripts_path = variant_transcripts_output_path,
                genes_path = genes_output_path,
                base_vat_table_name = effective_vat_table_name,
                go = flatten([PrepVtAnnotationJson.done, PrepGenesAnnotationJson.done]),
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call DeduplicateVatInBigQuery {
            input:
                input_vat_table_name = BigQueryLoadJson.vat_table,
                output_vat_table_name = effective_vat_table_name,
                vat_schema = MakeSubpopulationFilesAndReadSchemaFiles.vat_schema_json_file,
                project_id = project_id,
                dataset_name = dataset_name,
                cloud_sdk_docker = effective_cloud_sdk_docker,
        }

        call GvsCreateVATFilesFromBigQuery.GvsCreateVATFilesFromBigQuery {
            input:
                project_id = project_id,
                git_branch_or_tag = git_branch_or_tag,
                git_hash = GetToolVersions.git_hash,
                dataset_name = dataset_name,
                vat_table_name = DeduplicateVatInBigQuery.vat_table,
                output_path = effective_output_path,
                merge_vcfs_disk_size_override = merge_vcfs_disk_size_override,
                cloud_sdk_docker = effective_cloud_sdk_docker,
                cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
        }
    }

    output {
        String vat_table_name = effective_vat_table_name
        String? cluster_name = GenerateSitesOnlyVcf.cluster_name
        File? dropped_sites_file = MergeTsvs.output_file
        File? detailed_dropped_sites_file = MergeDroppedSynonyms.output_file
        File? final_tsv_file = GvsCreateVATFilesFromBigQuery.final_tsv_file
        String recorded_git_hash = GetToolVersions.git_hash
    }
}


task ExcludeSitesFromSitesOnlyVcf {
    input {
        File sites_to_exclude
        File input_sites_only_vcf
        Int disk_size_gb = ceil(4 * (size(input_sites_only_vcf, "GiB"))) + 500
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        regions_to_exclude=$(cat ~{sites_to_exclude} | tr '\n' ',' | tr -d '[:space:]' | sed 's/,$//g')

        # Exclude sites from the sites-only VCF that are listed in the sites_to_exclude file.
        bcftools view --threads 4 --targets ^${regions_to_exclude} ~{input_sites_only_vcf} -O z -o excluded.vcf.gz

        # Index the resulting VCF.
        bcftools index --tbi --threads 4 excluded.vcf.gz
    >>>

    runtime {
        docker: variants_docker
        memory: "4 GB"
        preemptible: 2
        cpu: 1
        disks: "local-disk ${disk_size_gb} HDD"
    }

    output {
        File output_sites_only_vcf = "excluded.vcf.gz"
        File output_sites_only_vcf_index = "excluded.vcf.gz.tbi"
    }
}


task GenerateSitesOnlyVcf {
    input {
        Boolean use_tiny_dataproc_cluster
        String vds_path
        String workspace_project
        String workspace_bucket
        String region
        Boolean leave_cluster_running_at_end
        String hail_version
        File run_in_hail_cluster_script
        File hail_create_vat_inputs_script
        File create_vat_inputs_script
        File? hail_wheel
        String ancestry_file_path
        String? hail_temp_path
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction

        String cloud_sdk_docker
    }
    String prefix = "sites-only-vcf"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        python3 -m venv ./localvenv
        . ./localvenv/bin/activate

        pip3 install --upgrade pip

        if [[ ! -z "~{hail_wheel}" ]]
        then
            pip3 install ~{hail_wheel}
        else
            pip3 install hail~{'==' + hail_version}
        fi

        pip3 install --upgrade google-cloud-dataproc ijson

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | od -h -An | tr -d '[:space:]')-$(head -c2 < /dev/urandom | od -h -An | tr -d '[:space:]')"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        sites_only_vcf_filename="~{workspace_bucket}/~{prefix}-${hex}.sites-only.vcf.bgz"
        echo ${sites_only_vcf_filename} > sites_only_vcf_filename.txt

        if [[ -z "~{hail_temp_path}" ]]
        then
            hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"
        else
            hail_temp_path="~{hail_temp_path}"
        fi

        # construct a JSON of arguments for python script to be run in the hail cluster
        cat > script-arguments.json <<FIN
        {
            "vds_input_path": "~{vds_path}",
            "temp_path": "${hail_temp_path}",
            "ancestry_input_path": "~{ancestry_file_path}",
            "sites_only_output_path" : "${sites_only_vcf_filename}"
        }
        FIN

        # Run the hail python script to make a sites-only VCF from a VDS
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{hail_create_vat_inputs_script} \
            --secondary-script-path-list ~{create_vat_inputs_script} \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            ~{true='--use-tiny-dataproc-cluster' false='' use_tiny_dataproc_cluster} \
            --region ~{region} \
            --workspace-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end}
    >>>

    runtime {
        memory: "6.5 GB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: 0
        docker: cloud_sdk_docker
        bootDiskSizeGb: 10
    }

    output {
        String cluster_name = read_string("cluster_name.txt")
        File sites_only_vcf = read_string("sites_only_vcf_filename.txt")
    }
}


task MakeSubpopulationFilesAndReadSchemaFiles {
    input {
        File input_ancestry_file

        String schema_filepath = "/data/variant_annotation_table/schema/"
        String vat_schema_json_filename = "vat_schema.json"
        String variant_transcript_schema_json_filename = "variant_transcript_schema.json"
        String genes_schema_json_filename = "genes_schema.json"
        String vep_loftee_115_raw_schema_json_filename = "vep_loftee_115_raw.json"
        String vep_loftee_115_cooked_schema_json_filename = "vep_loftee_115_cooked.json"
        String variants_docker
    }
    String output_ancestry_filename =  "ancestry_mapping.tsv"
    String custom_annotations_template_filename =  "custom_annotations_template.tsv"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        cp ~{schema_filepath}* .

        ## the ancestry file is processed down to a simple mapping from sample to subpopulation
        python3 /app/extract_subpop.py \
            --input_path ~{input_ancestry_file} \
            --output_path ~{output_ancestry_filename} \
            --custom_annotations_template_path ~{custom_annotations_template_filename}
    >>>

    runtime {
        docker: variants_docker
        memory: "1 GB"
        preemptible: 3
        cpu: 1
        disks: "local-disk 100 HDD"
    }

    output {
        File vat_schema_json_file = vat_schema_json_filename
        File variant_transcript_schema_json_file = variant_transcript_schema_json_filename
        File genes_schema_json_file = genes_schema_json_filename
        File vep_loftee_raw_schema_json_file = vep_loftee_115_raw_schema_json_filename
        File vep_loftee_cooked_schema_json_file = vep_loftee_115_cooked_schema_json_filename

        File ancestry_mapping_list = output_ancestry_filename
        File custom_annotations_template_file = custom_annotations_template_filename
        String ancestry_file_path = input_ancestry_file
    }
}


task StripCustomAnnotationsFromSitesOnlyVCF {
    input {
        File input_vcf
        File custom_annotations_header
        String output_vcf_name
        String output_custom_annotations_filename
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    Int disk_size = ceil((size(input_vcf, "GB") + size(custom_annotations_header, "GB")) * 4) + 100

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        python3 /app/strip_custom_annotations_from_sites_only_vcf.py \
        --input_vcf ~{input_vcf} \
        --input_custom_annotations_tsv ~{custom_annotations_header} \
        --output_vcf ~{output_vcf_name} \
        --output_custom_annotations_tsv ~{output_custom_annotations_filename}

    >>>

    runtime {
        docker: variants_docker
        memory: "7 GiB"
        cpu: 2
        preemptible: 3
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File output_vcf = output_vcf_name
        File output_custom_annotations_file = output_custom_annotations_filename
        File monitoring_log = "monitoring.log"
    }
}


task RemoveDuplicatesFromSitesOnlyVCF {
    input {
        File sites_only_vcf
        File ref
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    Int disk_size = ceil(size(sites_only_vcf, "GB") * 5) + 100

    # separate multi-allelic sites into their own lines, remove deletions and filtered sites and make a sites only vcf
    # while extracting and calculating the an/ac/af & sc by subpopulation into a tsv
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        # custom function to prepend the current datetime to an echo statement
        echo_date () { echo "`date "+%Y/%m/%d %H:%M:%S"` $1"; }

        echo_date "VAT: Convert input to BCF format"
        bcftools convert --threads 4 -O b -o sites_only.bcf ~{sites_only_vcf}

        echo_date "VAT: Calculating number of sites with Ns"

        ## track the dropped variants with N's in the reference (Since Nirvana cant handle N as a base, drop them for now)
        # The "2>/dev/null" is to suppress the profuse and useless "pass=0" stderr output from bcftools view.
        bcftools view --threads 4 -i 'REF~"N"' -O u sites_only.bcf 2>/dev/null | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' > track_dropped.tsv

        echo_date "VAT: filter out sites with N's in the reference AND sites with AC=0"
        ## NOTE: Sites that were filtered out because of AC=0 are not recorded in the 'track_dropped.tsv' file, but can be
        ##       determined by examining the sites-only VCF provided to this WDL.
        bcftools view --threads 4 -e 'REF~"N" || AC=0' -O b sites_only.bcf -o filtered_sites_only.bcf 2>/dev/null
        rm sites_only.bcf

        echo_date "VAT: normalize, left align and split multi allelic sites to new lines, remove duplicate lines"
        ## note that normalization may create sites with 100 or more alt alleles
        bcftools norm --threads 4 -m- --check-ref w -f ~{ref} filtered_sites_only.bcf -O b -o normalized.bcf
        rm filtered_sites_only.bcf

        echo_date "VAT: detecting and removing duplicate rows from sites-only VCF"

        ## After normalization, sometimes duplicate variants appear but with different calculations. This is due to
        ## non-left aligned synonyms (and possibly one left-aligned synonym) appearing in the input data, each synonym
        ## having its own calculations.
        ## We now keep the duplicate with the highest AC value for each synonym cluster
        ## to locate the duplicates, we first make a file of just the first 5 columns
        bcftools query normalized.bcf -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' | sort | uniq -d > duplicates.tsv

        echo_date "VAT: done with duplicate detection"
        wc -l duplicates.tsv
        echo_date "VAT: Duplicates may have been found"

        # If there ARE dupes to process
        if [ -s duplicates.tsv ]; then
            echo_date "VAT: Processing duplicates to keep highest AC variant per synonym cluster"

            # Extract all duplicate lines from the VCF (not removing them yet)
            bcftools view --threads 4 normalized.bcf | grep -wFf duplicates.tsv > all_duplicates.vcf

            # Process each synonym cluster to identify minority synonyms to remove
            > minority_synonyms.tsv  # Initialize empty file

            while IFS=$'\t' read -r chrom pos id ref alt; do
                # Find all lines matching this synonym pattern and extract AC values using a python for loop that prints all but the highest AC line
                grep "^${chrom}\t${pos}\t${id}\t${ref}\t${alt}\t" all_duplicates.vcf | \
                python -c '
import sys, re
max_ac = -1
highest_line = ""
for line in sys.stdin:
    match = re.search(r"AC=(\d+)", line.strip())
    if match:
        # match group 1 is the AC number
        ac = int(match.group(1))
        current_line = line.strip()
        if ac > max_ac:
            if highest_line != "":
                print(highest_line)
            max_ac = ac
            highest_line = current_line
        else:
            print(current_line)
' >> temp_minorities.tsv
            done < duplicates.tsv

            # Collect all minority synonyms
            if [ -f temp_minorities.tsv ]; then
                cat temp_minorities.tsv >> minority_synonyms.tsv
                rm temp_minorities.tsv
            fi

            # Remove only the minority synonyms from the VCF
            if [ -s minority_synonyms.tsv ]; then
                echo_date "VAT: Removing minority synonyms"
                bcftools view --threads 4 normalized.bcf | grep -v -wFf minority_synonyms.tsv > deduplicated.vcf

                ## Copy minority synonyms to filtered_synonyms.tsv for output
                cp minority_synonyms.tsv filtered_synonyms.tsv
            else
                echo_date "VAT: No minority synonyms to remove"
                bcftools view --threads 4 normalized.bcf -o deduplicated.vcf
                ## Create empty filtered_synonyms.tsv if no minorities found
                touch filtered_synonyms.tsv
            fi

            rm -f all_duplicates.vcf minority_synonyms.tsv
        else
            # There are no duplicates to process
            echo_date "VAT: No duplicates found"
            bcftools view --threads 4 normalized.bcf -o deduplicated.vcf
            ## Create empty filtered_synonyms.tsv if no duplicates found
            touch filtered_synonyms.tsv
        fi
        rm normalized.bcf

        ## add duplicates to the file that's tracking dropped variants (original functionality)
        cat duplicates.tsv >> track_dropped.tsv
        rm duplicates.tsv ## clean up unneeded file

        echo_date "VAT: finished"
    >>>

    runtime {
        docker: variants_docker
        maxRetries: 3
        memory: "16 GB"
        preemptible: 3
        cpu: 8
        disks: "local-disk " + disk_size + " HDD"
    }

    output {
        File track_dropped = "track_dropped.tsv"
        File output_vcf = "deduplicated.vcf"
        File filtered_synonyms = "filtered_synonyms.tsv"
        File monitoring_log = "monitoring.log"
    }
}

task GenerateVepAndLofteeAnnotations {
    input {
        String vep_loftee_docker
        # TODO make a reference disk for this stuff, some of these references are quite large.
        File vep_cache
        File loftee_human_ancestor_fa_gz
        File loftee_human_ancestor_fa_gz_fai
        File loftee_human_ancestor_fa_gz_gzi
        File loftee_gerp_scores
        File loftee_phylo_csf_database
        File input_vcf
        File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
        # The memory headroom left for other processes including the Batch agent.
        Float overhead_memory_mib = 1.6 * 1024
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "MEM_SIZE is ${MEM_SIZE}"
        echo "MEM_UNIT is ${MEM_UNIT}"

        if [[ -z "${MEM_UNIT:-}" ]]
        then
            echo "MEM_UNIT environment variable unexpectedly not set." 1>&2
            exit 1
        elif [[ ${MEM_UNIT} == "GB" ]]
        then
            vep_memory_kib=$(python -c "from math import floor; print(int(floor(((${MEM_SIZE} * 1024) - ~{overhead_memory_mib}) * 1024)))")
        else
            echo "Unexpected memory unit: ${MEM_UNIT}" 1>&2
            exit 1
        fi

        echo "overhead_memory_mib is ~{overhead_memory_mib}"
        echo "vep_memory_kib is ${vep_memory_kib}"

        bash ~{monitoring_script} > monitoring.log &

        if { grep -E -v '^#' ~{input_vcf} 2>&1 > /dev/null; }
        then
            # Only copy these references if there are actually data lines in the VCF to be processed,
            # Most of the shards in 20/X/Y integration runs don't have any work to do and don't need
            # to localize the references.
            #
            # gcloud storage cp ~{vep_cache} ~{vep_cache} ~{loftee_human_ancestor_fa_gz} ~{loftee_human_ancestor_fa_gz_fai} ~{loftee_human_ancestor_fa_gz_gzi} ~{loftee_gerp_scores} ~{loftee_phylo_csf_database} .
            #
            # TODO yeah that would be nice but here's no gcloud on the VEP + LOFTEE image. These references
            # *really* should be on a reference disk.
            tar xzf ~{vep_cache}

            LOFTEE_PATH=/opt/vep/src/loftee-1.0.4_GRCh38
            args=(

                # Some logging please.
                --verbose
                --warning_file warnings.txt

                # Explicitly turn off forking as LOFTEE might not deal well with that.
                --fork 1

                # Breaks out data into their own columns that otherwise would be nested (semicolon delimited) in the "Extra" column.
                --tab

                # Force writing versions on Ensembl transcripts for VAT compatibility.
                --transcript_version

                # Emit HGNC symbols and IDs.
                --symbol

                # Basic LOFTEE plugin setup
                --plugin LoF,loftee_path:$LOFTEE_PATH,gerp_bigwig:~{loftee_gerp_scores},human_ancestor_fa:~{loftee_human_ancestor_fa_gz},conservation_file:~{loftee_phylo_csf_database},check_complete_cds:false
                --dir_plugins $LOFTEE_PATH

                # Basic VEP cache setup
                --cache
                --offline
                --dir_cache .

                # For GERP (Genomic Evolutionary Rate Profiling) score output.
                --custom file=~{loftee_gerp_scores},short_name=GERP,format=bigwig,num_records=all

                # Input and output files
                --input_file ~{input_vcf}
                --output_file vep_loftee_raw_output.txt
            )

            # Unfortunately we don't seem to be able to limit the amount of memory the Perl process uses. There are
            # ways of limiting memory in Docker containers if one has access to Docker daemon configuration or in the
            # way the `docker run` command is invoked (cgroups, ulimit, etc.), but unfortunately we don't have those
            # options when running in GCP Batch. If we try to do a `ulimit -H -v <value>` when running under GCP Batch
            # it has no effect.
            #
            # Run the vep process with these default memory settings and a very limited number of retries (maybe zero).
            # Some shards will likely fail, which in turn will fail the workflow. After this happens, edit this wdl to
            # increase the `memory` runtime attribute and rerun with call caching enabled. Repeat as necessary.

            vep "${args[@]}"

        else
            echo "No data found for processing in VCF, exit 0."
            touch "vep_loftee_raw_output.txt"
        fi

    >>>

    runtime {
        preemptible: 2
        maxRetries: 2
        noAddress: true
        docker: vep_loftee_docker
        memory: "8 GB"
        disks: "local-disk 500 HDD"
    }

    output {
        File output_file = "vep_loftee_raw_output.txt"
        File monitoring_log = "monitoring.log"
        File? warnings = "warnings.txt"
        Boolean done = true
    }
}

task BigQueryLoadRawVepAndLofteeAnnotations {
    # Loads data into BQ in the same "raw" tab-delimited format that was emitted by VEP.
    input {
        String variants_docker
        Array[File] vep_loftee_raw_output
        String project_id
        String dataset_name
        String raw_data_table
        File raw_data_table_schema
        Int runtime_cpu
        String runtime_memory
        String runtime_disks
        Int runtime_preemptible
    }

    meta {
        volatile: true
    }

    parameter_meta {
        vep_loftee_raw_output: {
            localization_optional: true
        }
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{raw_data_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating raw VEP + LOFTEE table ~{dataset_name}.~{raw_data_table}"

            # 3 day TTL for this table
            DATE=$((3 * 24 * 60 * 60))
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{raw_data_table} ~{raw_data_table_schema}
        fi

        num_rows=$(bq --apilog=false show --project_id=~{project_id} --format json ~{dataset_name}.~{raw_data_table} | jq -r .numRows)
        if ((num_rows != 0))
        then
            echo "Found preexisting table with data, not adding more raw data."
        else
            echo "Raw data table is empty, copying VEP output to be loaded."
            for file in ~{sep=' ' vep_loftee_raw_output}
            do
                gcloud storage cp $file .
                filename=$(basename $file)
                if [ ! -e load_file.txt ]
                then
                    # Do a wee bit of processing of the raw output to create a load file for raw VEP + LOFTEE data
                    # - Remove lines beginning with '##'.
                    # - Remove the leading '#' from the one line that should be left with a single leading '#' so
                    #   the line can serve as a TSV header.
                    sed -E '/^##/d' $filename | sed -E 's/^#//' > load_file.txt
                fi
                set +o errexit
                # In our integration tests when running on chromosomes 20/X/Y (the default), many of these files will be
                # empty, so turn off errexit temporarily.
                grep -E -v '^#' $filename >> load_file.txt
                set -o errexit
            done

            bq --apilog=false load --project_id=~{project_id} --source_format=CSV --field_delimiter='\t' \
                --skip_leading_rows=1 --null_marker="-" ~{dataset_name}.~{raw_data_table} load_file.txt

            echo "VEP + LOFTEE raw data loading complete."
        fi
    >>>

    runtime {
        docker: variants_docker
        cpu: runtime_cpu
        memory: runtime_memory
        disks: runtime_disks
        preemptible: runtime_preemptible
    }

    output {
        Boolean done = true
    }
}

task BigQueryCookVepAndLofteeRawAnnotations {
    # Reformat the "raw" data ("cooks" it) into a more directly queryable BigQuery format. This involves:
    # - Creating a VID field by parsing information in the `Uploaded_variation` field.
    # - Splitting nested fields into BigQuery REPEATEDs.
    # - Extracting an HGNC ID.
    # - Splitting and castng a nested GERP field into an array of floating point numbers.
    # - Squashing any duplicate rows resulting from deletions spanning shards.
    input {
        # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Boolean go
        String variants_docker
        String project_id
        String dataset_name
        String raw_data_table
        String cooked_data_table
        File cooked_data_table_schema
    }
    meta {
        volatile: true
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{cooked_data_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]
        then
            echo 'Creating "cooked" VEP + LOFTEE table ~{dataset_name}.~{cooked_data_table}'

            # 3 day TTL for this table
            DATE=$((3 * 24 * 60 * 60))
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{cooked_data_table} ~{cooked_data_table_schema}
        fi

        num_rows=$(bq --apilog=false show --project_id=~{project_id} --format json ~{dataset_name}.~{cooked_data_table} | jq -r .numRows)
        if [ $num_rows -ne 0]
        then
          echo "Found preexisting table with data, not adding more cooked data."
        else

        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{cooked_data_table} --replace \
           --project_id=~{project_id} '

        SELECT * EXCEPT(row_number) FROM (
            SELECT
            -- Make a VID-compatible string from the data in Uploaded_variation.
            -- VEP appears to use a different convention for the encoding of indel positions than what is used in GVS:
            -- VEP indel positions are based on the first *discrepant* base and not the first base mentioned, which in the
            -- GVS convention agrees between reference and allele. Correct for that in the VID-building code below to
            -- subtract 1 if the variant is an indel.
            REGEXP_EXTRACT(Uploaded_variation, "^chr([^_]+)") || "-" ||
                -- A Location specified with a "-" range is an indel. Single-base deletions are a special case with a single
                -- position, but like all deletions they have a NULL Allele so look for that as well.
                IF ((Location LIKE "%-%") OR (Allele is NULL),
                    -- If this is an indel decrement the position by one for VAT compatibility.
                    CAST((CAST(REGEXP_EXTRACT(Uploaded_variation, "_(\\d+)") AS INT64) - 1) AS STRING),
                    -- Else SNPs use position without adjustment.
                    REGEXP_EXTRACT(Uploaded_variation, "_(\\d+)")) ||
                "-" || REGEXP_EXTRACT(Uploaded_variation, "_([ACGT]+)/") || "-" ||
                REGEXP_EXTRACT(Uploaded_variation, "([ACGT]+)$") AS vid,
            Uploaded_variation,
            Location,
            Allele,
            Gene,
            Feature,
            Feature_type,
            Consequence,
            cDNA_position,
            CDS_position,
            Protein_position,
            Amino_acids,
            Codons,
            Existing_variation,
            IMPACT,
            DISTANCE,
            STRAND,
            -- FLAGS can be multi-valued so SPLIT to make this REPEATED.
            SPLIT(FLAGS, ",") AS FLAGS,
            SYMBOL as HGNC_SYMBOL,
            SYMBOL_SOURCE,
            -- HGNC IDs are formatted like HGNC:1234; we only want the number part.
            CAST(SPLIT(HGNC_ID, ":")[OFFSET(1)] AS INTEGER) AS HGNC_ID,
            SOURCE,
            LoF,
            -- These three appear to sometimes be multi-valued so SPLIT to make them REPEATEDs.
            SPLIT(LoF_filter, ",") AS LoF_filter,
            SPLIT(LoF_flags, ",") AS LoF_flags,
            SPLIT(LoF_info, ",") AS LoF_info,
            -- Split and cast the GERP string to REPEATED FLOAT64s.
            (
            SELECT
            ARRAY_AGG(SAFE_CAST(s AS FLOAT64))
            FROM
            UNNEST(SPLIT(GERP, ",")) AS s
            ) AS GERP,

            -- Use the ROW_NUMBER() magic to squash duplicates. A small number of deletions span interval boundaries
            -- and are assigned to two different VEP processing shards. This duplicate data would cause problems when
            -- we try to assign back
            ROW_NUMBER()
            -- The expression below uses Uploaded_variation rather than vid because BigQuery claims to not be able to
            -- find the vid identifier. Uploaded_variation contains equivalent information to vid in a different format.
            OVER (PARTITION BY Uploaded_variation, Feature)
            ROW_NUMBER

            FROM
            ~{project_id}.~{dataset_name}.~{raw_data_table}
            )

        WHERE ROW_NUMBER = 1

        '
        fi

    >>>

    runtime {
        docker: variants_docker
        memory: "7 GB"
        disks: "local-disk 1000 HDD"
    }

    output {
        Boolean done = true
        String cooked_table_name = cooked_data_table
    }
}


task AnnotateVCF {
    input {
        File input_vcf
        String output_annotated_file_name
        File custom_annotations_file
        String cromwell_root

        # Mentioning this path in the inputs section of the task combined with checking the 'Use reference disks' option
        # in Terra UI tells Cromwell to arrange for the Nirvana reference disk to be attached to this VM.
        #@ except: UnusedInput
        File summon_reference_disk =
            "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/MITOMAP_20200819.nsa.idx"

        String variants_nirvana_docker

        File omim_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/OMIM_20220516.nga"
        File cosmic_gene_fusion_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/COSMIC_GeneFusions_94.gfj"
        File primate_ai_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/PrimateAI_0.2.nsa"
        File primate_ai_annotations_idx = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/PrimateAI_0.2.nsa.idx"
        File splice_ai_annotations = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/SpliceAi_1.3.nsa"
        File splice_ai_annotations_idx = "gs://gcp-public-data--broad-references/hg38/v0/Nirvana/3.18.1_2024-03-06/SupplementaryAnnotation/GRCh38/SpliceAi_1.3.nsa.idx"
        Boolean use_reference_disk
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String annotation_json_name = output_annotated_file_name + ".json.gz"
    String gene_annotation_json_name = output_annotated_file_name + ".genes.json.gz"
    String positions_annotation_json_name = output_annotated_file_name + ".positions.json.gz"
    String nirvana_location = "/Nirvana/Nirvana.dll"
    String custom_creation_location = "/Nirvana/SAUtils.dll"
    String jasix_location = "/Nirvana/Jasix.dll"
    String path = "/Cache/GRCh38/Both"
    String path_supplementary_annotations = "/SupplementaryAnnotation/GRCh38"
    String path_reference = "/References/Homo_sapiens.GRCh38.Nirvana.dat"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        set +o errexit
        cat ~{custom_annotations_file} | grep -v '^#' > content_check_file.txt
        set -o errexit

        if [ ! -s content_check_file.txt ]; then
            echo "Found NO custom annotations in ~{custom_annotations_file} skipping annotation of input VCF"
            echo "Creating empty annotation jsons for subsequent tasks"
            touch ~{gene_annotation_json_name}
            touch ~{positions_annotation_json_name}
            exit 0
        fi

        if [[ "~{use_reference_disk}" == "true" ]]
        then
            # There's an issue with how the projects/broad-dsde-cromwell-dev/global/images/nirvana-3-18-1-references-2023-01-03
            # disk image was built: while all the reference files do exist on the image they are not at the expected
            # locations. The following code works around this issue and should continue to work even after a corrected
            # version of the Nirvana reference image is deployed into Terra.

            # Find where the reference disk should have been mounted on this VM.  Note this is referred to as a "candidate
            # mount point" because we do not actually confirm this is a reference disk until the following code block.
            if [[ -e /mnt/disks/cromwell_root/gcs_delocalization.sh ]]
            then
              # GCP Batch mount points
              CANDIDATE_MOUNT_POINT=$(lsblk | grep -v cromwell_root | sed -E -n 's!.*(/mnt/disks/[a-f0-9]+).*!\1!p')
            elif [[ -e /cromwell_root/gcs_delocalization.sh ]]
            then
              # PAPI mount points
              CANDIDATE_MOUNT_POINT=$(lsblk | grep -v cromwell_root | sed -E -n 's!.*(/mnt/[a-f0-9]+).*!\1!p')
            else
              >&2 echo "Could not find a mounted volume that looks like a reference disk, exiting."
              exit 1
            fi

            # Find one particular reference under the mount path. Note this is not the same reference as was specified in the
            # `inputs` section, so this would only be present if the volume we're looking at is in fact a reference disk.
            REFERENCE_FILE="Homo_sapiens.GRCh38.Nirvana.dat"
            REFERENCE_PATH=$(find ${CANDIDATE_MOUNT_POINT} -name "${REFERENCE_FILE}")
            if [[ -z ${REFERENCE_PATH} ]]; then
                >&2 echo "Could not find reference file '${REFERENCE_FILE}' under candidate reference disk mount point '${CANDIDATE_MOUNT_POINT}', exiting."
                exit 1
            fi

            # Take the parent of the parent directory of this file as root of the locally mounted references:
            DATA_SOURCES_FOLDER="$(dirname $(dirname ${REFERENCE_PATH}))"
        else
            DATA_SOURCES_FOLDER=~{cromwell_root}/nirvana_references
            mkdir ${DATA_SOURCES_FOLDER}

            # Download the references
            dotnet /Nirvana/Downloader.dll --ga GRCh38 --out ${DATA_SOURCES_FOLDER}

            # As of 2024-01-24 OMIM is no longer included among the bundle of annotation resources pulled down by the
            # Nirvana downloader. As this annotation set is currently central for our VAT logic, special-case link in
            # the OMIM .nsa bundle we downloaded back when we made the Delta reference disk:
            ln ~{omim_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            # Similarly, the following annotations were removed from the latest Nirvana annotations (3.18.1), but we
            # re-add them as desired by Lee
            ln ~{cosmic_gene_fusion_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{primate_ai_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{primate_ai_annotations_idx} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{splice_ai_annotations} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
            ln ~{splice_ai_annotations_idx} ${DATA_SOURCES_FOLDER}/SupplementaryAnnotation/GRCh38/
        fi

        # =======================================
        echo "Creating custom annotations"
        mkdir customannotations_dir
        CUSTOM_ANNOTATIONS_FOLDER="$PWD/customannotations_dir"

        # Add AC/AN/AF as custom annotations
        ## use --skip-ref once you are on a version of nirvana later than 3.14 (once they have created a docker image for it)
        dotnet ~{custom_creation_location} customvar \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -i ~{custom_annotations_file} \
            -o $CUSTOM_ANNOTATIONS_FOLDER

        # =======================================
        # Create Nirvana annotations:

        dotnet ~{nirvana_location} \
            -i ~{input_vcf} \
            -c $DATA_SOURCES_FOLDER~{path} \
            --sd $DATA_SOURCES_FOLDER~{path_supplementary_annotations} \
            --sd $CUSTOM_ANNOTATIONS_FOLDER \
            -r $DATA_SOURCES_FOLDER~{path_reference} \
            -o ~{output_annotated_file_name}

        # https://illumina.github.io/NirvanaDocumentation/introduction/parsing-json#jasix
        # Parse out the Genes section into a separate annotated json
        dotnet  ~{jasix_location} \
            --in ~{annotation_json_name} \
            --section genes \
            --out ~{gene_annotation_json_name}

        # Parse out the Positions section into a separate annotated json
        dotnet  ~{jasix_location} \
        --in ~{annotation_json_name} \
        --section positions \
        --out ~{positions_annotation_json_name}

    >>>

    runtime {
        docker: variants_nirvana_docker
        memory: "128 GB"
        cpu: 4
        preemptible: 1
        maxRetries: 1
        disks: "local-disk 2000 HDD"
    }

    output {
        File genes_annotation_json = "~{gene_annotation_json_name}"
        File positions_annotation_json = "~{positions_annotation_json_name}"
        File monitoring_log = "monitoring.log"
    }
}

task PrepVtAnnotationJson {
    input {
        File positions_annotation_json
        String output_file_suffix
        String output_path
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_vt_json = "vat_vt_bq_load." + output_file_suffix

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Kick off the monitoring script
        bash ~{monitoring_script} > monitoring.log &

        if [ ! -s ~{positions_annotation_json} ]; then
            echo "Input annotations json file is empty. Skipping further prep"
            touch ~{output_vt_json}
            exit 0
        fi

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_vt_bqloadjson_from_annotations.py \
            --annotated_json ~{positions_annotation_json} \
            --output_vt_json ~{output_vt_json}

        gsutil cp ~{output_vt_json} '~{output_path}'

    >>>

    runtime {
        docker: variants_docker
        memory: "16 GB"
        preemptible: 2
        cpu: 1
        disks: "local-disk 500 HDD"
    }

    output {
        Boolean done = true
        File monitoring_log = "monitoring.log"
    }
}

task PrepGenesAnnotationJson {
    input {
        File genes_annotation_json
        String output_file_suffix
        String output_path
        String variants_docker
    }

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

    String output_genes_json = "vat_genes_bq_load." + output_file_suffix

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        bash ~{monitoring_script} > monitoring.log &

        if [ ! -s ~{genes_annotation_json} ]; then
            echo "Input annotations json file is empty. Skipping further prep"
            touch ~{output_genes_json}
            exit 0
        fi

        ## the annotation jsons are split into the specific VAT schema
        python3 /app/create_genes_bqloadjson_from_annotations.py \
            --annotated_json ~{genes_annotation_json} \
            --output_genes_json ~{output_genes_json}

        gsutil cp ~{output_genes_json} '~{output_path}'

    >>>

    runtime {
        docker: variants_docker
        memory: "7 GB"
        preemptible: 3
        cpu: 1
        disks: "local-disk 500 HDD"
    }

    output {
        Boolean done = true
        File monitoring_log = "monitoring.log"
    }
}

task LoadManeDataIntoBigQuery {
    meta {
        # volatile: true so that call caching does not prevent this task from loading MANE data
        # into the current dataset. Without this, a cache hit would return the correct output
        # string ("mane_annotations") but skip the actual bq load, leaving the table absent.
        volatile: true
    }

    input {
        String project_id
        String dataset_name
        String mane_table_name
        File mane_data_file
        String cloud_sdk_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Remove the leading comment character on the first line so BigQuery will name the columns all nice.
        sed -i 's/^\#NCBI_GeneID/NCBI_GeneID/' ~{mane_data_file}

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{mane_table_name} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Loading MANE annotations into table ~{dataset_name}.~{mane_table_name}"
            bq --apilog=false load --project_id=~{project_id} --source_format=CSV --field_delimiter='\t' --skip_leading_rows=1 --autodetect ~{dataset_name}.~{mane_table_name} ~{mane_data_file}
        else
            echo "Found existing MANE annotations table ~{dataset_name}.~{mane_table_name}. Using it"
        fi
    >>>

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        preemptible: 3
        cpu: 1
        disks: "local-disk 100 HDD"
    }

    output {
        String mane_table = mane_table_name
        Boolean done = true
    }
}

task BigQueryLoadJson {
    meta {
        # since the WDL will not see the updated data (its getting put in a gcp bucket)
        volatile: true
    }

    input {
        String base_vat_table_name
        File vat_schema
        File variant_transcript_schema
        File genes_schema
        String mane_table_name
        String? vep_loftee_cooked_table_name
        Boolean run_vep_loftee_update
        String project_id
        String dataset_name
        String variant_transcripts_path
        String genes_path
        # Intentionally unused: this input exists solely to enforce task ordering - the upstream task's `done` output
        # is passed here to prevent this task from running until the upstream task has completed.
        #@ except: UnusedInput
        Array[Boolean] go
        String cloud_sdk_docker
    }

    # This is the name of the vat table. Due to sharding (VS-1191) there may be some duplicated entries.
    # So we create it here, and then deduplicate it in a later step
    String vat_table_name = base_vat_table_name + "_w_dups"

    String variant_transcript_table = base_vat_table_name + "_variants"
    String genes_table = base_vat_table_name + "_genes"

    String variant_transcripts_wildcarded_path = variant_transcripts_path + '*'
    String genes_wildcarded_path = genes_path + '*'

    String bq_labels = "--label service:gvs --label team:variants --label managedby:create_vat"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        if [[ "~{run_vep_loftee_update}" == "true" && -z "~{vep_loftee_cooked_table_name}" ]]; then
            echo "Error: 'run_vep_loftee_update' is true but 'vep_loftee_cooked_table_name' is not defined." 1>&2
            exit 1
        fi

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{variant_transcript_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{variant_transcript_table} ~{variant_transcript_schema}
        fi

        echo "Loading variant_transcript data into a pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        echo ~{variant_transcripts_wildcarded_path}
        bq --apilog=false load --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON ~{dataset_name}.~{variant_transcript_table} ~{variant_transcripts_wildcarded_path}

        echo "Adding the Mane SELECT annotation data to the pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} \
        'UPDATE `~{dataset_name}.~{variant_transcript_table}` vtt SET vtt.mane_select_name = mane.name FROM `~{dataset_name}.~{mane_table_name}` mane WHERE vtt.transcript = mane.Ensembl_nuc AND mane.MANE_status = "MANE Select" AND vtt.transcript is not null;'

        echo "Adding the Mane Plus Clinical annotation data to the pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} \
        'UPDATE `~{dataset_name}.~{variant_transcript_table}` vtt SET vtt.mane_plus_clinical_name = mane.name FROM `~{dataset_name}.~{mane_table_name}` mane WHERE vtt.transcript = mane.Ensembl_nuc AND mane.MANE_status = "MANE Plus Clinical" AND vtt.transcript is not null;'

        if [[ "~{run_vep_loftee_update}" == "true" ]]; then
        echo "Adding VEP + LOFTEE annotation data to the pre-vat table ~{dataset_name}.~{variant_transcript_table}"
        bq --apilog=false --project_id=~{project_id} query --format=csv --use_legacy_sql=false ~{bq_labels} '

        UPDATE `~{dataset_name}.~{variant_transcript_table}` vtt SET

        vtt.hgnc_symbol = vep.hgnc_symbol,
        vtt.hgnc_id     = vep.hgnc_id,
        vtt.LoF         = vep.LoF,
        vtt.LoF_filter  = vep.LoF_filter,
        vtt.LoF_flags   = vep.LoF_flags,
        vtt.LoF_info    = vep.LoF_info,
        vtt.GERP        = vep.GERP

        FROM `~{dataset_name}.~{vep_loftee_cooked_table_name}` vep WHERE

        vtt.transcript is not null AND
        vep.Feature_type is not null AND
        vtt.vid = vep.vid AND
        -- Do not consider version numbers when matching on transcripts. In Quickstart about 25% of the transcripts are
        -- mismatched on version number, with VEP having newer versions.
        SPLIT(vtt.transcript, ".")[OFFSET(0)] = SPLIT(vep.Feature, ".")[OFFSET(0)]

        '
        else
            echo "Skipping VEP + LOFTEE update because run_vep_loftee_update is false."
        fi

        set +o errexit
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{genes_table} > /dev/null
        BQ_SHOW_RC=$?
        set -o errexit

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating a pre-vat table ~{dataset_name}.~{genes_table}"
            bq --apilog=false mk --expiration=$DATE --project_id=~{project_id}  ~{dataset_name}.~{genes_table} ~{genes_schema}
        fi

        echo "Loading genes data into a pre-vat table ~{dataset_name}.~{genes_table}"
        echo ~{genes_wildcarded_path}
        bq --apilog=false load  --project_id=~{project_id} --source_format=NEWLINE_DELIMITED_JSON  ~{dataset_name}.~{genes_table} ~{genes_wildcarded_path}

        set +e
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{vat_table_name} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating the vat table ~{dataset_name}.~{vat_table_name}"
        else
            echo "Dropping and recreating the vat table ~{dataset_name}.~{vat_table_name}"
            bq --apilog=false rm -t -f --project_id=~{project_id} ~{dataset_name}.~{vat_table_name}
        fi

        CLUSTERING_STRING="--clustering_fields=contig"
        bq --apilog=false mk ${CLUSTERING_STRING} --expiration=$DATE --project_id=~{project_id} ~{dataset_name}.~{vat_table_name} ~{vat_schema}
        echo "Loading data into it"


        # Now we run a giant query in BQ to get this all in the right table and join the genes properly
        # Note the genes table join includes the group by to avoid the duplicates that get created from genes that span shards
        # Commented out columns in the query are to be added in the next release
        # We want the vat creation query to overwrite the destination table because if new data has been put into the pre-vat tables
        # and this workflow has been run an additional time, we dont want duplicates being appended from the original run

        # bq query --max_rows check: ok selecting into a table
        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{vat_table_name} --replace --project_id=~{project_id} \
        'SELECT
            v.vid,
            v.transcript,
            v.contig,
            v.position,
            v.ref_allele,
            v.alt_allele,
            v.gvs_all_ac,
            v.gvs_all_an,
            v.gvs_all_af,
            v.gvs_all_sc,
            v.gvs_max_af,
            v.gvs_max_ac,
            v.gvs_max_an,
            v.gvs_max_sc,
            v.gvs_max_subpop,
            v.gvs_afr_ac,
            v.gvs_afr_an,
            v.gvs_afr_af,
            v.gvs_afr_sc,
            v.gvs_amr_ac,
            v.gvs_amr_an,
            v.gvs_amr_af,
            v.gvs_amr_sc,
            v.gvs_eas_ac,
            v.gvs_eas_an,
            v.gvs_eas_af,
            v.gvs_eas_sc,
            v.gvs_eur_ac,
            v.gvs_eur_an,
            v.gvs_eur_af,
            v.gvs_eur_sc,
            v.gvs_mid_ac,
            v.gvs_mid_an,
            v.gvs_mid_af,
            v.gvs_mid_sc,
            v.gvs_oth_ac,
            v.gvs_oth_an,
            v.gvs_oth_af,
            v.gvs_oth_sc,
            v.gvs_sas_ac,
            v.gvs_sas_an,
            v.gvs_sas_af,
            v.gvs_sas_sc,
            v.gene_symbol,
            v.transcript_source,
            v.aa_change,
            v.consequence,
            v.dna_change_in_transcript,
            v.variant_type,
            v.exon_number,
            v.intron_number,
            v.genomic_location,
            # v.hgvsc AS splice_distance
            v.dbsnp_rsid,
            v.gene_id,
            # v.entrez_gene_id,
            # g.hgnc_gene_id,
            g.gene_omim_id,
            CASE WHEN ( v.transcript is not null and v.is_canonical_transcript is not True)
            THEN False WHEN ( v.transcript is not null and v.is_canonical_transcript is True) THEN True END AS is_canonical_transcript,
            v.gnomad_all_af,
            v.gnomad_all_ac,
            v.gnomad_all_an,
            v.gnomad_failed_filter,
            v.gnomad_max_af,
            v.gnomad_max_ac,
            v.gnomad_max_an,
            v.gnomad_max_subpop,
            v.gnomad_afr_ac,
            v.gnomad_afr_an,
            v.gnomad_afr_af,
            v.gnomad_amr_ac,
            v.gnomad_amr_an,
            v.gnomad_amr_af,
            v.gnomad_asj_ac,
            v.gnomad_asj_an,
            v.gnomad_asj_af,
            v.gnomad_eas_ac,
            v.gnomad_eas_an,
            v.gnomad_eas_af,
            v.gnomad_fin_ac,
            v.gnomad_fin_an,
            v.gnomad_fin_af,
            v.gnomad_nfe_ac,
            v.gnomad_nfe_an,
            v.gnomad_nfe_af,
            v.gnomad_sas_ac,
            v.gnomad_sas_an,
            v.gnomad_sas_af,
            v.gnomad_oth_ac,
            v.gnomad_oth_an,
            v.gnomad_oth_af,
            v.revel,
            v.splice_ai_acceptor_gain_score,
            v.splice_ai_acceptor_gain_distance,
            v.splice_ai_acceptor_loss_score,
            v.splice_ai_acceptor_loss_distance,
            v.splice_ai_donor_gain_score,
            v.splice_ai_donor_gain_distance,
            v.splice_ai_donor_loss_score,
            v.splice_ai_donor_loss_distance,
            g.omim_phenotypes_id,
            g.omim_phenotypes_name,
            v.clinvar_classification,
            v.clinvar_last_updated,
            v.clinvar_phenotype,
            v.clinvar_rcv_ids,
            v.clinvar_rcv_classifications,
            v.clinvar_rcv_num_stars,
            v.mane_select_name,
            v.mane_plus_clinical_name,
            v.hgnc_symbol,
            v.hgnc_id,
            v.LoF,
            v.LoF_filter,
            v.LoF_flags,
            v.LoF_info,
            v.GERP
        FROM `~{dataset_name}.~{variant_transcript_table}` as v
            left join
        (SELECT gene_symbol, ANY_VALUE(gene_omim_id) AS gene_omim_id, ANY_VALUE(omim_phenotypes_id) AS omim_phenotypes_id, ANY_VALUE(omim_phenotypes_name) AS omim_phenotypes_name FROM `~{dataset_name}.~{genes_table}` group by gene_symbol) as g
        on v.gene_symbol = g.gene_symbol'
    >>>

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        preemptible: 3
        cpu: 1
        disks: "local-disk 1000 HDD"
    }

    output {
        String vat_table = vat_table_name
        Boolean done = true
    }
}

task DeduplicateVatInBigQuery {
    meta {
        # since the WDL will not see the updated data (it's getting put in a gcp bucket)
        volatile: true
    }

    input {
        String input_vat_table_name
        String output_vat_table_name
        File vat_schema

        String project_id
        String dataset_name
        String cloud_sdk_docker
    }


    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        echo "project_id = ~{project_id}" > ~/.bigqueryrc

        DATE=86400 ## 24 hours in seconds

        set +e
        bq --apilog=false show --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name} > /dev/null
        BQ_SHOW_RC=$?
        set -e

        if [ $BQ_SHOW_RC -ne 0 ]; then
            echo "Creating the final vat table ~{dataset_name}.~{output_vat_table_name}"
        else
            bq --apilog=false rm -t -f --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name}
        fi
        bq --apilog=false mk --project_id=~{project_id} ~{dataset_name}.~{output_vat_table_name} ~{vat_schema}
        echo "Loading data into it"

        # Now we query the original VAT table and recreate it, but remove any rows that appear twice.

        # bq query --max_rows check: ok selecting into a table
        bq --apilog=false query --nouse_legacy_sql --destination_table=~{dataset_name}.~{output_vat_table_name} --replace --project_id=~{project_id} \
        ' SELECT * EXCEPT(row_number) FROM (
            SELECT
                *,
                row_number()
                    over (partition by vid, transcript)
                    row_number
            FROM
                `~{dataset_name}.~{input_vat_table_name}`
            )
            where row_number = 1'
    >>>

    runtime {
        docker: cloud_sdk_docker
        memory: "3 GB"
        preemptible: 3
        cpu: 1
        disks: "local-disk 100 HDD"
    }

    output {
        String vat_table = output_vat_table_name
        Boolean done = true
    }
}
