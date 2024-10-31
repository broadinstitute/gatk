version 1.0

import "../GvsUtils.wdl" as Utils
import "../GvsExtractAvroFilesForHail.wdl" as ExtractAvroFilesForHail
import "../GvsCreateVDS.wdl" as CreateVds
import "GvsQuickstartVcfIntegration.wdl" as QuickstartVcfIntegration

workflow GvsQuickstartHailIntegration {
    input {
        String git_branch_or_tag
        String? git_hash
        Boolean is_wgs
        File? interval_list
        Boolean use_compressed_references = false
        Boolean extract_do_not_filter_override
        String dataset_suffix = "hail"
        Boolean use_default_dockers = false
        Boolean bgzip_output_vcfs = false

        String? basic_docker
        String? cloud_sdk_docker
        String? cloud_sdk_slim_docker
        String? variants_docker
        String? gatk_docker

        File? gatk_override
        String? hail_version
        String expected_output_prefix
        String sample_id_column_name
        String vcf_files_column_name
        String vcf_index_files_column_name
        String? sample_set_name ## NOTE: currently we only allow the loading of one sample set at a time

        String? workspace_bucket
        String? workspace_id
        String? submission_id

        Int? maximum_alternate_alleles
    }

    String project_id = "gvs-internal"

    if (!defined(workspace_bucket) || !defined(workspace_id) || !defined(submission_id) ||
        !defined(git_hash) || !defined(basic_docker) || !defined(cloud_sdk_docker) || !defined(cloud_sdk_slim_docker) ||
        !defined(variants_docker) || !defined(gatk_docker) || !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
    String effective_git_hash = select_first([git_hash, GetToolVersions.git_hash])
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_workspace_id = select_first([workspace_id, GetToolVersions.workspace_id])
    String effective_submission_id = select_first([submission_id, GetToolVersions.submission_id])

    call QuickstartVcfIntegration.GvsQuickstartVcfIntegration {
        input:
            git_branch_or_tag = git_branch_or_tag,
            git_hash = git_hash,
            drop_state = "ZERO",
            use_compressed_references = use_compressed_references,
            extract_do_not_filter_override = extract_do_not_filter_override,
            load_vcf_headers = true,
            dataset_suffix = dataset_suffix,
            use_default_dockers = use_default_dockers,
            check_expected_cost_and_table_size_outputs = false,
            gatk_override = gatk_override,
            bgzip_output_vcfs = bgzip_output_vcfs,
            is_wgs = is_wgs,
            interval_list = interval_list,
            expected_output_prefix = expected_output_prefix,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            variants_docker = effective_variants_docker,
            gatk_docker = effective_gatk_docker,
            workspace_bucket = effective_workspace_bucket,
            workspace_id = effective_workspace_id,
            submission_id = effective_submission_id,
            maximum_alternate_alleles = maximum_alternate_alleles,
    }

    call ExtractAvroFilesForHail.GvsExtractAvroFilesForHail {
        input:
            go = GvsQuickstartVcfIntegration.done,
            git_branch_or_tag = git_branch_or_tag,
            git_hash = git_hash,
            project_id = project_id,
            dataset_name = GvsQuickstartVcfIntegration.dataset_name,
            filter_set_name = GvsQuickstartVcfIntegration.filter_set_name,
            scatter_width = 10,
            call_set_identifier = git_branch_or_tag,
            basic_docker = effective_basic_docker,
            cloud_sdk_docker = effective_cloud_sdk_docker,
            variants_docker = effective_variants_docker,
    }

    String vds_path = GvsExtractAvroFilesForHail.avro_path + "/gvs_export.vds"

    call CreateVds.GvsCreateVDS {
        input:
            git_branch_or_tag = git_branch_or_tag,
            hail_version = effective_hail_version,
            avro_path = GvsExtractAvroFilesForHail.avro_path,
            vds_destination_path = vds_path,
            cluster_prefix = "vds-cluster",
            gcs_subnetwork_name = "subnetwork",
            region = "us-central1",
            basic_docker = effective_basic_docker,
            variants_docker = effective_variants_docker,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            cluster_max_age_minutes = 120,
            cluster_max_idle_minutes = 60,
            leave_cluster_running_at_end = false,
    }

    call TieOutVds {
        input:
            go = GvsCreateVDS.done,
            git_branch_or_tag = git_branch_or_tag,
            vds_path = vds_path,
            tieout_vcfs = GvsQuickstartVcfIntegration.output_vcfs,
            tieout_vcf_indexes = GvsQuickstartVcfIntegration.output_vcf_indexes,
            tieout_vcf_suffix = if (bgzip_output_vcfs) then ".bgz" else ".gz",
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            hail_version = effective_hail_version,
    }

    output {
        Array[File] output_vcfs = GvsQuickstartVcfIntegration.output_vcfs
        Array[File] output_vcf_indexes = GvsQuickstartVcfIntegration.output_vcf_indexes
        Float total_vcfs_size_mb = GvsQuickstartVcfIntegration.total_vcfs_size_mb
        File manifest = GvsQuickstartVcfIntegration.manifest
        String vds_output_path = GvsExtractAvroFilesForHail.avro_path + "/gvs_export.vds"
        String recorded_git_hash = effective_git_hash
        Boolean done = true
        Boolean used_tighter_gcp_quotas = GvsQuickstartVcfIntegration.used_tighter_gcp_quotas
    }
}


task TieOutVds {
    input {
        Boolean go
        String git_branch_or_tag
        String vds_path
        Array[File] tieout_vcfs
        Array[File] tieout_vcf_indexes
        String tieout_vcf_suffix
        String cloud_sdk_slim_docker
        String hail_version
    }
    parameter_meta {
        tieout_vcfs: {
             localization_optional: true
         }
        tieout_vcf_indexes: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Copy the versions of the Hail import and tieout scripts for this branch from GitHub.
        script_url_prefix="https://raw.githubusercontent.com/broadinstitute/gatk/~{git_branch_or_tag}/scripts/variantstore/wdl/extract"
        for script in hail_gvs_import.py hail_join_vds_vcfs.py gvs_vds_tie_out.py import_gvs.py
        do
            curl --silent --location --remote-name "${script_url_prefix}/${script}"
        done

        # Create a manifest of VCFs and indexes to bulk download with `gcloud storage cp`.
        touch vcf_manifest.txt
        # This is extremely noisy and not interesting, turn off xtrace.
        set +o xtrace
        for file in ~{sep=' ' tieout_vcfs} ~{sep=' ' tieout_vcf_indexes}
        do
            echo $file >> vcf_manifest.txt
        done
        # xtrace back on
        set -o xtrace

        # Copy VCFs and indexes to the current directory.
        cat vcf_manifest.txt | gcloud storage cp -I .

        export REFERENCES_PATH=$PWD/references
        mkdir -p ${REFERENCES_PATH}

        gcloud storage cp 'gs://hail-common/references/Homo_sapiens_assembly38.fasta*' ${REFERENCES_PATH}

        # Versions of Hail near 0.2.117 demand Java 8 or Java 11, and refuse to run on Java 17. (This is because Google Dataproc is still on Java 11)
        # Temurin Java 8
        apt-get -qq install wget apt-transport-https gnupg
        wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
        echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
        apt-get -qq update
        apt -qq install -y temurin-8-jdk

        export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
        pip install --upgrade pip
        pip install hail==~{hail_version}

        export WORK=$PWD/work
        mkdir ${WORK}

        export TEMP_PATH=$WORK/temp
        mkdir ${TEMP_PATH}

        export VDS_PATH=$WORK/gvs_export.vds
        mkdir ${VDS_PATH}

        gcloud storage cp -r ~{vds_path}  ${WORK}


        export JOINED_MATRIX_TABLE_PATH=${WORK}/joined.mt

        python3 ./hail_join_vds_vcfs.py --vds-path ${VDS_PATH} --joined-matrix-table-path ${JOINED_MATRIX_TABLE_PATH} *.vcf~{tieout_vcf_suffix}

        pip install pytest
        ln -s ${WORK}/joined.mt .
        pytest ./gvs_vds_tie_out.py
    >>>
    runtime {
        # `slim` here to be able to use Java
        docker: cloud_sdk_slim_docker
        maxRetries: 2
        disks: "local-disk 2000 HDD"
        memory: "30 GiB"
    }
    output {
        Boolean done = true
    }
}
