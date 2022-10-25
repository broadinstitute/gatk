version 1.0

import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes
import "GvsPopulateAltAllele.wdl" as PopulateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsExtractAvroFilesForHail.wdl" as ExtractAvroFilesForHail

workflow GvsUnifiedHail {
    input {
        # Begin GvsAssignIds
        String dataset_name
        String project_id
        String call_set_identifier

        Array[String] external_sample_names
        String call_set_identifier

        File? gatk_override
        # End GvsAssignIds

        # Begin GvsImportGenomes
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"

        # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
        String drop_state = "NONE"

        # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
        # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
        # those should work fine in most cases.
        Int? load_data_batch_size
        Int? load_data_preemptible_override
        Int? load_data_maxretries_override
        # End GvsImportGenomes

        # Begin GvsCreateFilterSet
        String filter_set_name = call_set_identifier
        Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
        Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]

        Int? INDEL_VQSR_max_gaussians_override = 4
        Int? INDEL_VQSR_mem_gb_override
        Int? SNP_VQSR_max_gaussians_override = 6
        Int? SNP_VQSR_mem_gb_override
        # End GvsCreateFilterSet

        # Begin CreateVds
        File hail_wheel = "gs://gvs-internal-scratch/hail-wheels/2022-10-18/0.2.102-964bee061eb0/hail-0.2.102-py3-none-any.whl"
        # End CreateVds
    }

#    call AssignIds.GvsAssignIds as AssignIds {
#        input:
#            dataset_name = dataset_name,
#            project_id = project_id,
#            external_sample_names = external_sample_names,
#            assign_ids_gatk_override = gatk_override
#    }
#
#    call ImportGenomes.GvsImportGenomes {
#        input:
#            go = AssignIds.done,
#            dataset_name = dataset_name,
#            project_id = project_id,
#            external_sample_names = external_sample_names,
#            input_vcfs = input_vcfs,
#            input_vcf_indexes = input_vcf_indexes,
#            interval_list = interval_list,
#            load_data_preemptible_override = load_data_preemptible_override,
#            load_data_maxretries_override = load_data_maxretries_override,
#            load_data_gatk_override = gatk_override,
#            load_data_batch_size = load_data_batch_size,
#            drop_state = drop_state
#    }
#
#    call PopulateAltAllele.GvsPopulateAltAllele {
#        input:
#            call_set_identifier = call_set_identifier,
#            go = GvsImportGenomes.done,
#            dataset_name = dataset_name,
#            project_id = project_id
#    }
#
#    call CreateFilterSet.GvsCreateFilterSet {
#        input:
#            go = GvsPopulateAltAllele.done,
#            dataset_name = dataset_name,
#            project_id = project_id,
#            call_set_identifier = call_set_identifier,
#            filter_set_name = filter_set_name,
#            indel_recalibration_annotation_values = indel_recalibration_annotation_values,
#            snp_recalibration_annotation_values = snp_recalibration_annotation_values,
#            interval_list = interval_list,
#            gatk_override = gatk_override,
#            INDEL_VQSR_max_gaussians_override = INDEL_VQSR_max_gaussians_override,
#            INDEL_VQSR_mem_gb_override = INDEL_VQSR_mem_gb_override,
#            SNP_VQSR_max_gaussians_override = SNP_VQSR_max_gaussians_override,
#            SNP_VQSR_mem_gb_override = SNP_VQSR_mem_gb_override
#    }
#
#    call ExtractAvroFilesForHail.GvsExtractAvroFilesForHail {
#        input:
#            go = GvsCreateFilterSet.done,
#            project_id = project_id,
#            dataset = dataset_name,
#            filter_set_name = filter_set_name,
#            scatter_width = 10,
#    }

    call CreateVds {
        input:
            # hail_gvs_import_script = GvsExtractAvroFilesForHail.hail_gvs_import_script,
            hail_gvs_import_script = "gs://fc-a1621719-20ea-471d-a0ef-a41383dc76bd/submissions/e7e3bf72-c849-4845-b0c2-f324569afb69/GvsQuickstartHailIntegration/e23f1bd3-29b5-43b2-bd53-f22f1e1509d5/call-GvsUnifiedHail/GvsUnifiedHail/a783025e-e639-42e7-981c-cc7ff4d923de/call-GvsExtractAvroFilesForHail/GvsExtractAvroFilesForHail/94a25fb5-5913-4432-9c2a-f1efdc0160f8/call-GenerateHailScripts/hail_gvs_import.py",
            hail_wheel = hail_wheel,
    }

    output {
        Boolean done = true
        # String vds_output_path = GvsExtractAvroFilesForHail.vds_output_path
    }
}

task CreateVds {
    input {
        File hail_gvs_import_script
        File hail_wheel
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Temurin Java 8
        apt-get -qq install wget apt-transport-https gnupg
        wget -O - https://packages.adoptium.net/artifactory/api/gpg/key/public | apt-key add -
        echo "deb https://packages.adoptium.net/artifactory/deb $(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | tee /etc/apt/sources.list.d/adoptium.list
        apt-get -qq update
        apt -qq install -y temurin-8-jdk

        pip install ~{hail_wheel}
        export PYSPARK_SUBMIT_ARGS='--driver-memory 16g --executor-memory 16g pyspark-shell'
        curl --silent --show-error --location https://broad.io/install-gcs-connector | python3

        python3 ~{hail_gvs_import_script}
    >>>
    runtime {
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:406.0.0-slim"
        disks: "local-disk 500 HDD"
        memory: "30 GiB"
    }
    output {
        Boolean done = true
    }
}
