version 1.0

import "GvsUtils.wdl" as Utils
import "GvsPrepareBulkImport.wdl" as PrepareBulkImport
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes

workflow GvsBulkIngestGenomes {
    input {
        # Begin GvsPrepareBulkImport
        String project_id
        String workspace_name
        String workspace_bucket
        String samples_table_name = "sample"
        String sample_id_column_name = "sample_id"
        String vcf_files_column_name = "hg38_reblocked_v2_vcf" ## TODO should these stay hardcoded!?!?!
        String vcf_index_files_column_name = "hg38_reblocked_v2_vcf_index"
        # End GvsPrepareBulkImport

        # Begin GvsAssignIds
        String dataset_name
        String project_id
        String call_set_identifier

        ## Array[String] external_sample_names TODO no longer an input param?

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
    }

    call PrepareBulkImport {
        input:
            project_id = project_id,
            workspace_name = workspace_name,
            workspace_bucket  = workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name
    }

    call ArrayStringify {

    }

    call AssignIds {
        input:
            go = ArrayStringify.done,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = ArrayStringify.external_sample_names, ## Do I need the array string format here?
            samples_are_controls = false ## interesting that this isn't the DEFAULT default
    }

    call ImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = project_id,

            external_sample_names = ArrayStringify.external_sample_names, ## Do I need the array string format here?
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,

            skip_loading_vqsr_fields = false

            # set to "NONE" to ingest all the reference data into GVS for VDS (instead of VCF) output
            drop_state = "NONE"

            interval_list = interval_list,

            # The larger the `load_data_batch_size` the greater the probability of preemptions and non-retryable
            # BigQuery errors so if specifying this adjust preemptible and maxretries accordingly. Or just take the defaults,
            # those should work fine in most cases.
            load_data_batch_size = load_data_batch_size,
            load_data_maxretries_override = load_data_maxretries_override,
            load_data_preemptible_override = load_data_preemptible_override,
            load_data_gatk_override = gatk_override,

    }

    output {
        Boolean done = true
        Array[File] load_data_stderrs = LoadData.stderr
    }
}

task ArrayStringify {
    input {
        String project_id
        String workspace_name
        String workspace_bucket
    }


    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        ## stuff

    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-1-20-FOFN"
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
         Array[String] external_sample_names = external_sample_names
    }
}
