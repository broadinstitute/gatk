version 1.0

import "GvsUtils.wdl" as Utils
import "GvsPrepareBulkImport.wdl" as PrepareBulkImport
import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes

workflow GvsBulkIngestGenomes {
    input {
        # Begin GvsPrepareBulkImport
        String terra_project_id
        String workspace_name
        String workspace_bucket
        String samples_table_name
        String sample_id_column_name # is this just <samples_table_name>_id ?
        String vcf_files_column_name
        String vcf_index_files_column_name
        # End GvsPrepareBulkImport

        # Begin GvsAssignIds  ## TODO what should the behavior be if we have run this once already and several samples did not load properly?
        String dataset_name
        String bq_project_id
        String call_set_identifier

        ## Array[String] external_sample_names ## TODO no longer an input param?

        File? gatk_override
        # End GvsAssignIds

        # Begin GvsImportGenomes
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

        ## TODO do we want to add Alt Allele Table creation?
    }

    call PrepareBulkImport.GvsPrepareBulkImport as PrepareBulkImport {
        input:
            project_id = terra_project_id,
            workspace_name = workspace_name,
            workspace_bucket  = workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name
    }

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            ## go = PrepareBulkImport.done, ## TODO we're gonna want to add this, right?
            dataset_name = dataset_name,
            project_id = bq_project_id,
            external_sample_names = read_lines(PrepareBulkImport.sampleFOFN),
            samples_are_controls = false ## TODO this shouldn't always be false tho
    }

    call ImportGenomes.GvsImportGenomes as ImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = bq_project_id,

            external_sample_names = read_lines(PrepareBulkImport.sampleFOFN),
            input_vcfs = read_lines(PrepareBulkImport.vcfFOFN),
            input_vcf_indexes = read_lines(PrepareBulkImport.vcfIndexFOFN),

            is_rate_limited_beta_customer = true,

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
        ## Boolean done = true
        Array[File] load_data_stderrs = ImportGenomes.load_data_stderrs
    }
}
