version 1.0

import "GvsAssignIds.wdl" as AssignIds
import "GvsImportGenomes.wdl" as ImportGenomes
import "GvsCreateAltAllele.wdl" as CreateAltAllele
import "GvsCreateFilterSet.wdl" as CreateFilterSet
import "GvsPrepareRangesCallset.wdl" as PrepareRangesCallset
import "GvsExtractCallset.wdl" as ExtractCallset

workflow GvsJointVariantCalling {
    input {
        String dataset_name
        String project_id
        Array[String] external_sample_names
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        String callset_identifier
        String extract_output_gcs_dir
    }

    String extract_output_file_base_name = callset_identifier ## TODO make sure there are no spaces here!??!
    String extract_table_prefix = callset_identifier ## TODO make sure there are no spaces here!??!

    call AssignIds.GvsAssignIds as AssignIds {
        input:
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names
    }

    call ImportGenomes.GvsImportGenomes {
        input:
            go = AssignIds.done,
            dataset_name = dataset_name,
            project_id = project_id,
            external_sample_names = external_sample_names,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes
    }

    call CreateAltAllele.GvsCreateAltAllele {
        input:
            go = GvsImportGenomes.done,
            dataset_name = dataset_name,
            project_id = project_id
    }

    call CreateFilterSet.GvsCreateFilterSet {
        input:
            go = GvsCreateAltAllele.done,
            dataset_name = dataset_name,
            project_id = project_id,
            filter_set_name = callset_identifier
    }

    call PrepareRangesCallset.GvsPrepareCallset {
        input:
            go = GvsCreateFilterSet.done,
            project_id = project_id,
            dataset_name = dataset_name,
            extract_table_prefix = callset_identifier
    }

    call ExtractCallset.GvsExtractCallset {
        input:
            go = GvsPrepareCallset.done,
            dataset_name = dataset_name,
            project_id = project_id,
            extract_table_prefix = callset_identifier,
            filter_set_name = callset_identifier,
            output_file_base_name = callset_identifier,
            output_gcs_dir = extract_output_gcs_dir
    }

    output {
        Array[File] output_vcfs = GvsExtractCallset.output_vcfs
        Array[File] output_vcf_indexes = GvsExtractCallset.output_vcf_indexes
        Float total_vcfs_size_mb = GvsExtractCallset.total_vcfs_size_mb
        File manifest = GvsExtractCallset.manifest
        Boolean done = true
    }
}
