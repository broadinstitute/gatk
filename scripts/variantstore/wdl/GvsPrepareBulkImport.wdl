version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsPrepareBulkImport {
    input {
        String project_id
        String workspace_name
        String workspace_namespace
        String workspace_bucket
        String samples_table_name
        String user_defined_sample_id_column_name
        String vcf_files_column_name
        String vcf_index_files_column_name
        String? sample_set_name
    }

    call GenerateFOFNsFromDataTables {
        input:
            google_project_id = project_id,
            workspace_name = workspace_name,
            workspace_namespace = workspace_namespace,
            workspace_bucket  = workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = user_defined_sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name
    }

    output {
        File sampleFOFN = GenerateFOFNsFromDataTables.sampleFOFN
        File vcfFOFN = GenerateFOFNsFromDataTables.vcfFOFN
        File vcfIndexFOFN = GenerateFOFNsFromDataTables.vcfIndexFOFN
        File errorRows = GenerateFOFNsFromDataTables.errors
    }
}

task GenerateFOFNsFromDataTables {
    ## In order to get the <entity>_ids in the sample_set for an inclusion list, we use Terra Notebook Utils
    ## This also allows us to validate that the requested sample_set exists
    input {
        String google_project_id
        String workspace_name
        String workspace_namespace
        String workspace_bucket
        String samples_table_name
        String sample_id_column_name ## NOTE: if the user has specified a different sample name column for GVS, it needs to be used independently of the sample_set info
        String vcf_files_column_name
        String vcf_index_files_column_name
        String? sample_set_name
    }

    String sample_names_file_name = "sample_names.txt"
    String vcf_files_name = "vcf_files.txt"
    String vcf_index_files_name = "vcf_index_files.txt"
    String error_file_name = "errors.txt"

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        export GOOGLE_PROJECT='~{google_project_id}'
        export WORKSPACE_NAMESPACE='~{workspace_namespace}'
        export WORKSPACE_NAME='~{workspace_name}'
        export WORKSPACE_BUCKET='~{workspace_bucket}'

        python3 /app/generate_FOFNs_for_import.py \
            --data_table_name ~{samples_table_name} \
            --sample_id_column_name ~{sample_id_column_name} \
            --vcf_files_column_name ~{vcf_files_column_name} \
            --vcf_index_files_column_name ~{vcf_index_files_column_name} \
            ~{"--sample_set_name " + sample_set_name} \
            --sample_names_file_name ~{sample_names_file_name} \
            --vcf_files_name ~{vcf_files_name} \
            --vcf_index_files_name ~{vcf_index_files_name} \
            --error_file_name ~{error_file_name}

        ## Validate by testing file lengths and failing if they are not all the same
        sample_count=$(wc -l < ~{sample_names_file_name})
        vcf_count=$(wc -l < ~{vcf_files_name})
        index_count=$(wc -l < ~{vcf_index_files_name})

        if [[ $sample_count -eq $vcf_count && $sample_count -eq $index_count ]]; then
            echo $sample_count
        else
            echo "Error: mismatched sample / VCF / index counts"
            echo "sample count: $sample_count"
            echo "vcf count: $vcf_count"
            echo "index count: $index_count"
            exit 1
        fi

    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-06-23-alpine"
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        File sampleFOFN = sample_names_file_name
        File vcfFOFN = vcf_files_name
        File vcfIndexFOFN = vcf_index_files_name
        File errors = error_file_name
    }
}
