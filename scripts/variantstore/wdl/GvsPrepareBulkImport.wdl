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

    call GenerateImportFofnFromDataTable {
        input:
            google_project_id = project_id,
            workspace_name = workspace_name,
            workspace_namespace = workspace_namespace,
            workspace_bucket  = workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = user_defined_sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name,
            sample_set_name = sample_set_name,
    }

    call SplitImportFofn {
        input:
            import_fofn = GenerateImportFofnFromDataTable.output_fofn,
    }

    output {
        Array[String] external_sample_names = SplitImportFofn.external_sample_names
        Array[String] vcf_file_names = SplitImportFofn.vcf_file_names
        Array[String] vcf_index_file_names = SplitImportFofn.vcf_index_file_names
        File errorRows = GenerateImportFofnFromDataTable.errors
    }
}

task GenerateImportFofnFromDataTable {
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

    String output_fofn_name = "output.tsv"
    String error_file_name = "errors.txt"

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        PS4='\D{+%F %T} \w $ '

        export GOOGLE_PROJECT='~{google_project_id}'
        export WORKSPACE_NAMESPACE='~{workspace_namespace}'
        export WORKSPACE_NAME='~{workspace_name}'
        export WORKSPACE_BUCKET='~{workspace_bucket}'

        python3 /app/generate_fofn_for_import.py \
            --data-table-name ~{samples_table_name} \
            --sample-id-column-name ~{sample_id_column_name} \
            --vcf-files-column-name ~{vcf_files_column_name} \
            --vcf-index-files-column-name ~{vcf_index_files_column_name} \
            ~{"--sample-set-name " + sample_set_name} \
            --output-file-name ~{output_fofn_name} \
            --error-file-name ~{error_file_name}
    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-07-07-alpine-13fa284e9"
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        File output_fofn = output_fofn_name
        File errors = error_file_name
    }
}

task SplitImportFofn {
    input {
        File import_fofn
    }

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail
        PS4='\D{+%F %T} \w $ '

        cut -f 1 ~{import_fofn} > sample_names.txt
        cut -f 2 ~{import_fofn} > vcf_file_names.txt
        cut -f 3 ~{import_fofn} > vcf_index_file_names.txt
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-07-07-alpine-13fa284e9"
        memory: "3 GB"
        disks: "local-disk 200 HDD"
        cpu: 1
    }

    output {
        Array[String] external_sample_names = read_lines("sample_names.txt")
        Array[String] vcf_file_names = read_lines("vcf_file_names.txt")
        Array[String] vcf_index_file_names = read_lines("vcf_index_file_names.txt")
    }
}