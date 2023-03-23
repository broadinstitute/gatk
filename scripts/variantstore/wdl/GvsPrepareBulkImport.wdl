version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsPrepareBulkImport {
    input {
        String project_id
        String workspace_name
        String workspace_bucket
        String samples_table_name
        String sample_id_column_name
        String vcf_files_column_name
        String vcf_index_files_column_name
    }

    call GenerateFOFNsFromDataTables {
        input:
            project_id = project_id,
            workspace_name = workspace_name,
            workspace_bucket  = workspace_bucket,
            samples_table_name = samples_table_name,
            sample_id_column_name = sample_id_column_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name
    }

    output {
        File sampleFOFN = GenerateFOFNsFromDataTables.sampleFOFN
        File vcfFOFN = GenerateFOFNsFromDataTables.vcfFOFN
        File vcfIndexFOFN = GenerateFOFNsFromDataTables.vcfIndexFOFN
        File errorRows = GenerateFOFNsFromDataTables.errors
    }
}

# This copies the entirety of the data tables over.  We'll want a different task to make FOFNs for just a sample set
task GenerateFOFNsFromDataTables {
    input {
        String project_id
        String workspace_name
        String workspace_bucket
        String samples_table_name = "sample"
        String sample_id_column_name = "sample_id"
        String vcf_files_column_name = "hg38_reblocked_v2_vcf"
        String vcf_index_files_column_name = "hg38_reblocked_v2_vcf_index"
    }

    String sample_names_file_name = "sample_names.txt"
    String vcf_files_name = "vcf_files.txt"
    String vcf_index_files_name = "vcf_index_files.txt"
    String error_file_name = "errors.txt"

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        export GOOGLE_PROJECT='~{project_id}'
        export WORKSPACE_NAMESPACE='~{project_id}'
        export WORKSPACE_NAME='~{workspace_name}'
        export WORKSPACE_BUCKET='~{workspace_bucket}'

        python3 /app/generate_FOFNs_for_import.py \
            --data_table_name ~{samples_table_name} \
            --sample_id_column_name ~{sample_id_column_name} \
            --vcf_files_column_name ~{vcf_files_column_name} \
            --vcf_index_files_column_name ~{vcf_index_files_column_name}

    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-03-23-bulk-ingest"
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
