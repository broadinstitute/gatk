version 1.0

import "GvsUtils.wdl" as Utils

workflow GvePrepareBulkImport {
    input {
        Boolean go = true
        String samples_table_name = "sample"
        String vcf_files_column_name = "hg38_reblocked_v2_vcf"
        String vcf_index_files_column_name = "hg38_reblocked_v2_vcf_index"
    }

    call GenerateFOFNsFromDataTables {
        input:
            samples_table_name = samples_table_name,
            vcf_files_column_name = vcf_files_column_name,
            vcf_index_files_column_name = vcf_index_files_column_name
    }


    output {
        Boolean done = false
    }
}

task GenerateFOFNsFromDataTables {
    input {
        Boolean go = true
        String samples_table_name = "sample"
        String vcf_files_column_name = "hg38_reblocked_v2_vcf"
        String vcf_index_files_column_name = "hg38_reblocked_v2_vcf_index"
    }

    String sample_names_file_name = "sample_names.txt"
    String vcf_files_name = "vcf_files.txt"
    String vcf_index_files_name = "vcf_index_files.txt"

    meta {
        # because this is being used to determine the current state of the GVS database, never use call cache
        volatile: false
    }

    command <<<
        set -o errexit -o nounset -o xtrace -o pipefail

        python3 /app/generate_FOFNs_for_import.py \
        --data_table_name ~{samples_table_name} \
        --vcf_files_column_name ~{vcf_files_column_name} \
        --vcf_index_files_column_name ~{vcf_index_files_column_name}}

    >>>
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore:2023-1-13-FOFN_1"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        cpu: 1
    }

    output {
        File sampleFOFN = sample_names_file_name
        File vcfFOFN = vcf_files_name
        File vcfIndexFOFN = vcf_index_files_name
    }
}
