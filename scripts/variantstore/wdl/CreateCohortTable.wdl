version 1.0

workflow CreateCohortTable {
   input {
        String project
        String dataset

        String? destination_cohort_table_name
        String? fq_cohort_sample_table
        String? fq_sample_mapping_table

        String? docker
    }

    String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

    call CreateCohortTableTask {
        input:
            project                         = project,
            dataset                         = dataset,
            destination_cohort_table_name   = destination_cohort_table_name,
            fq_cohort_sample_table          = fq_cohort_sample_table,
            fq_sample_mapping_table         = fq_sample_mapping_table,
            docker                          = docker_final
    }

}

task CreateCohortTableTask {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    input {
        String project
        String dataset

        String? destination_cohort_table_name
        String? fq_cohort_sample_table
        String? fq_sample_mapping_table

        String docker
    }

    String destination_cohort_table_name = if defined(destination_cohort_table_name) then "${destination_cohort_table_name}" else "exported_cohort_all_samples"
    String fq_cohort_sample_table = if defined(fq_cohort_sample_table) then "${fq_cohort_sample_table}" else "${project}.${dataset}.sample_info"
    String fq_sample_mapping_table = if defined(fq_sample_mapping_table) then "${fq_sample_mapping_table}" else "${project}.${dataset}.sample_info"

    command <<<
        set -e

        python scripts/variantstore/wdl/extract/create_cohort_data_table.py \
            --fq_petvet_dataset ${project}.${dataset} \
            --fq_temp_table_dataset ${project}.temp_tables \
            --fq_destination_dataset ${project}.${dataset} \
            --destination_table ${destination_cohort_table_name} \
            --fq_cohort_sample_names ${fq_cohort_sample_table} \
            --query_project ${project} \
            --fq_sample_mapping_table ${project}.${dataset}.sample_info
    >>>

    runtime {
        docker: docker
        memory: "10 GB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 15
        preemptible: 0
        cpu: 1
    }

 }



