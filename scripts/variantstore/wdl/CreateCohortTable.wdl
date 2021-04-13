version 1.0

workflow CreateCohortTable {
   input {
        String project
        String dataset

        String? docker
    }

    String docker_final = select_first([docker, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])

    call CreateCohortTableTask {
        input:
            project                         = project,
            dataset                         = dataset,

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

        String? query_project
        String? destination_project
        String? destination_dataset

        String? destination_cohort_table_name
        String? fq_cohort_sample_table
        String? fq_sample_mapping_table

        String? git_branch_for_script
        String docker
    }

    #### set defaults ####
    String query_project_final = if defined(query_project) then "${query_project}" else "${project}"
    String destination_project_final = if defined(destination_project) then "${destination_project}" else "${project}"
    String destination_dataset_final = if defined(destination_dataset) then "${destination_dataset}" else "${dataset}"

    String destination_cohort_table_name_final = if defined(destination_cohort_table_name) then "${destination_cohort_table_name}" else "exported_cohort_all_samples"
    String fq_cohort_sample_table_final = if defined(fq_cohort_sample_table) then "${fq_cohort_sample_table}" else "${project}.${dataset}.sample_info"
    String fq_sample_mapping_table_final = if defined(fq_sample_mapping_table) then "${fq_sample_mapping_table}" else "${project}.${dataset}.sample_info"

    String git_branch_for_script_final = if defined(git_branch_for_script) then "${git_branch_for_script}" else "master"

    command <<<
        set -e

        # TODO access this from docker directly
        wget -L "https://raw.githubusercontent.com/broadinstitute/gatk/~{git_branch_for_script_final}/scripts/variantstore/wdl/extract/create_cohort_data_table.py"

        python create_cohort_data_table.py \
            --fq_petvet_dataset ~{project}.~{dataset} \
            --fq_temp_table_dataset ~{destination_project_final}.temp_tables \
            --fq_destination_dataset ~{destination_project_final}.~{destination_dataset_final} \
            --destination_table ~{destination_cohort_table_name_final} \
            --fq_cohort_sample_names ~{fq_cohort_sample_table_final} \
            --query_project ~{query_project_final} \
            --fq_sample_mapping_table ~{fq_sample_mapping_table_final}
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



