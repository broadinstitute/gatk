version 1.0

workflow NgsCohortExtract {
   input {
        Int max_chrom_id = 24
        
        # bug in cromwell, can't support large integers...
        # https://github.com/broadinstitute/cromwell/issues/2685
        String chrom_offset = "1000000000000"
       
        File reference
        File reference_index
        File reference_dict
        
        String fq_sample_table
        String fq_cohort_extract_table
        String query_project
    
        String output_file_base_name
        File? gatk_override
    }
    
    scatter(i in range(max_chrom_id)) {
        call ExtractTask {
            input:
                gatk_override            = gatk_override,
                reference                = reference,
                reference_index          = reference_index,
                reference_dict           = reference_dict,
                fq_sample_table          = fq_sample_table,
                chrom_offset             = chrom_offset,
                chrom_id                 = i+1,
                fq_cohort_extract_table  = fq_cohort_extract_table,
                read_project_id          = query_project,
                output_file              = "${output_file_base_name}_${i}.vcf.gz"
        }
    }
}

################################################################################
task ExtractTask {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    input {
        # ------------------------------------------------
        # Input args:
        File reference
        File reference_index
        File reference_dict
    
        String fq_sample_table

        # bug in cromwell, can't support large integers...
        # https://github.com/broadinstitute/cromwell/issues/2685
        String chrom_offset
        Int chrom_id

        String fq_cohort_extract_table
        String read_project_id
        String output_file
        
        # Runtime Options:
        File? gatk_override
        
        Int? local_sort_max_records_in_ram = 10000000
    }


    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        df -h
        min_location=$(echo "~{chrom_id} * ~{chrom_offset}" | bc)
        max_location=$(echo "( ~{chrom_id} + 1 ) * ~{chrom_offset}" | bc)

        gatk --java-options "-Xmx4g" \
            ExtractCohort \
                --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
                -R "~{reference}" \
                -O "~{output_file}" \
                --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
                --sample-table ~{fq_sample_table} \
                --cohort-extract-table ~{fq_cohort_extract_table} \
                --min-location ${min_location} --max-location ${max_location} \
                --project-id ~{read_project_id}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "7 GB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 15
        preemptible: 3
        cpu: 2
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File output_vcf = "~{output_file}"
        File output_vcf_index = "~{output_file}.tbi"
    }
 }
 



