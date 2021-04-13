version 1.0

workflow NgsCohortExtract {
   input {

        File wgs_intervals
        Int scatter_count
       
        File reference
        File reference_index
        File reference_dict
        
        String fq_sample_table
        String fq_cohort_extract_table
        String query_project
        String? fq_filter_set_info_table
        String? fq_filter_set_site_table
        String? fq_filter_set_tranches_table
        String? filter_set_name
        File? excluded_intervals
        Boolean? emit_pls = false
    
        String output_file_base_name
        File? gatk_override
    }

    call SplitIntervals {
      input:
          intervals = wgs_intervals,
          ref_fasta = reference,
          ref_fai = reference_index,
          ref_dict = reference_dict,
          scatter_count = scatter_count
    }
    
    scatter(i in range(scatter_count) ) {
        call ExtractTask {
            input:
                gatk_override            = gatk_override,
                reference                = reference,
                reference_index          = reference_index,
                reference_dict           = reference_dict,
                fq_sample_table          = fq_sample_table,
                intervals                = SplitIntervals.interval_files[i],
                fq_cohort_extract_table  = fq_cohort_extract_table,
                read_project_id          = query_project,
                fq_filter_set_info_table = fq_filter_set_info_table,
                fq_filter_set_site_table = fq_filter_set_site_table,
                fq_filter_set_tranches_table = fq_filter_set_tranches_table,
                filter_set_name          = filter_set_name,
                excluded_intervals       = excluded_intervals,
                emit_pls                 = emit_pls,
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

        File intervals

        String fq_cohort_extract_table
        String read_project_id
        String output_file
        String? fq_filter_set_info_table
        String? fq_filter_set_site_table
        String? fq_filter_set_tranches_table
        String? filter_set_name
        File? excluded_intervals
        
        Boolean? emit_pls

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

        gatk --java-options "-Xmx9g" \
            ExtractCohort \
                --mode GENOMES --ref-version 38 --query-mode LOCAL_SORT \
                -R "~{reference}" \
                -O local.vcf.gz \
                --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
                --sample-table ~{fq_sample_table} \
                --cohort-extract-table ~{fq_cohort_extract_table} \
                -L ~{intervals} \
                --project-id ~{read_project_id} \
                ~{true='--emit-pls' false='' emit_pls} \
                ~{"--filter-set-info-table " + fq_filter_set_info_table} \
                ~{"--filter-set-site-table " + fq_filter_set_site_table} \
                ~{"--tranches-table " + fq_filter_set_tranches_table} \
                ~{"--filter-set-name " + filter_set_name}
        
        # XL does not currently work with ExtractCohort
        gatk --java-options "-Xmx9g" \
            SelectVariants \
            ~{"-XL " + excluded_intervals} \
            -V local.vcf.gz -O ~{output_file}
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "10 GB"
        disks: "local-disk 100 HDD"
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

 task SplitIntervals {
    input {
        File intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count
        String? split_intervals_extra_args

        File? gatk_override
    }

    parameter_meta {
        intervals: {
            localization_optional: true
        }
        ref_fasta: {
            localization_optional: true
        }
        ref_fai: {
            localization_optional: true
        }
        ref_dict: {
            localization_optional: true
        }
     }
	
     command {
         set -e
         export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

         mkdir interval-files
         gatk --java-options "-Xmx5g" SplitIntervals \
             -R ~{ref_fasta} \
             ~{"-L " + intervals} \
             -scatter ~{scatter_count} \
             -O interval-files \
             ~{split_intervals_extra_args}
         cp interval-files/*.interval_list .
     }

     runtime {
         docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
         bootDiskSizeGb: 15
         memory: "6 GB"
         disks: "local-disk 10 HDD"
         preemptible: 3
         cpu: 1
     }

     output {
         Array[File] interval_files = glob("*.interval_list")
     }
 }



