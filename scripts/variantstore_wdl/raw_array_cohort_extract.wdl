version 1.0

workflow RawArrayCohortExtract {
   input {
        Int number_of_partitions = 2
        Int max_probe_id = 1914822
        
        Int probes_per_partition = ceil ( max_probe_id / number_of_partitions)
        
        File reference
        File reference_index
        File reference_dict
    
        String? fq_probe_info_table
        File? probe_info_file
        
        String fq_dataset
        Int max_tables
        String fq_destination_dataset
        String query_project
        String fq_cohort_mapping_table
        File cohort_sample_names_file
        Int ttl = 24
        
        String output_file_base_name
        String? gatk_override
    }
    
    call CreateExtractTable {
        input:
            fq_dataset                = fq_dataset,
            max_tables                = max_tables,
            fq_destination_dataset    = fq_destination_dataset,
            query_project             = query_project,
            fq_sample_mapping_table   = fq_cohort_mapping_table,
            cohort_sample_names_file  = cohort_sample_names_file,
            ttl                       = ttl,
            number_of_partitions      = number_of_partitions,
            probes_per_partition      = probes_per_partition
    }
  
    
    scatter(i in range(number_of_partitions)) {
        call ExtractTask {
            input:
                gatk_override         = gatk_override,
                reference             = reference,
                reference_index       = reference_index,
                reference_dict        = reference_dict,
                fq_probe_info_table   = fq_probe_info_table,
                probe_info_file       = probe_info_file,
                min_probe_id          = 1 + i * probes_per_partition,
                max_probe_id          = (i+1) * probes_per_partition,
                fq_cohort_mapping_table     = fq_cohort_mapping_table,
                cohort_extract_table  = CreateExtractTable.cohort_extract_table,
                read_project_id       = query_project,
                output_file           = "${output_file_base_name}_${i}.vcf.gz"
        }
    }

    call MergeVCFs { 
       input:
           input_vcfs = ExtractTask.output_vcf,
           input_vcfs_indexes = ExtractTask.output_vcf_index,
           output_vcf_name = "${output_file_base_name}.vcf.gz",
           preemptible_tries = 3
    }
    
    output {
        File output_vcf = MergeVCFs.output_vcf
        File output_vcf_idx = MergeVCFs.output_vcf_index
    }
}

################################################################################
task CreateExtractTable {
    # indicates that this task should NOT be call cached
    meta {
       volatile: true
    }

    # ------------------------------------------------
    # Input args:
    input {
        String fq_dataset
        Int max_tables
        String fq_destination_dataset
        String query_project
        String fq_sample_mapping_table
        File cohort_sample_names_file
        Int ttl
        Int number_of_partitions
        Int probes_per_partition
    }

    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e

        uuid=$(cat /proc/sys/kernel/random/uuid | sed s/-/_/g)
        export_table="~{fq_destination_dataset}.${uuid}"
        echo "Exporting to ${export_table}"
        
        python /app/raw_array_cohort_extract.py \
          --dataset ~{fq_dataset} \
          --max_tables ~{max_tables} \
          --fq_destination_table ${export_table} \
          --query_project ~{query_project} \
          --fq_sample_mapping_table ~{fq_sample_mapping_table} \
          --cohort_sample_names_file ~{cohort_sample_names_file} \
          --ttl ~{ttl} \
          --number_of_partitions ~{number_of_partitions} \
          --probes_per_partition ~{probes_per_partition}
          
        echo ${export_table} > cohort_extract_table.txt

    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/variantstore-export:091920"
        memory: "3 GB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 15
        preemptible: 3
        cpu: 1
    }

    # Outputs:
    output {
        String cohort_extract_table = read_string("cohort_extract_table.txt")
    }    
}

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
    
        String? fq_probe_info_table 
        File? probe_info_file
        String probe_info_clause = if defined(probe_info_file) then "--probe-info-csv ${probe_info_file}" else "--probe-info-table ${fq_probe_info_table}"

        Int min_probe_id
        Int max_probe_id
        String fq_cohort_mapping_table
        String cohort_extract_table
        String read_project_id
        String output_file
        
        # Runtime Options:
        File? gatk_override
    }


    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        df -h

        gatk --java-options "-Xmx4g" \
            ArrayExtractCohort \
                -R "~{reference}" \
                -O "~{output_file}" \
                ~{probe_info_clause} \
                --read-project-id "~{read_project_id}" \
                --cohort-sample-table "~{fq_cohort_mapping_table}" \
                --use-compressed-data "false" \
                --cohort-extract-table "~{cohort_extract_table}" \
                --local-sort-max-records-in-ram "5000000" \
                --min-probe-id ~{min_probe_id} --max-probe-id ~{max_probe_id}

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
 
 task MergeVCFs {
   meta {
     volatile: true
   }

   input {
     Array[File] input_vcfs
     Array[File] input_vcfs_indexes
     String output_vcf_name
     Int preemptible_tries
   }

   Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

   # Using MergeVcfs instead of GatherVcfs so we can create indices
   # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
   command {
     java -Xms2000m -jar /usr/gitc/picard.jar \
       MergeVcfs \
       INPUT=~{sep=' INPUT=' input_vcfs} \
       OUTPUT=~{output_vcf_name}
   }
   runtime {
     docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
     preemptible: preemptible_tries
     memory: "3 GiB"
     disks: "local-disk ~{disk_size} HDD"
   }
   output {
     File output_vcf = "~{output_vcf_name}"
     File output_vcf_index = "~{output_vcf_name}.tbi"
   }
 }



