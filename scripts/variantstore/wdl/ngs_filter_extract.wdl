version 1.0

workflow NgsFilterExtract {
   input {
        Int max_chrom_id = 24
        
        # bug in cromwell, can't support large integers...
        # https://github.com/broadinstitute/cromwell/issues/2685
        String chrom_offset = "1000000000000"
       
        File reference
        File reference_index
        File reference_dict
        
        String fq_sample_table
        String fq_alt_allele_table
        String query_project
    
        String filter_set_name
        String fq_filter_set_info_table
        String fq_filter_set_tranches_table

        String output_file_base_name
        File? gatk_override

        File dbsnp_vcf
        File dbsnp_vcf_index

        File? snps_model
        File? indels_model

        Array[String] snp_recalibration_tranche_values
        Array[String] snp_recalibration_annotation_values
        Array[String] indel_recalibration_tranche_values
        Array[String] indel_recalibration_annotation_values

        File? excluded_sites_bed

        File hapmap_resource_vcf
        File hapmap_resource_vcf_index
        File omni_resource_vcf
        File omni_resource_vcf_index
        File one_thousand_genomes_resource_vcf
        File one_thousand_genomes_resource_vcf_index
        File mills_resource_vcf
        File mills_resource_vcf_index
        File axiomPoly_resource_vcf
        File axiomPoly_resource_vcf_index
        File dbsnp_resource_vcf = dbsnp_vcf
        File dbsnp_resource_vcf_index = dbsnp_vcf_index

        # Runtime attributes
        Int? small_disk_override
        Int small_disk = select_first([small_disk_override, "100"])
        Int? medium_disk_override
        Int medium_disk = select_first([medium_disk_override, "200"])
        Int? large_disk_override
        Int large_disk = select_first([large_disk_override, "300"])
        Int? huge_disk_override
        Int huge_disk = select_first([huge_disk_override, "400"])

        String? preemptible_tries_override
        Int preemptible_tries = select_first([preemptible_tries_override, "3"])

        Int? SNP_VQSR_machine_mem_gb
        Int? SNP_VQSR_downsampleFactor = 1
    }
    
    scatter(i in range(max_chrom_id)) {
        call ExtractFilterTask {
            input:
                gatk_override            = gatk_override,
                reference                = reference,
                reference_index          = reference_index,
                reference_dict           = reference_dict,
                fq_sample_table          = fq_sample_table,
                chrom_offset             = chrom_offset,
                chrom_id                 = i+1,
                fq_alt_allele_table      = fq_alt_allele_table,
                read_project_id          = query_project,
                output_file              = "${output_file_base_name}_${i}.vcf.gz"
        }
    }

    call MergeVCFs { 
       input:
           input_vcfs = ExtractFilterTask.output_vcf,
           input_vcfs_indexes = ExtractFilterTask.output_vcf_index,
           output_vcf_name = "${output_file_base_name}.vcf.gz",
           preemptible_tries = 3
    }
    
    call IndelsVariantRecalibrator {
        input:
        sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
        sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
        model_report = indels_model,
        recalibration_filename = filter_set_name + ".indels.recal",
        tranches_filename = filter_set_name + ".indels.tranches",
        recalibration_tranche_values = indel_recalibration_tranche_values,
        recalibration_annotation_values = indel_recalibration_annotation_values,
        excluded_sites_bed = excluded_sites_bed,
        mills_resource_vcf = mills_resource_vcf,
        mills_resource_vcf_index = mills_resource_vcf_index,
        axiomPoly_resource_vcf = axiomPoly_resource_vcf,
        axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = true,
        disk_size = large_disk
    }

    call SNPsVariantRecalibrator {
      input:
          sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
          sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
          model_report = snps_model,
          recalibration_filename = filter_set_name + ".snps.recal",
          tranches_filename = filter_set_name + ".snps.tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          excluded_sites_bed = excluded_sites_bed,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_resource_vcf,
          dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
          use_allele_specific_annotations = true,
          disk_size = large_disk,
          machine_mem_gb = SNP_VQSR_machine_mem_gb,
          downsampleFactor= SNP_VQSR_downsampleFactor
    }

   call UploadFilterSetToBQ {
     input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        snp_recal_file = SNPsVariantRecalibrator.recalibration,
        snp_recal_file_index = SNPsVariantRecalibrator.recalibration_index,
        snp_recal_tranches = SNPsVariantRecalibrator.tranches,

        indel_recal_file = IndelsVariantRecalibrator.recalibration,
        indel_recal_file_index = IndelsVariantRecalibrator.recalibration_index,
        indel_recal_tranches = IndelsVariantRecalibrator.tranches,

        fq_info_destination_table = fq_filter_set_info_table,
        fq_tranches_destination_table = fq_filter_set_tranches_table
   }

    output {
        File output_vcf = MergeVCFs.output_vcf
        File output_vcf_idx = MergeVCFs.output_vcf_index
    }
}

################################################################################
task ExtractFilterTask {
    # indicates that this task should NOT be call cached
    
    # TODO: should this be marked as volatile???
    #meta {
    #   volatile: true
    #}

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

        String fq_alt_allele_table
        String read_project_id
        String output_file
        
        # Runtime Options:
        File? gatk_override
        
        Int? local_sort_max_records_in_ram = 1000000
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
            ExtractFeatures \
                --ref-version 38  \
                -R "~{reference}" \
                -O "~{output_file}" \
                --local-sort-max-records-in-ram ~{local_sort_max_records_in_ram} \
                --sample-table ~{fq_sample_table} \
                --alt-allele-table ~{fq_alt_allele_table} \
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

task UploadFilterSetToBQ {
    # indicates that this task should NOT be call cached
    
    # TODO: should this be marked as volatile???
    #meta {
    #   volatile: true
    #}

    input {
        # ------------------------------------------------
        # Input args:

        # Runtime Options:
        File? gatk_override

        String filter_set_name

        File snp_recal_file
        File snp_recal_file_index
        File snp_recal_tranches

        File indel_recal_file
        File indel_recal_file_index
        File indel_recal_tranches

        String fq_info_destination_table
        String fq_tranches_destination_table
    }


    # ------------------------------------------------
    # Run our command:
    command <<<
        set -e
        set -o pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx4g" \
            CreateFilteringFiles \
            --ref-version 38 \
            --filter-set-name ~{filter_set_name} \
            -mode SNP \
            -V ~{snp_recal_file} \
            -O ~{filter_set_name}.snps.recal.tsv

        gatk --java-options "-Xmx4g" \
            CreateFilteringFiles \
            --ref-version 38 \
            --filter-set-name ~{filter_set_name} \
            -mode INDEL \
            -V ~{indel_recal_file} \
            -O ~{filter_set_name}.indels.recal.tsv

        # merge into a single file
        echo "Merging SNP + INDELs"
        cat ~{filter_set_name}.snps.recal.tsv ~{filter_set_name}.indels.recal.tsv | grep -v filter_set_name | grep -v "#"  > filter_set_load.tsv

        # BQ load likes a : instead of a . after the project 
        bq_info_table=$(echo ~{fq_info_destination_table} | sed s/\\./:/)

        echo "Loading Filter Set into BQ"
        bq load --skip_leading_rows 0 -F "tab" \
        ${bq_info_table} \
        filter_set_load.tsv \
        "filter_set_name:string,type:string,location:integer,ref:string,alt:string,vqslod:float,culprit:string,training_label:string,yng_status:string"

        echo "Merging Tranches"
        cat ~{snp_recal_tranches} ~{indel_recal_tranches} | grep -v targetTruthSensitivity | grep -v "#" | awk -v CALLSET=~{filter_set_name} '{ print CALLSET "," $0 }' > tranches_load.csv

        # BQ load likes a : instead of a . after the project 
        bq_tranches_table=$(echo ~{fq_tranches_destination_table} | sed s/\\./:/)

        echo "Loading Tranches into BQ"
        bq load --skip_leading_rows 0 -F "," \
        ${bq_tranches_table} \
        tranches_load.csv \
        "filter_set_name:string,target_truth_sensitivity:float,num_known:integer,num_novel:integer,known_ti_tv:float,novel_ti_tv:float,min_vqslod:float,filter_name:string,model:string,accessible_truth_sites:integer,calls_at_truth_sites:integer,truth_sensitivity:float"
    >>>

    # ------------------------------------------------
    # Runtime settings:
    runtime {
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
        memory: "7 GB"
        disks: "local-disk 200 HDD"
        bootDiskSizeGb: 15
        preemptible: 3
        cpu: 2
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File filter_set_load_tsv = "filter_set_load.tsv"
        File filter_tranche_load_tsv = "tranches_load.csv"
    }
 }

 task MergeVCFs {
    # TODO: should this be marked as volatile???
    #meta {
    #   volatile: true
    #}
 
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

task IndelsVariantRecalibrator {

  input {
    String recalibration_filename
    String tranches_filename
    File? model_report

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File? excluded_sites_bed
    File mills_resource_vcf
    File axiomPoly_resource_vcf
    File dbsnp_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 4

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms24g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      ~{"-XL " + excluded_sites_bed} \
      -O ~{recalibration_filename} \
      --output-model indels.model \
      --rscript-file indels.Rscript \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode INDEL \
      ~{"--input-model " + model_report} \
      --max-gaussians ~{max_gaussians} \
      -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "26 GiB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
    File model = "indels.model"
    File rscript = "indels.Rscript"
  }
}

task SNPsVariantRecalibrator {

  input {
    String recalibration_filename
    String tranches_filename
    File? model_report

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    File? excluded_sites_bed
    File hapmap_resource_vcf
    File omni_resource_vcf
    File one_thousand_genomes_resource_vcf
    File dbsnp_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 6
    Int? downsampleFactor = 1

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    Int? machine_mem_gb

  }

  Int auto_mem = ceil(2 * size([sites_only_variant_filtered_vcf,
                                hapmap_resource_vcf,
                                omni_resource_vcf,
                                one_thousand_genomes_resource_vcf,
                                dbsnp_resource_vcf],
                     "GiB"))
  Int machine_mem = select_first([machine_mem_gb, if auto_mem < 7 then 7 else auto_mem])
  Int java_mem = machine_mem - 5

  command <<<
    set -euo pipefail

    gatk --java-options -Xmx~{java_mem}g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      ~{"-XL " + excluded_sites_bed} \
      -O ~{recalibration_filename} \
      --output-model snps.model \
      --rscript-file snps.Rscript \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode SNP \
      --sample-every-Nth-variant ~{downsampleFactor} \
      ~{"--input-model " + model_report} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "~{machine_mem} GiB"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 0
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
    File model = "snps.model"
    File rscript = "snps.Rscript"
  }
}

 



