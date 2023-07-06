version 1.0

import "GvsUtils.wdl" as Utils

workflow JointVcfFiltering {
  input {
    String dataset_name
    String project_id
    String base_name

    String filter_set_name
    String filter_set_info_schema
    String fq_filter_set_info_destination_table

    Int num_samples_loaded
    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_idx
    Array[File] sites_only_variant_filtered_vcfs
    Array[File] sites_only_variant_filtered_vcf_idxs

    File? gatk_override

    Int? INDEL_VQSR_max_gaussians_override = 4
    Int? INDEL_VQSR_maximum_training_variants
    Int? INDEL_VQSR_mem_gb_override
    Int? SNP_VQSR_max_gaussians_override = 6
    Int? SNP_VQSR_mem_gb_override
    Int? SNP_VQSR_sample_every_nth_variant
    Int? SNP_VQSR_maximum_training_variants

    # This is the minimum number of samples where the SNP model will be created and applied in separate tasks
    # (SNPsVariantRecalibratorClassic vs. SNPsVariantRecalibratorCreateModel and SNPsVariantRecalibratorScattered)
    # For VQSR classic this is done with 20k but the 10K Stroke Anderson dataset would not work unscattered (at least
    # with the default VM memory settings) so this was adjusted down to 5K.
    Int snps_variant_recalibration_threshold = 5000
  }

  String fq_tranches_destination_table = "~{project_id}.~{dataset_name}.filter_set_tranches"

  Array[String] indel_recalibration_annotations = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
  Array[String] snp_recalibration_annotations = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]

  Array[String] snp_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]
  Array[String] indel_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]

  # reference files
  # Axiom - Used only for indels
  File axiomPoly_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  File axiomPoly_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"

  # DbSNP - BOTH SNPs and INDELs.
  File dbsnp_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  File dbsnp_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

  # HapMap - SNPs
  File hapmap_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"
  File hapmap_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi"

  # Mills - Indels
  File mills_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  File mills_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

  # Omni - SNPs
  File omni_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
  File omni_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"

  # 1000G - SNPs
  File one_thousand_genomes_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  File one_thousand_genomes_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"

  call IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcf,
      sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idx,
      recalibration_filename = base_name + ".indels.recal",
      tranches_filename = base_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotations,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_vcf,
      dbsnp_resource_vcf_index = dbsnp_vcf_index,
      use_allele_specific_annotations = true,
      disk_size = "1000",
      machine_mem_gb = INDEL_VQSR_mem_gb_override,
      max_gaussians = INDEL_VQSR_max_gaussians_override,
      maximum_training_variants = INDEL_VQSR_maximum_training_variants,
  }

  if (num_samples_loaded > snps_variant_recalibration_threshold) {
    call SNPsVariantRecalibratorCreateModel {
      input:
        sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcf,
        sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idx,
        recalibration_filename = base_name + ".snps.recal",
        tranches_filename = base_name + ".snps.tranches",
        model_report_filename = base_name + ".snps.model.report",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotations,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_vcf,
        dbsnp_resource_vcf_index = dbsnp_vcf_index,
        use_allele_specific_annotations = true,
        disk_size = "1000",
        machine_mem_gb = SNP_VQSR_mem_gb_override,
        max_gaussians = SNP_VQSR_max_gaussians_override,
        sample_every_nth_variant = SNP_VQSR_sample_every_nth_variant,
        maximum_training_variants = SNP_VQSR_maximum_training_variants
    }

    scatter (idx in range(length(sites_only_variant_filtered_vcfs))) {
      call SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
        input:
          sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcfs[idx],
          sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idxs[idx],
          recalibration_filename = base_name + ".snps." + idx + ".recal",
          tranches_filename = base_name + ".snps." + idx + ".tranches",
          model_report = SNPsVariantRecalibratorCreateModel.model_report,
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotations,
          hapmap_resource_vcf = hapmap_resource_vcf,
          hapmap_resource_vcf_index = hapmap_resource_vcf_index,
          omni_resource_vcf = omni_resource_vcf,
          omni_resource_vcf_index = omni_resource_vcf_index,
          one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
          one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
          dbsnp_resource_vcf = dbsnp_vcf,
          dbsnp_resource_vcf_index = dbsnp_vcf_index,
          use_allele_specific_annotations = true,
          disk_size = "1000",
          machine_mem_gb = SNP_VQSR_mem_gb_override
      }
    }

    call GatherTranches as SNPGatherTranches {
      input:
        tranches = SNPsVariantRecalibratorScattered.tranches,
        output_filename = base_name + ".snps.gathered.tranches",
        output_tranche_values = snp_recalibration_tranche_values,
        mode = "SNP",
        disk_size = "200",
        gatk_override = gatk_override
    }

    call Utils.MergeVCFs as MergeRecalibrationFiles {
      input:
        input_vcfs = SNPsVariantRecalibratorScattered.recalibration,
        gather_type = "CONVENTIONAL",
        output_vcf_name = "${base_name}.vrecalibration.vcf.gz",
        preemptible_tries = 3,
    }
  }

  if (num_samples_loaded <= snps_variant_recalibration_threshold) {
    call SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
        sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcf,
        sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idx,
        recalibration_filename = base_name + ".snps.recal",
        tranches_filename = base_name + ".snps.tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotations,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_vcf,
        dbsnp_resource_vcf_index = dbsnp_vcf_index,
        use_allele_specific_annotations = true,
        disk_size = "1000",
        machine_mem_gb = SNP_VQSR_mem_gb_override,
        max_gaussians = SNP_VQSR_max_gaussians_override,
    }
  }

  call Utils.PopulateFilterSetInfo {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      filter_schema = filter_set_info_schema,
      fq_filter_set_info_destination_table = fq_filter_set_info_destination_table,
      snp_recal_file = select_first([MergeRecalibrationFiles.output_vcf, SNPsVariantRecalibratorClassic.recalibration]),
      snp_recal_file_index = select_first([MergeRecalibrationFiles.output_vcf_index, SNPsVariantRecalibratorClassic.recalibration_index]),
      indel_recal_file = IndelsVariantRecalibrator.recalibration,
      indel_recal_file_index = IndelsVariantRecalibrator.recalibration_index,
      project_id = project_id,
      useClassic = true
  }

  call PopulateFilterSetTranches {
    input:
      project_id = project_id,
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      snp_recal_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
      indel_recal_tranches = IndelsVariantRecalibrator.tranches,
      fq_tranches_destination_table = fq_tranches_destination_table
  }

  output {
    File snps_variant_recalibration_file = select_first([MergeRecalibrationFiles.output_vcf, SNPsVariantRecalibratorClassic.recalibration])
    File snps_variant_recalibration_file_index = select_first([MergeRecalibrationFiles.output_vcf_index, SNPsVariantRecalibratorClassic.recalibration_index])
    File indels_variant_recalibration_file = IndelsVariantRecalibrator.recalibration
    File indels_variant_recalibration_file_index = IndelsVariantRecalibrator.recalibration_index
    Array[File] monitoring_logs = select_all(
                                  flatten(
                                   [
                                   [IndelsVariantRecalibrator.monitoring_log],
                                   [SNPsVariantRecalibratorCreateModel.monitoring_log],
                                   select_first([SNPsVariantRecalibratorScattered.monitoring_log, []]),
                                   [SNPGatherTranches.monitoring_log],
                                   [MergeRecalibrationFiles.monitoring_log],
                                   [SNPsVariantRecalibratorClassic.monitoring_log],
                                   [PopulateFilterSetInfo.monitoring_log],
                                   [PopulateFilterSetTranches.monitoring_log],
                                   ]))
  }
}

task SNPsVariantRecalibratorCreateModel {
  input {
    String recalibration_filename
    String tranches_filename
    String model_report_filename

    Array[String] recalibration_tranche_values
    Array[String] recalibration_annotation_values

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

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
    Int sample_every_nth_variant = 1
    Int maximum_training_variants = 2500000
    Int? machine_mem_gb

    Int disk_size
  }

  Int machine_mem = select_first([machine_mem_gb, 100])
  Int java_mem = machine_mem - 5

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -euo pipefail

    bash ~{monitoring_script} > monitoring.log &

    gatk --java-options -Xms~{java_mem}g \
    VariantRecalibrator \
    -V ~{sites_only_variant_filtered_vcf} \
    -O ~{recalibration_filename} \
    --tranches-file ~{tranches_filename} \
    --trust-all-polymorphic \
    -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
    -an ~{sep=' -an ' recalibration_annotation_values} \
    ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
    -mode SNP \
    --sample-every-Nth-variant ~{sample_every_nth_variant} \
    --maximum-training-variants ~{maximum_training_variants} \
    --output-model ~{model_report_filename} \
    --max-gaussians ~{max_gaussians} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
    -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
    -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "~{machine_mem} GiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
  }

  output {
    File model_report = "~{model_report_filename}"
    File monitoring_log = "monitoring.log"
  }
}

task GatherTranches {
  input {
    Array[File] tranches
    Array[String] output_tranche_values
    String output_filename
    String mode
    Int disk_size
    File? gatk_override
  }

  parameter_meta {
    tranches: {
                localization_optional: true
              }
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -euo pipefail

    bash ~{monitoring_script} > monitoring.log &

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    tranches_fofn=~{write_lines(tranches)}

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf tranches
    mkdir tranches
    RETRY_LIMIT=5

    count=0
    until cat $tranches_fofn | gsutil -m cp -L cp.log -c -I tranches/; do
    sleep 1
    ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
    echo 'Could not copy all the tranches from the cloud' && exit 1
    fi

    cat $tranches_fofn | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list

    gatk --java-options -Xms6g \
    GatherTranches \
    --input inputs.list \
    --mode ~{mode} \
    -tranche ~{sep=' -tranche ' output_tranche_values} \
    --output ~{output_filename}
  >>>

  runtime {
    memory: "7.5 GiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_07_02_e90d90e00615dcbd9a71d4301fdc04fe2fe155fc"
  }

  output {
    File tranches_file = "~{output_filename}"
    File monitoring_log = "monitoring.log"
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

    File mills_resource_vcf
    File axiomPoly_resource_vcf
    File dbsnp_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf_index
    Boolean use_allele_specific_annotations
    Int max_gaussians = 4
    Int maximum_training_variants = 2500000

    Int disk_size
    Int? machine_mem_gb
  }

  Int machine_mem = select_first([machine_mem_gb, 30])
  Int java_mem = machine_mem - 5

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -euo pipefail

    bash ~{monitoring_script} > monitoring.log &

    gatk --java-options -Xmx~{java_mem}g \
    VariantRecalibrator \
    -V ~{sites_only_variant_filtered_vcf} \
    -O ~{recalibration_filename} \
    --output-model indels.model \
    --tranches-file ~{tranches_filename} \
    --trust-all-polymorphic \
    -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
    -an ~{sep=' -an ' recalibration_annotation_values} \
    ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
    -mode INDEL \
    ~{"--input-model " + model_report} \
    --max-gaussians ~{max_gaussians} \
    --maximum-training-variants ~{maximum_training_variants} \
    -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
    -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "~{machine_mem} GiB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
    File model = "indels.model"
    File monitoring_log = "monitoring.log"
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

    Int disk_size
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
  String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -euo pipefail

    bash ~{monitoring_script} > monitoring.log &

    MODEL_REPORT=~{model_report}

    gatk --java-options -Xmx~{java_mem}g \
    VariantRecalibrator \
    -V ~{sites_only_variant_filtered_vcf} \
    -O ~{recalibration_filename} \
    --output-model snps.model \
    --tranches-file ~{tranches_filename} \
    --trust-all-polymorphic \
    -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
    -an ~{sep=' -an ' recalibration_annotation_values} \
    ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
    -mode SNP \
    --sample-every-Nth-variant 1 \
    ~{model_report_arg} \
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
    docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
    File model = "snps.model"
    File monitoring_log = "monitoring.log"
  }
}

task PopulateFilterSetTranches {
  input {
    String project_id

    File? gatk_override

    String filter_set_name
    String fq_tranches_destination_table

    File snp_recal_tranches
    File indel_recal_tranches
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -eo pipefail

    bash ~{monitoring_script} > monitoring.log &

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    cat ~{snp_recal_tranches} ~{indel_recal_tranches} | grep -v targetTruthSensitivity | grep -v "#" | awk -v CALLSET=~{filter_set_name} '{ print CALLSET "," $0 }' > ~{filter_set_name}.tranches_load.csv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_tranches_destination_table} | sed s/\\./:/)

    echo "Loading combined tranches CSV into ~{fq_tranches_destination_table}"
    bq --apilog=false load --project_id=~{project_id} --skip_leading_rows 0 -F "," \
    --schema "filter_set_name:string,target_truth_sensitivity:float,num_known:integer,num_novel:integer,known_ti_tv:float,novel_ti_tv:float,min_vqslod:float,filter_name:string,model:string,accessible_truth_sites:integer,calls_at_truth_sites:integer,truth_sensitivity:float" \
    ${bq_table} \
    ~{filter_set_name}.tranches_load.csv > status_load_filter_set_tranches
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_07_02_e90d90e00615dcbd9a71d4301fdc04fe2fe155fc"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_tranches = read_string("status_load_filter_set_tranches")
    File monitoring_log = "monitoring.log"
  }
}

