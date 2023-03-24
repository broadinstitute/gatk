version 1.0

import "GvsWarpTasks.wdl" as Tasks
import "GvsUtils.wdl" as Utils
import "../../vcf_site_level_filtering_wdl/JointVcfFiltering.wdl" as VQSRLite

workflow GvsCreateFilterSet {
  input {
    Boolean go = true
    String dataset_name
    String project_id
    String call_set_identifier

    String filter_set_name
    Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
    Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File? gatk_override

    Boolean use_classic_VQSR = true
    Int? INDEL_VQSR_max_gaussians_override = 4
    Int? INDEL_VQSR_maximum_training_variants
    Int? INDEL_VQSR_mem_gb_override
    Int? SNP_VQSR_max_gaussians_override = 6
    Int? SNP_VQSR_mem_gb_override
    Int? SNP_VQSR_sample_every_nth_variant
    Int? SNP_VQSR_maximum_training_variants
    # This is the minimum number of samples where the SNP model will be created and applied in separate tasks
    # (SNPsVariantRecalibratorClassic vs. SNPsVariantRecalibratorCreateModel and SNPsVariantRecalibratorScattered)
    # For WARP classic this is done with 20k but the 10K Stroke Anderson dataset would not work unscattered (at least
    # with the default VM memory settings) so this was adjusted down to 5K.
    Int snps_variant_recalibration_threshold = 5000
  }

  Array[String] snp_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]

  # reference files
  File axiomPoly_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  File axiomPoly_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"
  File dbsnp_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  File dbsnp_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"
  File hapmap_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"
  File hapmap_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi"
  File mills_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  File mills_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
  File omni_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
  File omni_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"
  File one_thousand_genomes_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  File one_thousand_genomes_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"
  File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  # fully-qualified table names
  String fq_sample_table = "~{project_id}.~{dataset_name}.sample_info"
  String fq_alt_allele_table = "~{project_id}.~{dataset_name}.alt_allele"
  String fq_info_destination_table = "~{project_id}.~{dataset_name}.filter_set_info"
  String fq_info_destination_table_vqsr_lite = "~{project_id}.~{dataset_name}.filter_set_info_vqsr_lite"
  String fq_tranches_destination_table = "~{project_id}.~{dataset_name}.filter_set_tranches"
  String fq_filter_sites_destination_table = "~{project_id}.~{dataset_name}.filter_set_sites"

  String fq_info_destination_table_schema =           "filter_set_name:string,type:string,location:integer,ref:string,alt:string,vqslod:float,culprit:string,training_label:string,yng_status:string"
  String fq_info_destination_table_vqsr_lite_schema = "filter_set_name:string,type:string,location:integer,ref:string,alt:string,calibration_sensitivity:float,score:float,training_label:string,yng_status:string"

  call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
    input:
      project_id = project_id,
      fq_table = fq_sample_table
  }

  call Utils.GetNumSamplesLoaded {
    input:
      fq_sample_table = fq_sample_table,
      project_id = project_id,
      sample_table_timestamp = SamplesTableDatetimeCheck.last_modified_timestamp,
  }

  Int scatter_count = if GetNumSamplesLoaded.num_samples < 100 then 20
                      else if GetNumSamplesLoaded.num_samples < 1000 then 100
                           else if GetNumSamplesLoaded.num_samples < 10000 then 200
                                else if GetNumSamplesLoaded.num_samples < 100000 then 500 else 1000

  call Utils.SplitIntervals {
    input:
      intervals = interval_list,
      ref_fasta = reference,
      ref_fai = reference_index,
      ref_dict = reference_dict,
      scatter_count = scatter_count,
      gatk_override = gatk_override
  }

  call Utils.GetBQTableLastModifiedDatetime as AltAlleleTableDatetimeCheck {
    input:
      project_id = project_id,
      fq_table = fq_alt_allele_table
  }

  scatter(i in range(length(SplitIntervals.interval_files))) {
    call ExtractFilterTask {
      input:
        gatk_override              = gatk_override,
        reference                  = reference,
        reference_index            = reference_index,
        reference_dict             = reference_dict,
        fq_sample_table            = fq_sample_table,
        sample_table_timestamp     = SamplesTableDatetimeCheck.last_modified_timestamp,
        intervals                  = SplitIntervals.interval_files[i],
        fq_alt_allele_table        = fq_alt_allele_table,
        alt_allele_table_timestamp = AltAlleleTableDatetimeCheck.last_modified_timestamp,
        excess_alleles_threshold   = 1000000,
        output_file                = "${filter_set_name}_${i}.vcf.gz",
        project_id                 = project_id,
        dataset_id                 = dataset_name,
        call_set_identifier        = call_set_identifier
    }
  }

  call Utils.MergeVCFs as MergeVCFs {
    input:
      input_vcfs = ExtractFilterTask.output_vcf,
      output_vcf_name = "${filter_set_name}.vcf.gz",
      preemptible_tries = 3,
  }

  # From this point, the paths diverge depending on whether they're using classic VQSR or VQSR-Lite
  # The first branch here is VQSR-Lite, and the second is classic VQSR
  if (!use_classic_VQSR) {
    call VQSRLite.JointVcfFiltering as JointVcfFiltering {
      input:
        vcf = ExtractFilterTask.output_vcf,
        vcf_index = ExtractFilterTask.output_vcf_index,
        sites_only_vcf = MergeVCFs.output_vcf,
        sites_only_vcf_index = MergeVCFs.output_vcf_index,
        basename = filter_set_name,
        gatk_docker = "us.gcr.io/broad-gatk/gatk:4.3.0.0",
        extract_interval_list = interval_list,
        score_interval_list = interval_list,
        snp_annotations   = "-A AS_QD -A AS_MQRankSum -A AS_ReadPosRankSum -A AS_FS -A AS_MQ -A AS_SOR",
        indel_annotations = "-A AS_QD -A AS_MQRankSum -A AS_ReadPosRankSum -A AS_FS -A AS_MQ -A AS_SOR",
        use_allele_specific_annotations = true,
        monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    }

    call Utils.MergeVCFs as MergeINDELScoredVCFs {
      input:
        input_vcfs = JointVcfFiltering.indels_variant_scored_vcf,
        gather_type = "CONVENTIONAL",
        output_vcf_name = "${filter_set_name}.indel.vrecalibration.gz",
        preemptible_tries = 3,
    }

    call Utils.MergeVCFs as MergeSNPScoredVCFs {
      input:
        input_vcfs = JointVcfFiltering.snps_variant_scored_vcf,
        gather_type = "CONVENTIONAL",
        output_vcf_name = "${filter_set_name}.snp.vrecalibration.gz",
        preemptible_tries = 3,
    }

    # These calls to SelectVariants are being added for two reasons
    # 1) The snps_variant_scored_vcf and indels_variant_scored_vcf output by JointVcfFiltering contains ALL variants,
    #     but are currently ONLY annotating SNPs and INDELs respectively.
    # 2) Those output VCFs also contain filtered sites (sites at which the FILTER field set to anything other than '.' or 'PASS')
    #     which we don't want to put into the filter_set_info_vqsr_lite table.
    call Utils.SelectVariants as CreateFilteredScoredSNPsVCF {
      input:
        input_vcf = MergeSNPScoredVCFs.output_vcf,
        input_vcf_index = MergeSNPScoredVCFs.output_vcf_index,
        type_to_include = "SNP",
        exclude_filtered = true,
        output_basename = "${filter_set_name}.filtered.scored.snps"
    }

    call Utils.SelectVariants as CreateFilteredScoredINDELsVCF {
      input:
        input_vcf = MergeINDELScoredVCFs.output_vcf,
        input_vcf_index = MergeINDELScoredVCFs.output_vcf_index,
        type_to_include = "INDEL",
        exclude_filtered = true,
        output_basename = "${filter_set_name}.filtered.scored.indels"
    }

    call PopulateFilterSetInfo {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        snp_recal_file = CreateFilteredScoredSNPsVCF.output_vcf,
        snp_recal_file_index = CreateFilteredScoredSNPsVCF.output_vcf_index,
        indel_recal_file = CreateFilteredScoredINDELsVCF.output_vcf,
        indel_recal_file_index = CreateFilteredScoredINDELsVCF.output_vcf_index,
        fq_info_destination_table = fq_info_destination_table_vqsr_lite,
        filter_schema = fq_info_destination_table_vqsr_lite_schema,
        project_id = project_id,
        useClassic = false
    }

    call PopulateFilterSetSites {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
        sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
        fq_filter_sites_destination_table = fq_filter_sites_destination_table,
        project_id = project_id
    }

    call UberMonitor as UberMonitorLite {
      input:
        inputs = flatten([[JointVcfFiltering.extract_variant_anotations_snps_monitoring_log],
                  [JointVcfFiltering.extract_variant_anotations_indels_monitoring_log],
                  [JointVcfFiltering.train_variant_anotation_model_snps_monitoring_log],
                  [JointVcfFiltering.train_variant_anotation_model_indels_monitoring_log],
                  JointVcfFiltering.score_variant_annotations_snps_monitoring_log,
                  JointVcfFiltering.score_variant_annotations_indels_monitoring_log,
                  [MergeSNPScoredVCFs.monitoring_log],
                  [MergeINDELScoredVCFs.monitoring_log],
                  [CreateFilteredScoredSNPsVCF.monitoring_log],
                  [CreateFilteredScoredINDELsVCF.monitoring_log],
                  [PopulateFilterSetInfo.monitoring_log],
                  [PopulateFilterSetSites.monitoring_log]])
    }
  }

  if (use_classic_VQSR) {

    call Tasks.IndelsVariantRecalibrator {
      input:
        sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
        sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
        recalibration_filename = filter_set_name + ".indels.recal",
        tranches_filename = filter_set_name + ".indels.tranches",
        recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"],
        recalibration_annotation_values = indel_recalibration_annotation_values,
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

    if (GetNumSamplesLoaded.num_samples > snps_variant_recalibration_threshold) {
      call Tasks.SNPsVariantRecalibratorCreateModel {
        input:
          sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
          sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
          recalibration_filename = filter_set_name + ".snps.recal",
          tranches_filename = filter_set_name + ".snps.tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
          model_report_filename = filter_set_name + ".snps.model.report",
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

      scatter (idx in range(length(ExtractFilterTask.output_vcf))) {
        call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
          input:
            sites_only_variant_filtered_vcf = ExtractFilterTask.output_vcf[idx],
            sites_only_variant_filtered_vcf_index = ExtractFilterTask.output_vcf_index[idx],
            recalibration_filename = filter_set_name + ".snps." + idx + ".recal",
            tranches_filename = filter_set_name + ".snps." + idx + ".tranches",
            recalibration_tranche_values = snp_recalibration_tranche_values,
            recalibration_annotation_values = snp_recalibration_annotation_values,
            model_report = SNPsVariantRecalibratorCreateModel.model_report,
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

      call Tasks.GatherTranches as SNPGatherTranches {
        input:
          tranches = SNPsVariantRecalibratorScattered.tranches,
          output_filename = filter_set_name + ".snps.gathered.tranches",
          output_tranche_values = snp_recalibration_tranche_values,
          mode = "SNP",
          disk_size = "200",
          gatk_override = gatk_override
      }

      call Utils.MergeVCFs as MergeRecalibrationFiles {
        input:
          input_vcfs = SNPsVariantRecalibratorScattered.recalibration,
          gather_type = "CONVENTIONAL",
          output_vcf_name = "${filter_set_name}.vrecalibration.gz",
          preemptible_tries = 3,
      }
      Array[Array[File]] m1_logs = flatten([[IndelsVariantRecalibrator.monitoring_log],
                                            [SNPsVariantRecalibratorCreateModel.monitoring_log],
                                            SNPsVariantRecalibratorScattered.monitoring_log,
                                            [SNPGatherTranches.monitoring_log],
                                            [MergeRecalibrationFiles.monitoring_log]])
    }

    if (GetNumSamplesLoaded.num_samples <= snps_variant_recalibration_threshold) {
      call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
        input:
          sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
          sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
          recalibration_filename = filter_set_name + ".snps.recal",
          tranches_filename = filter_set_name + ".snps.tranches",
          recalibration_tranche_values = snp_recalibration_tranche_values,
          recalibration_annotation_values = snp_recalibration_annotation_values,
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
      Array[Array[File]] m1_logs = flatten([[SNPsVariantRecalibratorClassic.monitoring_log]])
    }

    call PopulateFilterSetInfo as PopulateFilterSetInfoClassic {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        snp_recal_file = select_first([MergeRecalibrationFiles.output_vcf, SNPsVariantRecalibratorClassic.recalibration]),
        snp_recal_file_index = select_first([MergeRecalibrationFiles.output_vcf_index, SNPsVariantRecalibratorClassic.recalibration_index]),
        indel_recal_file = IndelsVariantRecalibrator.recalibration,
        indel_recal_file_index = IndelsVariantRecalibrator.recalibration_index,
        fq_info_destination_table = fq_info_destination_table,
        filter_schema = fq_info_destination_table_schema,
        project_id = project_id,
        useClassic = true
    }

    call PopulateFilterSetSites as PopulateFilterSetSitesClassic {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
        sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
        fq_filter_sites_destination_table = fq_filter_sites_destination_table,
        project_id = project_id
    }

    call PopulateFilterSetTranches {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        snp_recal_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
        indel_recal_tranches = IndelsVariantRecalibrator.tranches,
        fq_tranches_destination_table = fq_tranches_destination_table,
        project_id = project_id
    }

    call UberMonitor as UberMonitorClassic {
      input:
        inputs = flatten([m1_logs,
                         [PopulateFilterSetInfoClassic.monitoring_log],
                         [PopulateFilterSetSitesClassic.monitoring_log],
                         [PopulateFilterSetTranches.monitoring_log]])
    }
  }


  output {
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_idx = MergeVCFs.output_vcf_index
    Boolean done = true
  }
}

################################################################################

task ExtractFilterTask {
  input {
    String project_id
    String dataset_id
    String call_set_identifier

    File reference
    File reference_index
    File reference_dict

    String fq_sample_table
    String sample_table_timestamp

    File intervals

    String fq_alt_allele_table
    String alt_allele_table_timestamp

    String cost_observability_tablename = "cost_observability"

    String output_file
    Int? excess_alleles_threshold

    # Runtime Options:
    File? gatk_override
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  String intervals_name = basename(intervals)


  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    df -h

    gatk --java-options "-Xmx4g" ExtractFeatures \
      --ref-version 38  \
      -R "~{reference}" \
      -O "~{output_file}" \
      --local-sort-max-records-in-ram 1000000 \
      --sample-table ~{fq_sample_table} \
      --alt-allele-table ~{fq_alt_allele_table} \
      ~{"--excess-alleles-threshold " + excess_alleles_threshold} \
      -L ~{intervals} \
      --dataset-id ~{dataset_id} \
      --project-id ~{project_id} \
      --cost-observability-tablename ~{cost_observability_tablename} \
      --call-set-identifier ~{call_set_identifier} \
      --wdl-step GvsCreateFilterSet \
      --wdl-call ExtractFilterTask \
      --shard-identifier ~{intervals_name}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_03_01_b01183576153cf000e17dea32144d332cb7b79a9"
    memory: "7 GB"
    disks: "local-disk 10 HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    maxRetries: 3
    cpu: 2
  }

  output {
    File output_vcf = "~{output_file}"
    File output_vcf_index = "~{output_file}.tbi"
  }
}

task PopulateFilterSetInfo {
  input {
    String filter_set_name
    String filter_schema
    String fq_info_destination_table
    Boolean useClassic = true

    File snp_recal_file
    File snp_recal_file_index
    File indel_recal_file
    File indel_recal_file_index

    String project_id

    File? gatk_override
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -eo pipefail

    bash ~{monitoring_script} > monitoring.log &

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    echo "Creating SNPs recalibration file"
    gatk --java-options "-Xmx1g" \
      CreateFilteringFiles \
      --ref-version 38 \
      --filter-set-name ~{filter_set_name} \
      -mode SNP \
      --classic ~{useClassic} \
      -V ~{snp_recal_file} \
      -O ~{filter_set_name}.snps.recal.tsv

    echo "Creating INDELs racalibration file"
    gatk --java-options "-Xmx1g" \
      CreateFilteringFiles \
      --ref-version 38 \
      --filter-set-name ~{filter_set_name} \
      -mode INDEL \
      --classic ~{useClassic} \
      -V ~{indel_recal_file} \
      -O ~{filter_set_name}.indels.recal.tsv

    # merge into a single file
    echo "Merging SNP + INDELs"
    cat ~{filter_set_name}.snps.recal.tsv ~{filter_set_name}.indels.recal.tsv | grep -v filter_set_name | grep -v "#"  > ~{filter_set_name}.filter_set_load.tsv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_info_destination_table} | sed s/\\./:/)

    echo "Loading combined TSV into ~{fq_info_destination_table}"
    bq load --project_id=~{project_id} --skip_leading_rows 0 -F "tab" \
      --range_partitioning=location,0,26000000000000,6500000000 \
      --clustering_fields=location \
      --schema "~{filter_schema}" \
      ${bq_table} \
      ~{filter_set_name}.filter_set_load.tsv > status_load_filter_set_info
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_03_01_b01183576153cf000e17dea32144d332cb7b79a9"
    memory: "3500 MB"
    disks: "local-disk 250 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_info = read_string("status_load_filter_set_info")
    File monitoring_log = "monitoring.log"
  }
}

task PopulateFilterSetSites {
  input {
    String filter_set_name
    String fq_filter_sites_destination_table

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    String project_id

    File? gatk_override
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -eo pipefail

    bash ~{monitoring_script} > monitoring.log &

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    echo "Generating filter set sites TSV"
    gatk --java-options "-Xmx1g" \
    CreateSiteFilteringFiles \
    --ref-version 38 \
    --filter-set-name ~{filter_set_name} \
    -V ~{sites_only_variant_filtered_vcf} \
    -O ~{filter_set_name}.filter_sites_load.tsv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_filter_sites_destination_table} | sed s/\\./:/)

    echo "Loading filter set sites TSV into ~{fq_filter_sites_destination_table}"
    bq load --project_id=~{project_id} --skip_leading_rows 1 -F "tab" \
    --range_partitioning=location,0,26000000000000,6500000000 \
    --clustering_fields=location \
    --schema "filter_set_name:string,location:integer,filters:string" \
    ${bq_table} \
    ~{filter_set_name}.filter_sites_load.tsv > status_load_filter_set_sites
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_03_01_b01183576153cf000e17dea32144d332cb7b79a9"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_sites = read_string("status_load_filter_set_sites")
    File monitoring_log = "monitoring.log"
  }
}

task PopulateFilterSetTranches {
  input {
    File? gatk_override

    String filter_set_name
    String fq_tranches_destination_table

    File snp_recal_tranches
    File indel_recal_tranches

    String project_id
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
    bq load --project_id=~{project_id} --skip_leading_rows 0 -F "," \
    --schema "filter_set_name:string,target_truth_sensitivity:float,num_known:integer,num_novel:integer,known_ti_tv:float,novel_ti_tv:float,min_vqslod:float,filter_name:string,model:string,accessible_truth_sites:integer,calls_at_truth_sites:integer,truth_sensitivity:float" \
    ${bq_table} \
    ~{filter_set_name}.tranches_load.csv > status_load_filter_set_tranches
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_03_01_b01183576153cf000e17dea32144d332cb7b79a9"
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

task UberMonitor {
  input {
    Array[File] inputs
  }

  command <<<
    set -e

    python3 /app/uber_monitor.py \
    --input ~{sep=" " inputs}
  >>>

  # ------------------------------------------------
  # Runtime settings:
  runtime {
    docker: "us.gcr.io/broad-dsde-methods/variantstore:ah_var_store_2023_03_24"
    memory: "1 GB"
    preemptible: 3
    cpu: "1"
    disks: "local-disk 100 HDD"
  }
  # ------------------------------------------------
  # Outputs:
  output {
    File out = stdout()
  }
}
