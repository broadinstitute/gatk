version 1.0

import "GvsWarpTasks.wdl" as Tasks
import "GvsUtils.wdl" as Utils

workflow JointVcfFiltering {
  input {
    String base_name
    Int num_samples_loaded
    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_idx
    Array[File] sites_only_variant_filtered_vcfs
    Array[File] sites_only_variant_filtered_vcf_idxs

    Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
    Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
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

  Array[String] snp_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]

  # reference files
  # Axiom - Used only for indels
  # Classic: known=false,training=true,truth=false
  File axiomPoly_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
  File axiomPoly_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi"

  # DbSNP - BOTH SNPs and INDELs. But used only as known in classic (which isn't used in Lite and so dropped in lite)
  # Classic: known=true,training=false,truth=false
  File dbsnp_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
  File dbsnp_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

  # HapMap - SNPs
  # Classic: known=false,training=true,truth=true
  File hapmap_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz"
  File hapmap_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi"

  # Mills - Indels
  # Classic: known=false,training=true,truth=true
  File mills_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
  File mills_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

  # Omni - SNPs
  # Classic: known=false,training=true,truth=true
  File omni_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz"
  File omni_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi"

  # 1000G - SNPs
  # Classic: known=false,training=true,truth=false
  File one_thousand_genomes_resource_vcf = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
  File one_thousand_genomes_resource_vcf_index = "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi"

  call Tasks.IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcf,
      sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idx,
      recalibration_filename = base_name + ".indels.recal",
      tranches_filename = base_name + ".indels.tranches",
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

  if (num_samples_loaded > snps_variant_recalibration_threshold) {
    call Tasks.SNPsVariantRecalibratorCreateModel {
      input:
        sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcf,
        sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idx,
        recalibration_filename = base_name + ".snps.recal",
        tranches_filename = base_name + ".snps.tranches",
        model_report_filename = base_name + ".snps.model.report",
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
        sample_every_nth_variant = SNP_VQSR_sample_every_nth_variant,
        maximum_training_variants = SNP_VQSR_maximum_training_variants
    }

    scatter (idx in range(length(sites_only_variant_filtered_vcfs))) {
      call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
        input:
          sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcfs[idx],
          sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idxs[idx],
          recalibration_filename = base_name + ".snps." + idx + ".recal",
          tranches_filename = base_name + ".snps." + idx + ".tranches",
          model_report = SNPsVariantRecalibratorCreateModel.model_report,
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
          machine_mem_gb = SNP_VQSR_mem_gb_override
      }
    }

    call Tasks.GatherTranches as SNPGatherTranches {
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
    call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorClassic {
      input:
        sites_only_variant_filtered_vcf = sites_only_variant_filtered_vcf,
        sites_only_variant_filtered_vcf_index = sites_only_variant_filtered_vcf_idx,
        recalibration_filename = base_name + ".snps.recal",
        tranches_filename = base_name + ".snps.tranches",
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
  }

  output {
    File snps_variant_recalibration_file = select_first([MergeRecalibrationFiles.output_vcf, SNPsVariantRecalibratorClassic.recalibration])
    File snps_variant_recalibration_file_index = select_first([MergeRecalibrationFiles.output_vcf_index, SNPsVariantRecalibratorClassic.recalibration_index])
    File snps_variant_tranches_file = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches])
    File indels_variant_recalibration_file = IndelsVariantRecalibrator.recalibration
    File indels_variant_recalibration_file_index = IndelsVariantRecalibrator.recalibration_index
    File indels_variant_tranches_file = IndelsVariantRecalibrator.tranches
    Array[File] monitoring_logs = select_all(
                                  flatten(
                                   [
                                   [IndelsVariantRecalibrator.monitoring_log],
                                   [SNPsVariantRecalibratorCreateModel.monitoring_log],
                                   select_first([SNPsVariantRecalibratorScattered.monitoring_log, []]),
                                   [SNPGatherTranches.monitoring_log],
                                   [MergeRecalibrationFiles.monitoring_log],
                                   [SNPsVariantRecalibratorClassic.monitoring_log]
                                   ]))
  }

}

