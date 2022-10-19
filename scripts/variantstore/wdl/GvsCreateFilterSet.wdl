version 1.0

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
  String fq_info_destination_table = "~{project_id}.~{dataset_name}.vqsr_lite_filter_set_info"
  String fq_filter_sites_destination_table = "~{project_id}.~{dataset_name}.vqsr_lite_filter_set_sites"

  call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
    input:
      query_project = project_id,
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
      query_project = project_id,
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
        query_project              = project_id,
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

  call VQSRLite.JointVcfFiltering as JointVcfFiltering {
    input:
      vcf = ExtractFilterTask.output_vcf,
      vcf_index = ExtractFilterTask.output_vcf_index,
      sites_only_vcf = MergeVCFs.output_vcf,
      sites_only_vcf_index = MergeVCFs.output_vcf_index,
      basename = filter_set_name,
      gatk_docker = "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:UG_feature_branch_v4",
      extract_interval_list = interval_list,
      score_interval_list = interval_list,
      snp_annotations = "-A AS_QD -A AS_MQRankSum -A AS_ReadPosRankSum -A AS_FS -A AS_MQ -A AS_SOR",
      indel_annotations = "-A AS_FS -A AS_ReadPosRankSum -A AS_MQRankSum -A AS_QD -A AS_SOR",
      use_allele_specific_annotations = true,
      gatk_override = "gs://gvs-internal-scratch/rsa/gatk-package-4.2.0.0-614-g971d82f-SNAPSHOT-local.jar"
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

  call PopulateFilterSetInfo {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      snp_recal_file = MergeSNPScoredVCFs.output_vcf,
      snp_recal_file_index = MergeSNPScoredVCFs.output_vcf_index,
      indel_recal_file = MergeINDELScoredVCFs.output_vcf,
      indel_recal_file_index = MergeINDELScoredVCFs.output_vcf_index,
      fq_info_destination_table = fq_info_destination_table,
      query_project = project_id
  }

  call PopulateFilterSetSites {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
      sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
      fq_filter_sites_destination_table = fq_filter_sites_destination_table,
      query_project = project_id
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
    String query_project
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
      --project-id ~{query_project} \
      --cost-observability-tablename ~{cost_observability_tablename} \
      --call-set-identifier ~{call_set_identifier} \
      --wdl-step GvsCreateFilterSet \
      --wdl-call ExtractFilterTask \
      --shard-identifier ~{intervals_name}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2022_10_17_2a8c210ac35094997603259fa1cd784486b92e42"
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
    String fq_info_destination_table

    File snp_recal_file
    File snp_recal_file_index
    File indel_recal_file
    File indel_recal_file_index

    String query_project

    File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    File? gatk_override
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -eo pipefail

    bash ~{monitoring_script} > monitoring.log &

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    echo "Creating SNPs reacalibration file"
    gatk --java-options "-Xmx1g" \
      CreateFilteringFiles \
      --ref-version 38 \
      --filter-set-name ~{filter_set_name} \
      -mode SNP \
      -V ~{snp_recal_file} \
      -O ~{filter_set_name}.snps.recal.tsv

    echo "Creating INDELs reacalibration file"
    gatk --java-options "-Xmx1g" \
      CreateFilteringFiles \
      --ref-version 38 \
      --filter-set-name ~{filter_set_name} \
      -mode INDEL \
      -V ~{indel_recal_file} \
      -O ~{filter_set_name}.indels.recal.tsv

    # merge into a single file
    echo "Merging SNP + INDELs"
    cat ~{filter_set_name}.snps.recal.tsv ~{filter_set_name}.indels.recal.tsv | grep -v filter_set_name | grep -v "#"  > ~{filter_set_name}.filter_set_load.tsv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_info_destination_table} | sed s/\\./:/)

    echo "Loading combined TSV into ~{fq_info_destination_table}"
    bq load --project_id=~{query_project} --skip_leading_rows 0 -F "tab" \
      --range_partitioning=location,0,26000000000000,6500000000 \
      --clustering_fields=location \
      --schema "filter_set_name:string,type:string,location:integer,ref:string,alt:string,vqslod:float,culprit:string,training_label:string,yng_status:string" \
      ${bq_table} \
      ~{filter_set_name}.filter_set_load.tsv > status_load_filter_set_info
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2022_10_17_2a8c210ac35094997603259fa1cd784486b92e42"
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

    String query_project

    File? gatk_override
  }
  meta {
    # Not `volatile: true` since there shouldn't be a need to re-run this if there has already been a successful execution.
  }

  command <<<
    set -eo pipefail

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
    bq load --project_id=~{query_project} --skip_leading_rows 1 -F "tab" \
    --range_partitioning=location,0,26000000000000,6500000000 \
    --clustering_fields=location \
    --schema "filter_set_name:string,location:integer,filters:string" \
    ${bq_table} \
    ~{filter_set_name}.filter_sites_load.tsv > status_load_filter_set_sites
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2022_10_17_2a8c210ac35094997603259fa1cd784486b92e42"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_sites = read_string("status_load_filter_set_sites")
  }
}
