version 1.0

import "GvsUtils.wdl" as Utils
import "GvsVQSRClassic.wdl" as VQSRClassic
import "../../vcf_site_level_filtering_wdl/JointVcfFiltering.wdl" as VQSRLite

workflow GvsCreateFilterSet {
  input {
    Boolean go = true
    String dataset_name
    String project_id
    String call_set_identifier

    String filter_set_name

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File? gatk_override

    Boolean use_classic_VQSR = true

    Int? INDEL_VQSR_CLASSIC_max_gaussians_override = 4
    Int? INDEL_VQSR_CLASSIC_mem_gb_override
    Int? SNP_VQSR_CLASSIC_max_gaussians_override = 6
    Int? SNP_VQSR_CLASSIC_mem_gb_override

    RuntimeAttributes? vqsr_lite_extract_runtime_attributes = {"command_mem_gb": 27}
    RuntimeAttributes? vqsr_lite_train_runtime_attributes = {"command_mem_gb": 27}
    RuntimeAttributes? vqsr_lite_score_runtime_attributes = {"command_mem_gb": 15}
  }

  File reference = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
  File reference_dict = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
  File reference_index = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

  # fully-qualified table names
  String fq_sample_table = "~{project_id}.~{dataset_name}.sample_info"
  String fq_alt_allele_table = "~{project_id}.~{dataset_name}.alt_allele"
  String fq_filter_sites_destination_table =    "~{project_id}.~{dataset_name}.filter_set_sites"
  String fq_filter_set_info_destination_table = "~{project_id}.~{dataset_name}.filter_set_info"
  String fq_tranches_destination_table =        "~{project_id}.~{dataset_name}.filter_set_tranches"

  String filter_set_info_destination_table_schema = "filter_set_name:string,type:string,location:integer,ref:string,alt:string,calibration_sensitivity:float,vqslod:float,culprit:string,training_label:string,yng_status:string"

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
        input_vcfs = ExtractFilterTask.output_vcf,
        input_vcf_idxs = ExtractFilterTask.output_vcf_index,
        sites_only_vcf = MergeVCFs.output_vcf,
        sites_only_vcf_idx = MergeVCFs.output_vcf_index,
        output_prefix = filter_set_name,
        annotations = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"],
        resource_args = "--resource:hapmap,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz --resource:omni,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz --resource:1000G,training=true,calibration=false gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:mills,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource:axiom,training=true,calibration=false gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        extract_extra_args = "-L ${interval_list} --use-allele-specific-annotations",
        score_extra_args = "-L ${interval_list} --use-allele-specific-annotations",
        extract_runtime_attributes = vqsr_lite_extract_runtime_attributes,
        train_runtime_attributes = vqsr_lite_train_runtime_attributes,
        score_runtime_attributes = vqsr_lite_score_runtime_attributes,
        gatk_docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0",
        gatk_override = gatk_override,
        monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"
    }

    call Utils.MergeVCFs as MergeScoredVCFs {
      input:
        input_vcfs = JointVcfFiltering.scored_vcfs,
        gather_type = "CONVENTIONAL",
        output_vcf_name = "${filter_set_name}.vrecalibration.gz",
        preemptible_tries = 3,
    }

    # These calls to SelectVariants are being added for two reasons
    # 1) The snps_variant_scored_vcf and indels_variant_scored_vcf output by JointVcfFiltering contains ALL variants,
    #     but are currently ONLY annotating SNPs and INDELs respectively.
    # 2) Those output VCFs also contain filtered sites (sites at which the FILTER field set to anything other than '.' or 'PASS')
    #     which we don't want to put into the filter_set_info table.
    call Utils.SelectVariants as CreateFilteredScoredSNPsVCF {
      input:
        input_vcf = MergeScoredVCFs.output_vcf,
        input_vcf_index = MergeScoredVCFs.output_vcf_index,
        type_to_include = "SNP",
        exclude_filtered = true,
        output_basename = "${filter_set_name}.filtered.scored.snps"
    }

    call Utils.SelectVariants as CreateFilteredScoredINDELsVCF {
      input:
        input_vcf = MergeScoredVCFs.output_vcf,
        input_vcf_index = MergeScoredVCFs.output_vcf_index,
        type_to_include = "INDEL",
        exclude_filtered = true,
        output_basename = "${filter_set_name}.filtered.scored.indels"
    }

    call PopulateFilterSetInfo {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        fq_filter_set_info_destination_table = fq_filter_set_info_destination_table,
        filter_schema = filter_set_info_destination_table_schema,
        snp_recal_file = CreateFilteredScoredSNPsVCF.output_vcf,
        snp_recal_file_index = CreateFilteredScoredSNPsVCF.output_vcf_index,
        indel_recal_file = CreateFilteredScoredINDELsVCF.output_vcf,
        indel_recal_file_index = CreateFilteredScoredINDELsVCF.output_vcf_index,
        project_id = project_id,
        useClassic = false
    }
  }

  if (use_classic_VQSR) {
    call VQSRClassic.JointVcfFiltering as VQSRClassic {
      input:
        base_name = filter_set_name,
        num_samples_loaded = GetNumSamplesLoaded.num_samples,
        sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
        sites_only_variant_filtered_vcf_idx = MergeVCFs.output_vcf_index,
        sites_only_variant_filtered_vcfs = ExtractFilterTask.output_vcf,
        sites_only_variant_filtered_vcf_idxs = ExtractFilterTask.output_vcf_index,
        INDEL_VQSR_max_gaussians_override = INDEL_VQSR_CLASSIC_max_gaussians_override,
        INDEL_VQSR_mem_gb_override = INDEL_VQSR_CLASSIC_mem_gb_override,
        SNP_VQSR_max_gaussians_override = SNP_VQSR_CLASSIC_max_gaussians_override,
        SNP_VQSR_mem_gb_override = SNP_VQSR_CLASSIC_mem_gb_override
    }

    call PopulateFilterSetInfo as PopulateFilterSetInfoClassic {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        fq_filter_set_info_destination_table = fq_filter_set_info_destination_table,
        filter_schema = filter_set_info_destination_table_schema,
        snp_recal_file = VQSRClassic.snps_variant_recalibration_file,
        snp_recal_file_index = VQSRClassic.snps_variant_recalibration_file_index,
        indel_recal_file = VQSRClassic.indels_variant_recalibration_file,
        indel_recal_file_index = VQSRClassic.indels_variant_recalibration_file_index,
        project_id = project_id,
        useClassic = true
    }

    call PopulateFilterSetTranches {
      input:
        gatk_override = gatk_override,
        filter_set_name = filter_set_name,
        snp_recal_tranches = VQSRClassic.snps_variant_tranches_file,
        indel_recal_tranches = VQSRClassic.indels_variant_tranches_file,
        fq_tranches_destination_table = fq_tranches_destination_table,
        project_id = project_id
    }
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

  call Utils.SummarizeTaskMonitorLogs as SummarizeItAll {
    input:
      inputs = select_all(
               flatten(
               [
               [SamplesTableDatetimeCheck.monitoring_log],
               [GetNumSamplesLoaded.monitoring_log],
               [SplitIntervals.monitoring_log],
               [AltAlleleTableDatetimeCheck.monitoring_log],
               ExtractFilterTask.monitoring_log,
               [MergeVCFs.monitoring_log],
               select_first([JointVcfFiltering.monitoring_logs, []]),
               [MergeScoredVCFs.monitoring_log],
               [CreateFilteredScoredSNPsVCF.monitoring_log],
               [CreateFilteredScoredINDELsVCF.monitoring_log],
               [PopulateFilterSetInfo.monitoring_log],
               select_first([VQSRClassic.monitoring_logs, []]),
               [PopulateFilterSetInfoClassic.monitoring_log],
               [PopulateFilterSetTranches.monitoring_log],
               [PopulateFilterSetSites.monitoring_log]
               ]
               )
               )
  }

  output {
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_idx = MergeVCFs.output_vcf_index
    File monitoring_summary = SummarizeItAll.monitoring_summary
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
  File monitoring_script = "gs://gvs_quickstart_storage/cromwell_monitoring_script.sh"

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    bash ~{monitoring_script} > monitoring.log &

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
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_06_07_7a8ddfc"
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
    File monitoring_log = "monitoring.log"
  }
}

task PopulateFilterSetInfo {
  input {
    String filter_set_name
    String filter_schema
    String fq_filter_set_info_destination_table
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
    bq_table=$(echo ~{fq_filter_set_info_destination_table} | sed s/\\./:/)

    echo "Loading combined TSV into ~{fq_filter_set_info_destination_table}"
    bq --apilog=false load --project_id=~{project_id} --skip_leading_rows 0 -F "tab" \
      --range_partitioning=location,0,26000000000000,6500000000 \
      --clustering_fields=location \
      --schema "~{filter_schema}" \
      ${bq_table} \
      ~{filter_set_name}.filter_set_load.tsv > status_load_filter_set_info
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_06_07_7a8ddfc"
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
    bq --apilog=false load --project_id=~{project_id} --skip_leading_rows 1 -F "tab" \
    --range_partitioning=location,0,26000000000000,6500000000 \
    --clustering_fields=location \
    --schema "filter_set_name:string,location:integer,filters:string" \
    ${bq_table} \
    ~{filter_set_name}.filter_sites_load.tsv > status_load_filter_set_sites
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_06_07_7a8ddfc"
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
    bq --apilog=false load --project_id=~{project_id} --skip_leading_rows 0 -F "," \
    --schema "filter_set_name:string,target_truth_sensitivity:float,num_known:integer,num_novel:integer,known_ti_tv:float,novel_ti_tv:float,min_vqslod:float,filter_name:string,model:string,accessible_truth_sites:integer,calls_at_truth_sites:integer,truth_sensitivity:float" \
    ${bq_table} \
    ~{filter_set_name}.tranches_load.csv > status_load_filter_set_tranches
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_2023_06_07_7a8ddfc"
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
