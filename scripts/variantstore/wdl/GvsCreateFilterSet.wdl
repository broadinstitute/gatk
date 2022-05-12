version 1.0

import "GvsWarpTasks.wdl" as Tasks
import "GvsUtils.wdl" as Utils

workflow GvsCreateFilterSet {
  input {
    Boolean go = true
    String dataset_name
    String project_id

    String filter_set_name
    Array[String] indel_recalibration_annotation_values = ["AS_FS", "AS_ReadPosRankSum", "AS_MQRankSum", "AS_QD", "AS_SOR"]
    Int scatter_count
    Array[String] snp_recalibration_annotation_values = ["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_MQ", "AS_SOR"]

    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File gatk_override = "gs://broad-dsp-spec-ops/scratch/bigquery-jointcalling/jars/ah_var_store_20220415/gatk-package-4.2.0.0-492-g1387d47-SNAPSHOT-local.jar"

    Int? INDEL_VQSR_max_gaussians_override = 4
    Int? INDEL_VQSR_mem_gb_override
    String? service_account_json_path
    Int? SNP_VQSR_max_gaussians_override = 6
    Int? SNP_VQSR_mem_gb_override
  }

  # this is the minimum number of samples where the SNP model will be created and applied in separate tasks
  # (SNPsVariantRecalibratorClassic vs. SNPsVariantRecalibratorCreateModel and SNPsVariantRecalibratorScattered)
  Int snps_variant_recalibration_threshold = 20000


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
  String fq_tranches_destination_table = "~{project_id}.~{dataset_name}.filter_set_tranches"
  String fq_filter_sites_destination_table = "~{project_id}.~{dataset_name}.filter_set_sites"

  call Utils.GetBQTableLastModifiedDatetime as SamplesTableDatetimeCheck {
    input:
      query_project = project_id,
      fq_table = fq_sample_table,
      service_account_json_path = service_account_json_path
  }

  call GetNumSamplesLoaded {
    input:
      fq_sample_table = fq_sample_table,
      fq_sample_table_lastmodified_timestamp = SamplesTableDatetimeCheck.last_modified_timestamp,
      service_account_json_path = service_account_json_path,
      project_id = project_id
  }

  call Utils.SplitIntervals {
    input:
      intervals = interval_list,
      ref_fasta = reference,
      ref_fai = reference_index,
      ref_dict = reference_dict,
      scatter_count = scatter_count,
      gatk_override = gatk_override
  }

  scatter(i in range(length(SplitIntervals.interval_files))) {
    call ExtractFilterTask {
      input:
        gatk_override             = gatk_override,
        reference                 = reference,
        reference_index           = reference_index,
        reference_dict            = reference_dict,
        fq_sample_table           = fq_sample_table,
        intervals                 = SplitIntervals.interval_files[i],
        fq_alt_allele_table       = fq_alt_allele_table,
        excess_alleles_threshold  = 1000000,
        output_file               = "${filter_set_name}_${i}.vcf.gz",
        service_account_json_path = service_account_json_path,
        query_project             = project_id,
        dataset_id                = dataset_name,
    }
  }

  call Utils.MergeVCFs as MergeVCFs {
    input:
      input_vcfs = ExtractFilterTask.output_vcf,
      output_vcf_name = "${filter_set_name}.vcf.gz",
      preemptible_tries = 3,
  }

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
  }

  call PopulateFilterSetInfo {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      snp_recal_file = select_first([MergeRecalibrationFiles.output_vcf, SNPsVariantRecalibratorClassic.recalibration]),
      snp_recal_file_index = select_first([MergeRecalibrationFiles.output_vcf_index, SNPsVariantRecalibratorClassic.recalibration_index]),
      indel_recal_file = IndelsVariantRecalibrator.recalibration,
      indel_recal_file_index = IndelsVariantRecalibrator.recalibration_index,
      fq_info_destination_table = fq_info_destination_table,
      service_account_json_path = service_account_json_path,
      query_project = project_id
  }

  call PopulateFilterSetSites {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      sites_only_variant_filtered_vcf = MergeVCFs.output_vcf,
      sites_only_variant_filtered_vcf_index = MergeVCFs.output_vcf_index,
      fq_filter_sites_destination_table = fq_filter_sites_destination_table,
      service_account_json_path = service_account_json_path,
      query_project = project_id
  }

  call PopulateFilterSetTranches {
    input:
      gatk_override = gatk_override,
      filter_set_name = filter_set_name,
      snp_recal_tranches = select_first([SNPGatherTranches.tranches_file, SNPsVariantRecalibratorClassic.tranches]),
      indel_recal_tranches = IndelsVariantRecalibrator.tranches,
      fq_tranches_destination_table = fq_tranches_destination_table,
      service_account_json_path = service_account_json_path,
      query_project = project_id
  }

  output {
    File output_vcf = MergeVCFs.output_vcf
    File output_vcf_idx = MergeVCFs.output_vcf_index
    Boolean done = true
  }
}

################################################################################
task GetNumSamplesLoaded {
  input {
    String fq_sample_table
    String fq_sample_table_lastmodified_timestamp
    String? service_account_json_path
    String project_id
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e

    if [ ~{has_service_account_file} = 'true' ]; then
    gsutil cp ~{service_account_json_path} local.service_account.json
    gcloud auth activate-service-account --key-file=local.service_account.json
    gcloud config set project ~{project_id}
    fi

    echo "project_id = ~{project_id}" > ~/.bigqueryrc
    bq query --location=US --project_id=~{project_id} --format=csv --use_legacy_sql=false \
    'SELECT COUNT(*) as num_rows FROM `~{fq_sample_table}` WHERE is_loaded = true' > num_rows.csv

    NUMROWS=$(python3 -c "csvObj=open('num_rows.csv','r');csvContents=csvObj.read();print(csvContents.split('\n')[1]);")

    [[ $NUMROWS =~ ^[0-9]+$ ]] && echo $NUMROWS || exit 1
  >>>

  output {
    Int num_samples = read_int(stdout())
  }

  runtime {
    docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:305.0.0"
    memory: "3 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
    cpu: 1
  }
}

task ExtractFilterTask {
  input {
    File reference
    File reference_index
    File reference_dict

    String fq_sample_table

    File intervals

    String fq_alt_allele_table
    String output_file
    Int? excess_alleles_threshold

    # Runtime Options:
    File? gatk_override
    String? service_account_json_path
    String query_project
    String dataset_id
  }


  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project}
    fi

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
      --project-id ~{query_project}
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    memory: "7 GB"
    disks: "local-disk 10 HDD"
    bootDiskSizeGb: 15
    preemptible: 3
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

    String? service_account_json_path
    String query_project

    File? gatk_override
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -eo pipefail

    if [ ~{has_service_account_file} = 'true' ]; then
      gsutil cp ~{service_account_json_path} local.service_account.json
      export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
      gcloud auth activate-service-account --key-file=local.service_account.json
      gcloud config set project ~{query_project}
    fi

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
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_info = read_string("status_load_filter_set_info")
  }
}

task PopulateFilterSetSites {
  input {
    String filter_set_name
    String fq_filter_sites_destination_table

    File sites_only_variant_filtered_vcf
    File sites_only_variant_filtered_vcf_index

    String? service_account_json_path
    String query_project

    File? gatk_override
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -eo pipefail

    if [ ~{has_service_account_file} = 'true' ]; then
    gsutil cp ~{service_account_json_path} local.service_account.json
    export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
    gcloud auth activate-service-account --key-file=local.service_account.json
    gcloud config set project ~{query_project}
    fi

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
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
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

task PopulateFilterSetTranches {
  input {
    File? gatk_override

    String filter_set_name
    String fq_tranches_destination_table

    File snp_recal_tranches
    File indel_recal_tranches

    String? service_account_json_path
    String query_project
  }

  String has_service_account_file = if (defined(service_account_json_path)) then 'true' else 'false'

  command <<<
    set -eo pipefail

    if [ ~{has_service_account_file} = 'true' ]; then
    gsutil cp ~{service_account_json_path} local.service_account.json
    export GOOGLE_APPLICATION_CREDENTIALS=local.service_account.json
    gcloud auth activate-service-account --key-file=local.service_account.json
    gcloud config set project ~{query_project}
    fi

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    cat ~{snp_recal_tranches} ~{indel_recal_tranches} | grep -v targetTruthSensitivity | grep -v "#" | awk -v CALLSET=~{filter_set_name} '{ print CALLSET "," $0 }' > ~{filter_set_name}.tranches_load.csv

    # BQ load likes a : instead of a . after the project
    bq_table=$(echo ~{fq_tranches_destination_table} | sed s/\\./:/)

    echo "Loading combined tranches CSV into ~{fq_tranches_destination_table}"
    bq load --project_id=~{query_project} --skip_leading_rows 0 -F "," \
    --schema "filter_set_name:string,target_truth_sensitivity:float,num_known:integer,num_novel:integer,known_ti_tv:float,novel_ti_tv:float,min_vqslod:float,filter_name:string,model:string,accessible_truth_sites:integer,calls_at_truth_sites:integer,truth_sensitivity:float" \
    ${bq_table} \
    ~{filter_set_name}.tranches_load.csv > status_load_filter_set_tranches
  >>>

  runtime {
    docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    memory: "3500 MB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: 15
    preemptible: 0
    cpu: 1
  }

  output {
    String status_load_filter_set_tranches = read_string("status_load_filter_set_tranches")
  }
}
