version 1.0

import "GvsUtils.wdl" as Utils
import "GvsExtractCohortFromSampleNames.wdl" as GvsExtractSubCohortVCFs


workflow GvsCalculatePrecisionAndSensitivity {
  input {
    String call_set_identifier
    String dataset_name
    String filter_set_name
    File interval_list = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list"
    File? vcf_eval_bed_file
    Array[String] chromosomes = ["chr20"]
    String project_id
    Array[String] sample_names

    Array[File] truth_vcfs
    Array[File] truth_vcf_indices
    Array[File] truth_beds

    Int? extract_scatter_count_override
    File ref_fasta = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

    String? basic_docker
    String? cloud_sdk_docker
    String? gatk_docker
    String? git_branch_or_tag
    String? gatk_override
    String? gotc_imputation_docker
    String? real_time_genomics_docker
    String? variants_docker
  }

  parameter_meta {
    call_set_identifier: "The name of the callset for which we are calculating precision and sensitivity."
    chromosomes: "The chromosome(s) on which to analyze precision and sensitivity. The default value for this is `['chr20']`."
    dataset_name: "The GVS BigQuery dataset name."
    filter_set_name: "The filter_set_name used to generate the callset."
    interval_list: "The intervals over which to extract VCFs for calculating precision and sensitivity."
    project_id: "The Google Project ID where the GVS lives."
    sample_names: "A list of the sample names that are controls and that will be used for the analysis. For every element on the list of sample names there must be a corresponding element on the list of `truth_vcfs`, `truth_vcf_indices`, and `truth_beds`."
    truth_vcfs: "A list of the VCFs that contain the truth data used for analyzing the samples in `sample_names`."
    truth_vcf_indices: "A list of the VCF indices for the truth data VCFs supplied above."
    truth_beds: "A list of the bed files for the truth data used for analyzing the samples in `sample_names`."
    ref_fasta: "The cloud path for the reference fasta sequence."
    vcf_eval_bed_file: "Optional bed file for EvaluateVcf; if passed, will be used instead of chromosomes."
  }

  String output_basename = call_set_identifier + "_PS"

  # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
  # no calling WDLs that might supply `git_hash`).
  call Utils.GetToolVersions {
    input:
      git_branch_or_tag = git_branch_or_tag,
  }

  String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
  String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_real_time_genomics_docker = select_first([real_time_genomics_docker, GetToolVersions.real_time_genomics_docker])
  String effective_gotc_imputation_docker = select_first([gotc_imputation_docker, GetToolVersions.gotc_imputation_docker])

  if ((length(sample_names) != length(truth_vcfs)) || (length(sample_names) != length(truth_vcf_indices)) || (length(sample_names) != length(truth_beds))) {
    call Utils.TerminateWorkflow {
      input:
        message = "The inputs 'sample_names', 'truth_vcfs', 'truth_vcf_indices', and 'truth_beds' must all contain the same number of elements",
        basic_docker = effective_basic_docker,
    }
  }

  call GvsExtractSubCohortVCFs.GvsExtractCohortFromSampleNames as GenerateControlVCFs {
    input:
      cohort_sample_names_array = sample_names,
      gvs_project = project_id,
      gvs_dataset = dataset_name,
      call_set_identifier = call_set_identifier,
      filter_set_name = filter_set_name,
      control_samples = true,
      interval_list = interval_list,
      extract_scatter_count_override = extract_scatter_count_override,
      cloud_sdk_docker = effective_cloud_sdk_docker,
      gatk_docker = effective_gatk_docker,
      gatk_override = gatk_override,
      git_branch_or_tag = GetToolVersions.git_hash,
      variants_docker = effective_variants_docker,
  }

  call GatherVcfs {
    input:
      input_vcfs = GenerateControlVCFs.output_vcfs,
      output_basename = output_basename,
      total_vcfs_size_mb = GenerateControlVCFs.total_vcfs_size_mb,
      gatk_docker = effective_gatk_docker,
  }

  scatter(i in range(length(sample_names))) {
    String sample_name = sample_names[i]
    String output_sample_basename = output_basename + "." + sample_name

    call SelectVariants {
      input:
        input_vcf = GatherVcfs.output_vcf,
        input_vcf_index = GatherVcfs.output_vcf_index,
        sample_name = sample_name,
        output_basename = output_sample_basename,
        gatk_docker = effective_gatk_docker,
    }

    call Add_AS_MAX_VQS_SCORE_ToVcf {
      input:
        input_vcf = SelectVariants.output_vcf,
        output_basename = output_sample_basename + ".maxas",
        variants_docker = effective_variants_docker,
    }

    call IsVQSRLite {
      input:
        input_vcf = Add_AS_MAX_VQS_SCORE_ToVcf.output_vcf,
        basic_docker = effective_basic_docker,
    }

    call BgzipAndTabix {
      input:
        input_vcf = Add_AS_MAX_VQS_SCORE_ToVcf.output_vcf,
        output_basename = output_sample_basename + ".maxas",
        gotc_imputation_docker = effective_gotc_imputation_docker,
    }

    call EvaluateVcf as EvaluateVcfFiltered {
      input:
        input_vcf = BgzipAndTabix.output_vcf,
        input_vcf_index = BgzipAndTabix.output_vcf_index,
        truth_vcf = truth_vcfs[i],
        truth_vcf_index = truth_vcf_indices[i],
        truth_bed = truth_beds[i],
        vcf_eval_bed_file = vcf_eval_bed_file,
        chromosomes = chromosomes,
        output_basename = sample_name + "-bq_roc_filtered",
        is_vqsr_lite = IsVQSRLite.is_vqsr_lite,
        ref_fasta = ref_fasta,
        real_time_genomics_docker = effective_real_time_genomics_docker,
    }

    call EvaluateVcf as EvaluateVcfAll {
      input:
        input_vcf = BgzipAndTabix.output_vcf,
        input_vcf_index = BgzipAndTabix.output_vcf_index,
        truth_vcf = truth_vcfs[i],
        truth_vcf_index = truth_vcf_indices[i],
        truth_bed = truth_beds[i],
        vcf_eval_bed_file = vcf_eval_bed_file,
        chromosomes = chromosomes,
        all_records = true,
        output_basename = sample_name + "-bq_all",
        is_vqsr_lite = IsVQSRLite.is_vqsr_lite,
        ref_fasta = ref_fasta,
        real_time_genomics_docker = effective_real_time_genomics_docker,
    }
  }

  call CollateReports {
    input:
      all_snp_reports = EvaluateVcfAll.snp_report,
      all_indel_reports = EvaluateVcfAll.indel_report,
      filtered_snp_reports = EvaluateVcfFiltered.snp_report,
      filtered_indel_reports = EvaluateVcfFiltered.indel_report,
      basic_docker = effective_basic_docker,
  }

  output {
    File report = CollateReports.report
    Array[Array[File]] filtered_eval_outputs = EvaluateVcfFiltered.outputs
    Array[Array[File]] all_eval_outputs = EvaluateVcfAll.outputs
    String recorded_git_hash = GetToolVersions.git_hash
  }
}

task GatherVcfs {
  input {
    Array[File] input_vcfs
    String output_basename
    Float total_vcfs_size_mb

    String gatk_docker
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil((3*total_vcfs_size_mb)/1024) + 500
  }

  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    CHR_VCFS_ARG=""
    for file in ~{sep=' ' input_vcfs}
    do
      if [ -s $file ]; then
        CHR_VCFS_ARG+=" --INPUT $file "
      fi
    done
    echo $CHR_VCFS_ARG

    # --REORDER_INPUT_BY_FIRST_VARIANT means that the vcfs supplied here need not be ordered by location.
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      GatherVcfs \
        --REORDER_INPUT_BY_FIRST_VARIANT \
        $CHR_VCFS_ARG \
        --OUTPUT ~{output_basename}.vcf.gz

    tabix ~{output_basename}.vcf.gz
  >>>

  runtime {
    docker: gatk_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 1
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task SelectVariants {
  input {
    File input_vcf
    File input_vcf_index
    String sample_name

    String output_basename

    String gatk_docker
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(2*size(input_vcf, "GiB")) + 500
  }
  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      SelectVariants \
        -V ~{input_vcf} \
        --sample-name ~{sample_name} \
        --select-type-to-exclude NO_VARIATION \
        -O ~{output_basename}.vcf.gz
  >>>

  runtime {
    docker: gatk_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 1
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task Add_AS_MAX_VQS_SCORE_ToVcf {
  input {
    File input_vcf
    String output_basename

    String variants_docker
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(2*size(input_vcf, "GiB")) + 500
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    python3 /app/add_max_as_vqs_score.py ~{input_vcf} > ~{output_basename}.vcf
  >>>
  runtime {
    docker: variants_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File output_vcf = "~{output_basename}.vcf"
  }
}

task IsVQSRLite {
  input {
    File input_vcf
    String basic_docker
  }

  String is_vqsr_lite_file = "is_vqsr_lite_file.txt"

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # See if there are any non-header lines that contain the string 'CALIBRATION_SENSITIVITY'. If so, grep will return 0 else 1
    set +o errexit
    grep -v '^#' ~{input_vcf} | grep CALIBRATION_SENSITIVITY > /dev/null
    if [[ $? -eq 0 ]]; then
      echo "true" > ~{is_vqsr_lite_file}
    else
      echo "false" > ~{is_vqsr_lite_file}
    fi
    set -o errexit
  >>>

  runtime {
    docker: basic_docker
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Boolean is_vqsr_lite = read_boolean(is_vqsr_lite_file)
  }
}

task BgzipAndTabix {
  input {
    File input_vcf
    String output_basename

    String gotc_imputation_docker
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(3 * size(input_vcf, "GiB")) + 500
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # note that bgzip has an option (-i) to index the bgzipped output, but this file is not a tabix file
    # note also that we use '-c' so that bgzip doesn't create the bgzipped file in place, rather it's in a location
    # where it's easy to output from the task.
    bgzip -c ~{input_vcf} > ~{output_basename}.vcf.gz
    tabix ~{output_basename}.vcf.gz
  >>>
  runtime {
    docker: gotc_imputation_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"

  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task EvaluateVcf {
  input {
    File input_vcf
    File input_vcf_index
    File truth_vcf
    File truth_vcf_index
    File truth_bed
    File? vcf_eval_bed_file
    Array[String] chromosomes

    Boolean all_records = false

    File ref_fasta

    String output_basename

    Boolean is_vqsr_lite

    String real_time_genomics_docker
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(2 * size(ref_fasta, "GiB")) + 500
  }

  String max_score_field_tag = if (is_vqsr_lite == true) then 'MAX_CALIBRATION_SENSITIVITY' else 'MAX_AS_VQSLOD'

  command <<<
    chromosomes=( ~{sep=' ' chromosomes} )

    echo "Creating .bed file to control which chromosomes should be evaluated."
    for i in "${chromosomes[@]}"
    do
    echo "$i	0	300000000" >> chromosomes.to.eval.txt
    done

    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    rtg format --output human_REF_SDF ~{ref_fasta}

    rtg vcfeval \
      --bed-regions ~{if defined(vcf_eval_bed_file) then vcf_eval_bed_file else "chromosomes.to.eval.txt"} \
      ~{if all_records then "--all-records" else ""} \
      --roc-subset snp,indel \
      --vcf-score-field=INFO.~{max_score_field_tag} \
      ~{if is_vqsr_lite then "--sort-order ascending" else "--sort-order descending"} \
      -t human_REF_SDF \
      -b ~{truth_vcf} \
      -e ~{truth_bed}\
      -c ~{input_vcf} \
      -o ~{output_basename}

    # Touch a file with the name of the sample in that directory, so that it's identifiable among the globbed outputs.
    touch ~{output_basename}/~{output_basename}

    touch snp_report.txt
    touch indel_report.txt
    for type in "snp" "indel"
      do
        d=$(cat ~{output_basename}/${type}_roc.tsv.gz | gunzip | tail -1 | cut -f2,3,5,6,7)
        echo -e "~{output_basename}\t$type\t$d" >> ${type}_report.txt
    done
  >>>
  runtime {
    docker: real_time_genomics_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"

  }
  output {
    File snp_report = "snp_report.txt"
    File indel_report = "indel_report.txt"
    Array[File] outputs = glob("~{output_basename}/*")
  }
}

task CollateReports {
  input {
    Array[File] all_snp_reports
    Array[File] all_indel_reports
    Array[File] filtered_snp_reports
    Array[File] filtered_indel_reports
    String basic_docker
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    echo "Sample  Type  TP  FP  FN  Precision Sensitivity"
    while read -r a;
    do
      cat $a
    done < ~{write_lines(all_snp_reports)}

    while read -r a;
    do
      cat $a
    done < ~{write_lines(all_indel_reports)}

    while read -r a;
    do
      cat $a
    done < ~{write_lines(filtered_snp_reports)}

    while read -r a;
    do
      cat $a
    done < ~{write_lines(filtered_indel_reports)}
  >>>

  runtime {
    docker: basic_docker
    disks: "local-disk 100 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    File report = stdout()
  }
}


task CountInputVcfs {
  input {
    File input_vcf_fofn
    String basic_docker
  }
  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    wc -l < ~{input_vcf_fofn} > num_vcfs.txt
  >>>
  output {
    Int num_vcfs = read_int("num_vcfs.txt")
  }
  runtime {
    docker: basic_docker
  }
}
