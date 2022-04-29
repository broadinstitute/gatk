version 1.0


task CheckSamplesUnique {
  input {
    File sample_name_map
    Int sample_num_threshold = 50
  }

  command {
    set -euo pipefail
    if [[ $(cut -f 1 ~{sample_name_map} | wc -l) -ne $(cut -f 1 ~{sample_name_map} | sort | uniq | wc -l) ]]
    then
      echo "Samples in the sample_name_map are not unique" 1>&2
      exit 1
    elif [[ $(cut -f 1 ~{sample_name_map} | wc -l) -lt ~{sample_num_threshold} ]]
    then
      echo "There are fewer than ~{sample_num_threshold} samples in the sample_name_map" 1>&2
      echo "Having fewer than ~{sample_num_threshold} samples means there likely isn't enough data to complete joint calling" 1>&2
      exit 1
    else
      echo true
    fi
  }

  output {
    Boolean samples_unique = read_boolean(stdout())
  }

  runtime {
    memory: "1 GiB"
    preemptible: 1
    disks: "local-disk 10 HDD"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
  }
}

task SplitIntervalList {

  input {
    File interval_list
    Int scatter_count
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    Boolean sample_names_unique_done
    Int disk_size
    String scatter_mode = "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW"
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    interval_list: {
      localization_optional: true
    }
  }

  command <<<
    gatk --java-options -Xms3g SplitIntervals \
      -L ~{interval_list} -O  scatterDir -scatter ~{scatter_count} -R ~{ref_fasta} \
      -mode ~{scatter_mode} --interval-merging-rule OVERLAPPING_ONLY
    >>>

  runtime {
    memory: "3.75 GiB"
    preemptible: 1
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
  }

  output {
    Array[File] output_intervals = glob("scatterDir/*")
  }
}

task ImportGVCFs {

  input {
    File sample_name_map
    File interval
    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String workspace_dir_name

    Int disk_size
    Int batch_size

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    set -euo pipefail

    rm -rf ~{workspace_dir_name}

    # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
    # using the --merge-input-intervals arg
    # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

    # The memory setting here is very important and must be several GiB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    gatk --java-options -Xms8g \
      GenomicsDBImport \
      --genomicsdb-workspace-path ~{workspace_dir_name} \
      --batch-size ~{batch_size} \
      -L ~{interval} \
      --sample-name-map ~{sample_name_map} \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf ~{workspace_dir_name}.tar ~{workspace_dir_name}
  >>>

  runtime {
    memory: "26 GiB"
    cpu: 4
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    docker: gatk_docker
    preemptible: 1
  }

  output {
    File output_genomicsdb = "~{workspace_dir_name}.tar"
  }
}

task GenotypeGVCFs {

  input {
    File workspace_tar
    File interval

    String output_vcf_filename

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    String dbsnp_vcf

    Int disk_size
    # This is needed for gVCFs generated with GATK3 HaplotypeCaller
    Boolean allow_old_rms_mapping_quality_annotation_data = false
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    interval: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    tar -xf ~{workspace_tar}
    WORKSPACE=$(basename ~{workspace_tar} .tar)

    gatk --java-options -Xms8g \
      GenotypeGVCFs \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      -D ~{dbsnp_vcf} \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --only-output-calls-starting-in-intervals \
      -V gendb://$WORKSPACE \
      -L ~{interval} \
      ~{true='--allow-old-rms-mapping-quality-annotation-data' false='' allow_old_rms_mapping_quality_annotation_data} \
      --merge-input-intervals
  >>>

  runtime {
    memory: "26 GiB"
    cpu: 2
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
  }
}

task GnarlyGenotyper {

  input {
    File workspace_tar
    File interval
    String output_vcf_filename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    String dbsnp_vcf

    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    interval: {
      localization_optional: true
    }
  }

  Int disk_size = ceil(size(workspace_tar, "GiB") + size(ref_fasta, "GiB") + size(dbsnp_vcf, "GiB") * 3)

  command <<<
    set -e

    tar -xf ~{workspace_tar}
    WORKSPACE=$( basename ~{workspace_tar} .tar)

    gatk --java-options -Xms8g \
      GnarlyGenotyper \
      -R ~{ref_fasta} \
      -O ~{output_vcf_filename} \
      --output-database-name annotationDB.vcf.gz \
      -D ~{dbsnp_vcf} \
      --only-output-calls-starting-in-intervals \
      -V gendb://$WORKSPACE \
      -L ~{interval} \
      -stand-call-conf 10 \
      --merge-input-intervals
  >>>

  runtime {
    memory: "26 GiB"
    cpu: 2
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_filename}"
    File output_vcf_index = "~{output_vcf_filename}.tbi"
    File output_database = "annotationDB.vcf.gz"
    File output_database_index = "annotationDB.vcf.gz.tbi"
  }
}

task HardFilterAndMakeSitesOnlyVcf {

  input {
    File vcf
    File vcf_index
    Float excess_het_threshold

    String variant_filtered_vcf_filename
    String sites_only_vcf_filename

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms3g \
      VariantFiltration \
      --filter-expression "ExcessHet > ~{excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ~{variant_filtered_vcf_filename} \
      -V ~{vcf}

    gatk --java-options -Xms3g \
      MakeSitesOnlyVcf \
      -I ~{variant_filtered_vcf_filename} \
      -O ~{sites_only_vcf_filename}
  >>>

  runtime {
    memory: "3.75 GiB"
    cpu: "1"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File variant_filtered_vcf = "~{variant_filtered_vcf_filename}"
    File variant_filtered_vcf_index = "~{variant_filtered_vcf_filename}.tbi"
    File sites_only_vcf = "~{sites_only_vcf_filename}"
    File sites_only_vcf_index = "~{sites_only_vcf_filename}.tbi"
  }
}

task IndelsVariantRecalibrator {

  input {
    String recalibration_filename
    String tranches_filename

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

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms24g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode INDEL \
      --max-gaussians ~{max_gaussians} \
      -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
      -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "26 GiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}

task SNPsVariantRecalibratorCreateModel {

  input {
    String recalibration_filename
    String tranches_filename
    Int downsampleFactor
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

    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms100g \
      VariantRecalibrator \
      -V ~{sites_only_variant_filtered_vcf} \
      -O ~{recalibration_filename} \
      --tranches-file ~{tranches_filename} \
      --trust-all-polymorphic \
      -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
      -an ~{sep=' -an ' recalibration_annotation_values} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      -mode SNP \
      --sample-every-Nth-variant ~{downsampleFactor} \
      --output-model ~{model_report_filename} \
      --max-gaussians ~{max_gaussians} \
      -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
      -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
      -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
      -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
  >>>

  runtime {
    memory: "104 GiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File model_report = "~{model_report_filename}"
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
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
    Int? machine_mem_gb

  }

  Int auto_mem = ceil(2 * size([sites_only_variant_filtered_vcf,
                              hapmap_resource_vcf,
                              omni_resource_vcf,
                              one_thousand_genomes_resource_vcf,
                              dbsnp_resource_vcf],
                      "GiB"))
  Int machine_mem = select_first([machine_mem_gb, if auto_mem < 7 then 7 else auto_mem])
  Int java_mem = machine_mem - 1


  String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

  command <<<
    set -euo pipefail

    MODEL_REPORT=~{model_report}

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
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File recalibration = "~{recalibration_filename}"
    File recalibration_index = "~{recalibration_filename}.idx"
    File tranches = "~{tranches_filename}"
  }
}

task GatherTranches {

  input {
    Array[File] tranches
    String output_filename
    String mode
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    tranches: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

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
      --output ~{output_filename}
  >>>

  runtime {
    memory: "7.5 GiB"
    cpu: "2"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File tranches = "~{output_filename}"
  }
}

task ApplyRecalibration {

  input {
    String recalibrated_vcf_filename
    File input_vcf
    File input_vcf_index
    File indels_recalibration
    File indels_recalibration_index
    File indels_tranches
    File snps_recalibration
    File snps_recalibration_index
    File snps_tranches
    Float indel_filter_level
    Float snp_filter_level
    Boolean use_allele_specific_annotations
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms5g \
      ApplyVQSR \
      -O tmp.indel.recalibrated.vcf \
      -V ~{input_vcf} \
      --recal-file ~{indels_recalibration} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      --tranches-file ~{indels_tranches} \
      --truth-sensitivity-filter-level ~{indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    gatk --java-options -Xms5g \
      ApplyVQSR \
      -O ~{recalibrated_vcf_filename} \
      -V tmp.indel.recalibrated.vcf \
      --recal-file ~{snps_recalibration} \
      ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
      --tranches-file ~{snps_tranches} \
      --truth-sensitivity-filter-level ~{snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  >>>

  runtime {
    memory: "7 GiB"
    cpu: "1"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File recalibrated_vcf = "~{recalibrated_vcf_filename}"
    File recalibrated_vcf_index = "~{recalibrated_vcf_filename}.tbi"
  }
}

task GatherVcfs {

  input {
    Array[File] input_vcfs
    String output_vcf_name
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    input_vcfs: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
    # This argument disables expensive checks that the file headers contain the same set of
    # genotyped samples and that files are in order by position of first record.
    gatk --java-options -Xms6g \
      GatherVcfsCloud \
      --ignore-safety-checks \
      --gather-type BLOCK \
      --input ~{sep=" --input " input_vcfs} \
      --output ~{output_vcf_name}

    tabix ~{output_vcf_name}
  >>>

  runtime {
    memory: "7 GiB"
    cpu: "1"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }
}

task SelectFingerprintSiteVariants {

  input {
    File input_vcf
    File haplotype_database
    String base_output_name
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    input_vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    function hdb_to_interval_list() {
      input=$1
      awk 'BEGIN{IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {print $1,$2,$2,"+","interval-"NR}' $1
    }

    hdb_to_interval_list ~{haplotype_database} > hdb.interval_list

    gatk --java-options -Xms6g \
      SelectVariants \
      --variant ~{input_vcf} \
      --intervals hdb.interval_list \
      --output ~{base_output_name}.vcf.gz
  >>>

  runtime {
    memory: "7.5 GiB"
    cpu: 1
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File output_vcf = "~{base_output_name}.vcf.gz"
    File output_vcf_index = "~{base_output_name}.vcf.gz.tbi"
  }
}

task CollectVariantCallingMetrics {

  input {
    File input_vcf
    File input_vcf_index
    String metrics_filename_prefix
    File dbsnp_vcf
    File dbsnp_vcf_index
    File interval_list
    File ref_dict
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options -Xms6g \
      CollectVariantCallingMetrics \
      --INPUT ~{input_vcf} \
      --DBSNP ~{dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ~{ref_dict} \
      --OUTPUT ~{metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ~{interval_list}
  >>>

  output {
    File detail_metrics_file = "~{metrics_filename_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "~{metrics_filename_prefix}.variant_calling_summary_metrics"
  }

  runtime {
    memory: "7.5 GiB"
    cpu: 2
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }
}

task GatherVariantCallingMetrics {

  input {
    Array[File] input_details
    Array[File] input_summaries
    String output_prefix
    Int disk_size
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    input_details: {
      localization_optional: true
    }
    input_summaries: {
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail

    input_details_fofn=~{write_lines(input_details)}
    input_summaries_fofn=~{write_lines(input_summaries)}

    # Jose says:
    # Cromwell will fall over if we have it try to localize tens of thousands of files,
    # so we manually localize files using gsutil.
    # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
    # PAPI doesn't do.

    # This is here to deal with the JES bug where commands may be run twice
    rm -rf metrics

    mkdir metrics
    RETRY_LIMIT=5

    count=0
    until cat $input_details_fofn | gsutil -m cp -L cp.log -c -I metrics/; do
      sleep 1
      ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
      echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    count=0
    until cat $input_summaries_fofn | gsutil -m cp -L cp.log -c -I metrics/; do
      sleep 1
      ((count++)) && ((count >= $RETRY_LIMIT)) && break
    done
    if [ "$count" -ge "$RETRY_LIMIT" ]; then
      echo 'Could not copy all the metrics from the cloud' && exit 1
    fi

    INPUT=$(cat $input_details_fofn | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("--INPUT metrics/%s ", $1)}')

    gatk --java-options -Xms2g \
      AccumulateVariantCallingMetrics \
      $INPUT \
      --OUTPUT ~{output_prefix}
  >>>

  runtime {
    memory: "3 GiB"
    cpu: "1"
    bootDiskSizeGb: 15
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 1
    docker: gatk_docker
  }

  output {
    File detail_metrics_file = "~{output_prefix}.variant_calling_detail_metrics"
    File summary_metrics_file = "~{output_prefix}.variant_calling_summary_metrics"
  }
}

task CrossCheckFingerprint {

  input {
    Array[File] gvcf_paths
    Array[File] vcf_paths
    File sample_name_map
    File haplotype_database
    String output_base_name
    Boolean scattered = false
    Array[String] expected_inconclusive_samples = []
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  parameter_meta {
    gvcf_paths: {
      localization_optional: true
    }
    vcf_paths: {
      localization_optional: true
    }
  }

  Int num_gvcfs = length(gvcf_paths)
  Int cpu = if num_gvcfs < 32 then num_gvcfs else 32
  # Compute memory to use based on the CPU count, following the pattern of
  # 3.75GiB / cpu used by GCP's pricing: https://cloud.google.com/compute/pricing
  Int memMb = round(cpu * 3.75 * 1024)
  Int disk = 100

  String output_name = output_base_name + ".fingerprintcheck"

  command <<<
    set -eu

    gvcfInputsList=~{write_lines(gvcf_paths)}
    vcfInputsList=~{write_lines(vcf_paths)}

    cp $gvcfInputsList gvcf_inputs.list
    cp $vcfInputsList vcf_inputs.list

    gatk --java-options -Xms~{memMb - 512}m \
      CrosscheckFingerprints \
      --INPUT gvcf_inputs.list \
      --SECOND_INPUT vcf_inputs.list \
      --HAPLOTYPE_MAP ~{haplotype_database} \
      --INPUT_SAMPLE_FILE_MAP ~{sample_name_map} \
      --CROSSCHECK_BY SAMPLE \
      --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
      --NUM_THREADS ~{cpu} \
      ~{true='--EXIT_CODE_WHEN_MISMATCH 0' false='' scattered} \
      --OUTPUT ~{output_name}

    if ~{scattered}; then
      # UNEXPECTED_MATCH is not possible with CHECK_SAME_SAMPLE
      matches=$(grep "EXPECTED_MATCH" ~{output_name} | wc -l)

      # check inconclusive samples
      expectedInconclusiveSamples=("~{sep='" "' expected_inconclusive_samples}")
      inconclusiveSamplesCount=0
      inconclusiveSamples=($(grep 'INCONCLUSIVE' ~{output_name} | cut -f 1))
      for sample in ${inconclusiveSamples[@]}; do
      if printf '%s\n' ${expectedInconclusiveSamples[@]} | grep -P '^'${sample}'$'; then
        inconclusiveSamplesCount=$((inconclusiveSamplesCount+1))
      fi
      done

      total_matches=$((inconclusiveSamplesCount + matches))
      if [[ ${total_matches} -eq ~{num_gvcfs} ]]; then
        >&2 echo "Found the correct number of matches (~{num_gvcfs}) for this shard"
      else
        >&2 echo "ERROR: Found $total_matches 'EXPECTED_MATCH' records, but expected ~{num_gvcfs}"
      exit 1
      fi
    fi
  >>>

  runtime {
    memory: memMb + " MiB"
    disks: "local-disk " + disk + " HDD"
    preemptible: 0
    docker: gatk_docker
  }

  output {
    File crosscheck_metrics = output_name
  }
}

task GatherPicardMetrics {

  input {
    Array[File] metrics_files
    String output_file_name
    Int disk_size
  }

  command {
    # Don't use this task to gather tens of thousands of files.
    # Cromwell can't handle it.

    # This cannot gather metrics with histograms

    head -n 7 ~{metrics_files[0]} > ~{output_file_name}

    for metrics_file in ~{sep=' ' metrics_files}; do
      sed -n '1,7d;p' $metrics_file | grep -v '^$' >> ~{output_file_name}
    done
  }

  output {
    File gathered_metrics = "~{output_file_name}"
  }

  runtime {
    cpu: 1
    memory: "3.75 GiB"
    preemptible: 1
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
  }
}

task GetFingerprintingIntervalIndices {

  input {
    Array[File] unpadded_intervals
    File haplotype_database
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.8.0"
  }

  command <<<
    set -xeo pipefail

    function rename_intervals(){
      interval_list=$1
      name=$2

      awk 'BEGIN{FS=IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {$5="'$name'"; print}' $interval_list
    }
    export -f rename_intervals

    function hdb_to_interval_list(){
      input=$1

      awk 'BEGIN{IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {print $1,$2,$2,"+","interval-"NR}' $1
    }

    function rename_scatter(){
      file=$1
      number=$(echo $file | sed -E 's|([0-9]+)-scattered\.interval.*|\1|')
      rename_intervals $file $number > scattered.renamed.$number.interval_list
    }
    export -f rename_scatter

    # rename the intervals within each interval_list according to the number in the name of the list

    cp ~{sep=' ' unpadded_intervals} ./

    cat ~{write_lines(unpadded_intervals)} | xargs -n1 basename | xargs -I{} bash -c 'rename_scatter $@' _ {}

    #find the first header
    find . -name "scattered.renamed.*.interval_list" | head -n1 | xargs cat | grep '^@' > all.interval_list

    # concatenate the resulting intervals (with no header)
    find . -name "scattered.renamed.*.interval_list"  | xargs cat | grep -v '^@' >> all.interval_list

    # convert the Haplotype_database to an interval_list
    hdb_to_interval_list ~{haplotype_database} > hdb.interval_list

    # find the intervals that overlap the haplotype_database
    gatk IntervalListTools \
      -ACTION OVERLAPS \
      -O all.sorted.interval_list \
      -I all.interval_list \
      -SI hdb.interval_list

    if grep -v '^@' all.sorted.interval_list; then
      grep -v '^@' all.sorted.interval_list | awk '{FS="\t"; print $5}' | uniq > indices.out
    else
      touch indices.out
    fi
  >>>

  output {
    Array[String] indices_to_fingerprint = read_lines("indices.out")
    File all_sorted_interval_list = "all.sorted.interval_list"
    File all_interval_list = "all.interval_list"
    File hdb_interval_list = "hdb.interval_list"
  }

  runtime {
    cpu: 2
    memory: "3.75 GiB"
    preemptible: 1
    bootDiskSizeGb: 15
    disks: "local-disk 10 HDD"
    docker: gatk_docker
  }
}

task PartitionSampleNameMap {

  input {
    File sample_name_map
    Int line_limit
  }

  command {

    cut -f 2 ~{sample_name_map} > sample_paths
    split -l ~{line_limit} -d sample_paths partition_

    # Let the OS catch up with creation of files for glob command
    sleep 1
  }

  output {
    Array[File] partitions = glob("partition_*")
  }

  runtime {
    memory: "1 GiB"
    preemptible: 1
    disks: "local-disk 10 HDD"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
  }
}