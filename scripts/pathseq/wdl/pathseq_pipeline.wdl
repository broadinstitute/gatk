########################################################################################################################
## PathSeq Pipeline WDL
########################################################################################################################
##
## Runs the PathSeq pipeline for detecting microbes in deep sequencing data.
##
## Important: Filtering metrics should only be generated (gather_filter_metrics=true) on a downsampled BAM.
##
## Important: For microbe-rich samples users should set downsampling=true and adjust downsample_reads so that the total
##   number of microbe reads is <10M. The total number of microbe reads can be estimated from a preliminary run of this
#    workflow with downsampling=true, filtering_only=false, and gather_filter_metrics=true.
##
## For further info see the GATK Documentation for the PathSeqPipelineSpark tool:
##   https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_spark_pathseq_PathSeqPipelineSpark.php
##
########################################################################################################################
##
## Input requirements :
## - Sequencing data in CRAM or BAM format
## - Host and microbe references files available in the GATK Resource Bundle (available on FTP):
##     https://software.broadinstitute.org/gatk/download/bundle
##
## - Input BAM files must comply with the following requirements:
## - - file must pass validation by ValidateSamFile
## - - all reads must have an RG tag
## - - one or more read groups all belong to a single sample (SM)
##
## Output:
## - non_host_paired_bam, non_host_unpaired_bam: BAM files containing quality-filtered non-host paired and unpaired reads
## - If filtering_only = false:
##   - final_bam : BAM file containing microbe-mapped reads and reads of unknown sequence
##   - taxonomy_scores : tab-separated value (.tsv) file of taxonomic abundance scores
##   - score_metrics_file : Picard-style metrics file for scoring stage
##   - non_host_mapped_reads : number of non-host reads mapped to microbe reference
##   - non_host_unmapped_reads : number of non-host reads that did not map to the microbe reference
## - If gather_filter_metrics = true:
##   - filter_metrics_file : Picard-style metrics file for filtering stage
##   - frac_after_prealigned_filter : estimated fraction of reads remaining after filtering prealigned host reads
##   - frac_after_qual_cpx_filter : estimated fraction of reads remaining after low-quality and low-complexity filtering
##   - frac_after_host_filter : estimated fraction of reads remaining after host read filtering
##   - frac_after_dedup : estimated fraction of reads remaining after deduplication
##   - frac_final_paired : estimated fraction of reads that were paired in the final output
##   - frac_final_unpaired : estimated fraction of reads that were unpaired in the final output
##   - frac_final_total : estimated fraction of reads in the final output
##   - frac_qual_cpx_filtered : estimated fraction of reads removed by low-quality/low-complexity filtering
##   - frac_host_filtered : estimated fraction of reads removed by host filtering
##   - frac_dup_filtered : estimated fraction of reads removed by deduplication
## - If downsample = true:
##   - total_reads_if_avail : total number of reads in the original input BAM
##
########################################################################################################################

workflow PathSeqPipeline {

  String sample_name
  File input_bam_or_cram
  # Required if a cram
  File? input_bam_or_cram_index

  # Required if the input is a cram
  File? cram_reference_fasta
  File? cram_reference_fasta_index
  File? cram_reference_dict

  # Set to true if host aligned. WARNING: Results in loss of EBV reads if in the aligned reference.
  Boolean is_host_aligned

  File? kmer_file
  File? filter_bwa_image
  # Required if filtering_only=true
  File? microbe_bwa_image
  File? microbe_dict
  File? taxonomy_file

  # If enabled, filter metrics will be estimated using a downsampled bam with this many reads (recommended)
  # If disabled, no filter metrics will be generated
  Boolean downsample = false
  Int downsample_reads = 1000000

  # Enable to only perform host filtering
  Boolean filtering_only = false

  # Only recommended if downsample = true
  Boolean gather_filter_metrics = false

  # This can be calculated from a downsampled run to help optimize disk allocation during filtering
  Float frac_non_host_reads = 1.0

  # Filtering options
  Array[String]? host_contigs_to_ignore
  Boolean? filter_duplicates
  Boolean? skip_quality_filters
  Int? host_min_identity
  Int? host_kmer_threshold
  Int? filter_bwa_seed_length
  Int? min_clipped_read_length
  Int? max_masked_bases
  Int? min_base_quality
  Int? quality_threshold
  Int? dust_mask_quality
  Int? dust_window
  Float? dust_t
  Int? max_adapter_mistmatches
  Int? min_adapter_length
  Int? filter_reads_per_partition
  Int? filter_bam_partition_size
  Boolean? skip_pre_bwa_repartition

  # Alignment options
  Int? microbe_min_seed_length
  Int? max_alternate_hits
  Int? bwa_score_threshold

  # Taxonomy scoring options
  Boolean? divide_by_genome_length
  Float? min_score_identity
  Float? identity_margin
  Boolean? not_normalized_by_kingdom
  Int? score_reads_per_partition_estimate

  # Runtime parameters
  String gatk_docker
  String samtools_docker = "biocontainers/samtools:v1.9-4-deb_cv1"
  String genomes_in_the_cloud_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1512499786"
  String linux_docker = "ubuntu:18.10"

  File? gatk4_jar_override

  Int? cram_to_bam_preemptible_attempts
  Int? downsample_preemptible_attempts
  Int? filter_preemptible_attempts
  Int? align_preemptible_attempts
  Int? score_preemptible_attempts
  Int? process_filter_metrics_preemptible_attempts
  Int? process_score_metrics_preemptible_attempts

  Int? cram_to_bam_cpu
  Int? filter_cpu
  Int? align_cpu
  Int? score_cpu

  Float? cram_to_bam_mem_gb
  Int? filter_mem_gb
  Int? align_mem_gb
  Int? score_mem_gb

  Boolean? filter_ssd
  Boolean? align_ssd
  Boolean? score_ssd

  Int? cram_to_bam_max_retries

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? cram_to_bam_additional_disk_gb
  Int? downsample_additional_disk_gb
  Int? filter_additional_disk_gb
  Int? align_additional_disk_gb
  Int? score_additional_disk_gb

  Boolean is_bam = basename(input_bam_or_cram, ".bam") + ".bam" == basename(input_bam_or_cram)

  # Convert to BAM if we have a CRAM
  if (!is_bam) {
    call CramToBam {
      input:
        cram_file = input_bam_or_cram,
        cram_index = select_first([input_bam_or_cram_index]),
        reference_fasta = select_first([cram_reference_fasta]),
        reference_index = select_first([cram_reference_fasta_index]),
        docker = samtools_docker,
        cpu = cram_to_bam_cpu,
        mem_gb = cram_to_bam_mem_gb,
        extra_disk_gb = cram_to_bam_additional_disk_gb,
        preemptible = cram_to_bam_preemptible_attempts,
        max_retries = cram_to_bam_max_retries,
    }
  }

  File bam_file = select_first([CramToBam.bam_file, input_bam_or_cram])
  File? bam_index = if defined(CramToBam.bam_index) then CramToBam.bam_index else input_bam_or_cram_index

  # Downsample bam for filter metrics estimation
  if (downsample) {
    call Downsample {
      input:
        input_bam_file=bam_file,
        input_bam_index_file=bam_index,
        downsampled_bam_filename="${sample_name}.downsampled.bam",
        reads_after_downsampling=downsample_reads,
        additional_disk_gb=downsample_additional_disk_gb,
        preemptible_tries=downsample_preemptible_attempts,
        docker=genomes_in_the_cloud_docker
    }
  }

  call PathSeqFilter {
    input:
      sample_name=sample_name,
      input_bam_or_cram=select_first([Downsample.output_bam_file, bam_file]),
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      frac_non_host_reads=frac_non_host_reads,
      gather_metrics=gather_filter_metrics,
      is_host_aligned=is_host_aligned,
      host_contigs_to_ignore=host_contigs_to_ignore,
      filter_duplicates=filter_duplicates,
      skip_quality_filters=skip_quality_filters,
      min_clipped_read_length=min_clipped_read_length,
      bam_partition_size=filter_bam_partition_size,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      max_masked_bases=max_masked_bases,
      min_base_quality=min_base_quality,
      quality_threshold=quality_threshold,
      dust_mask_quality=dust_mask_quality,
      dust_window=dust_window,
      dust_t=dust_t,
      host_kmer_threshold=host_kmer_threshold,
      max_adapter_mistmatches=max_adapter_mistmatches,
      min_adapter_length=min_adapter_length,
      filter_reads_per_partition=filter_reads_per_partition,
      skip_pre_bwa_repartition=skip_pre_bwa_repartition,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=filter_mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=filter_preemptible_attempts,
      additional_disk_gb=filter_additional_disk_gb,
      cpu=filter_cpu,
      use_ssd=filter_ssd
  }

  if (gather_filter_metrics) {
    call ProcessFilterMetrics {
      input:
        metrics_file=PathSeqFilter.filter_metrics,
        preemptible_tries=process_filter_metrics_preemptible_attempts,
        docker=linux_docker
    }
  }

  if (!filtering_only) {
    call PathSeqAlign {
      input:
        sample_name=sample_name,
        input_paired_bam=PathSeqFilter.paired_bam_out,
        input_unpaired_bam=PathSeqFilter.unpaired_bam_out,
        microbe_bwa_image=select_first([microbe_bwa_image]),
        microbe_dict=select_first([microbe_dict]),
        microbe_min_seed_length=microbe_min_seed_length,
        max_alternate_hits=max_alternate_hits,
        bwa_score_threshold=bwa_score_threshold,
        gatk4_jar_override=gatk4_jar_override,
        mem_gb=align_mem_gb,
        gatk_docker=gatk_docker,
        preemptible_attempts=align_preemptible_attempts,
        additional_disk_gb=align_additional_disk_gb,
        cpu=align_cpu,
        use_ssd=align_ssd
    }

    call PathSeqScore {
      input:
        sample_name=sample_name,
        input_paired_bam=PathSeqAlign.paired_bam_out,
        input_unpaired_bam=PathSeqAlign.unpaired_bam_out,
        taxonomy_file=select_first([taxonomy_file]),
        divide_by_genome_length=divide_by_genome_length,
        min_score_identity=min_score_identity,
        identity_margin=identity_margin,
        not_normalized_by_kingdom=not_normalized_by_kingdom,
        score_reads_per_partition_estimate=score_reads_per_partition_estimate,
        gatk4_jar_override=gatk4_jar_override,
        mem_gb=score_mem_gb,
        gatk_docker=gatk_docker,
        preemptible_attempts=score_preemptible_attempts,
        additional_disk_gb=score_additional_disk_gb,
        cpu=score_cpu,
        use_ssd=score_ssd
    }

    call ProcessScoreMetrics {
      input:
        metrics_file=PathSeqScore.score_metrics,
        preemptible_tries=process_score_metrics_preemptible_attempts,
        docker=linux_docker
    }
  }

  output {
    # Filtered non-host reads
    File non_host_paired_bam = PathSeqFilter.paired_bam_out
    File non_host_unpaired_bam = PathSeqFilter.unpaired_bam_out

    # Only generated if filtering_only=false
    File? final_bam = PathSeqScore.bam_out
    File? taxonomy_scores = PathSeqScore.scores
    File? score_metrics_file = PathSeqScore.score_metrics
    Int? non_host_mapped_reads = ProcessScoreMetrics.non_host_mapped_reads
    Int? non_host_unmapped_reads = ProcessScoreMetrics.non_host_unmapped_reads

    # Only generated if gather_filter_metrics=true
    File? filter_metrics_file = ProcessFilterMetrics.metrics_file_out
    Float? frac_after_prealigned_filter = ProcessFilterMetrics.frac_after_prealigned_filter
    Float? frac_after_qual_cpx_filter = ProcessFilterMetrics.frac_after_qual_cpx_filter
    Float? frac_after_host_filter = ProcessFilterMetrics.frac_after_host_filter
    Float? frac_after_dedup = ProcessFilterMetrics.frac_after_dedup
    Float? frac_non_host_paired = ProcessFilterMetrics.frac_final_paired
    Float? frac_non_host_unpaired = ProcessFilterMetrics.frac_final_unpaired
    Float? frac_non_host_total =ProcessFilterMetrics.frac_final_total
    Float? frac_qual_cpx_filtered = ProcessFilterMetrics.frac_qual_cpx_filtered
    Float? frac_host_filtered = ProcessFilterMetrics.frac_host_filtered
    Float? frac_dup_filtered = ProcessFilterMetrics.frac_dup_filtered

    # Only generated if downsample=true
    Int? total_reads = Downsample.total_reads
  }
}

task CramToBam {
  File cram_file
  File? cram_index
  File reference_fasta
  File? reference_index

  String docker

  Int? cpu = 4
  Float? mem_gb = 15
  Int? extra_disk_gb = 10
  Int? preemptible = 3
  Int? max_retries = 1

  String bam_file_name = basename(cram_file, ".cram") + ".bam"

  File cram_index_file = select_first([cram_index, cram_file + ".crai"])
  File reference_index_file = select_first([reference_index, reference_fasta + ".fai"])

  Float cram_inflate_ratio = 3.0
  Float cram_size = size(cram_file, "GiB")
  Float cram_index_size = size(cram_index_file, "GiB")
  Float bam_size = cram_inflate_ratio * cram_size
  Float bam_index_size = cram_index_size
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(reference_index_file, "GiB")
  Int vm_disk_size = ceil(cram_size + cram_index_size + bam_size + bam_index_size + ref_size + ref_index_size + extra_disk_gb)

  output {
    File bam_file = bam_file_name
    File bam_index = bam_file_name + ".bai"
  }
  command <<<

        set -Eeuo pipefail

        # covert cram to bam
        samtools view  -@ ${cpu} -h -b -T "${reference_fasta}" -o "${bam_file_name}" "${cram_file}"

        # index bam file
        samtools index -@ ${cpu} "${bam_file_name}"

  >>>
  runtime {
    cpu: 1
    memory: mem_gb + " GiB"
    disks: "local-disk " + vm_disk_size + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: preemptible
    maxRetries: max_retries
  }
}

# Downsamples BAM to a specified number of reads
task Downsample {
  File input_bam_file
  File? input_bam_index_file
  String downsampled_bam_filename

  Int reads_after_downsampling

  String docker
  Int additional_disk_gb = 20
  Int preemptible_tries = 3

  Int disk_size = ceil(size(input_bam_file, "GB")*2 + additional_disk_gb)

  command <<<
    set -euo pipefail
    if ${defined(input_bam_index_file)}; then
      NUM_READS=`samtools idxstats ${input_bam_file} | awk '{s+=$3+$4} END {print s}'`
    else
      NUM_READS=`samtools view -c ${input_bam_file}`
    fi
    P_DOWNSAMPLE=`python -c "print ${reads_after_downsampling}/float($NUM_READS)"`
    RESULT=`python -c "print $P_DOWNSAMPLE > 1"`
    if [ "$RESULT" == "True" ]
    then
        P_DOWNSAMPLE="1"
    fi
    echo $NUM_READS > num_reads.txt
    java -Xmx2000m -jar /usr/gitc/picard.jar \
      DownsampleSam \
      INPUT=${input_bam_file} \
      OUTPUT=${downsampled_bam_filename} \
      P=$P_DOWNSAMPLE \
      VALIDATION_STRINGENCY=SILENT
  >>>
  output {
    File output_bam_file = "${downsampled_bam_filename}"
    Int total_reads = read_int("num_reads.txt")
  }
  runtime {
    preemptible: "${preemptible_tries}"
    docker: docker
    memory: "3.75 GiB"
    cpu: "1"
    disks: "local-disk ${disk_size} HDD"
  }
}

task PathSeqFilter {

  # Inputs for this task
  String sample_name
  File input_bam_or_cram

  File? kmer_file
  File? filter_bwa_image

  Boolean is_host_aligned
  Array[String]? host_contigs_to_ignore
  Boolean gather_metrics
  # Optimizes disk space if provided
  Float frac_non_host_reads = 1.0

  Boolean? skip_quality_filters
  Boolean? skip_pre_bwa_repartition
  Boolean? filter_duplicates
  Int? host_min_identity
  Int? filter_bwa_seed_length
  Int? min_clipped_read_length
  Int? bam_partition_size
  Int? max_masked_bases
  Int? min_base_quality
  Int? quality_threshold
  Int? dust_mask_quality
  Int? dust_window
  Float? dust_t
  Int? host_kmer_threshold
  Int? max_adapter_mistmatches
  Int? min_adapter_length
  Int? filter_reads_per_partition

  String paired_bam_output_path = "${sample_name}.non_host.paired.bam"
  String unpaired_bam_output_path = "${sample_name}.non_host.unpaired.bam"
  String filter_metrics_output_path = "${sample_name}.pathseq.filter_metrics"

  File? gatk4_jar_override

  # Default to WARNING which will avoid excessive Spark logging
  String verbosity = "WARNING"

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 32
  Int preemptible_attempts = 3
  Float additional_disk_gb = 50
  Int cpu = 8
  Boolean use_ssd = false

  Int disk_size = ceil(((1.0 + frac_non_host_reads)*size(input_bam_or_cram, "GB")) + size(kmer_file, "GB") + size(filter_bwa_image, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil((machine_mem - size(filter_bwa_image, "MB")) * 0.8)

  command <<<
    set -euo pipefail
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    touch ${filter_metrics_output_path}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqFilterSpark \
      --input ${input_bam_or_cram} \
      --paired-output ${paired_bam_output_path} \
      --unpaired-output ${unpaired_bam_output_path} \
      --is-host-aligned ${is_host_aligned}  \
      --verbosity ${verbosity} \
      ${if defined(host_contigs_to_ignore) then "--ignore-alignment-contigs " else ""} ${sep=" --ignore-alignment-contigs " host_contigs_to_ignore} \
      ${if gather_metrics then "--filter-metrics ${filter_metrics_output_path}" else ""} \
      ${if defined(kmer_file) then "--kmer-file ${kmer_file}" else ""} \
      ${if defined(filter_bwa_image) then "--filter-bwa-image ${filter_bwa_image}" else ""} \
      ${if defined(bam_partition_size) then "--bam-partition-size ${bam_partition_size}" else ""} \
      ${if defined(skip_quality_filters) then "--skip-quality-filters ${skip_quality_filters}" else ""}\
      ${if defined(min_clipped_read_length) then "--min-clipped-read-length ${min_clipped_read_length}" else ""} \
      ${if defined(max_masked_bases) then "--max-masked-bases ${max_masked_bases}" else ""} \
      ${if defined(min_base_quality) then "--min-base-quality ${min_base_quality}" else ""} \
      ${if defined(quality_threshold) then "--quality-threshold ${quality_threshold}" else ""} \
      ${if defined(dust_mask_quality) then "--dust-mask-quality ${dust_mask_quality}" else ""} \
      ${if defined(dust_window) then "--dust-window ${dust_window}" else ""} \
      ${if defined(dust_t) then "--dust-t ${dust_t}" else ""} \
      ${if defined(host_kmer_threshold) then "--host-kmer-thresh ${host_kmer_threshold}" else ""} \
      ${if defined(max_adapter_mistmatches) then "--max-adapter-mismatches ${max_adapter_mistmatches}" else ""} \
      ${if defined(min_adapter_length) then "--min-adapter-length ${min_adapter_length}" else ""} \
      ${if defined(filter_bwa_seed_length) then "--filter-bwa-seed-length ${filter_bwa_seed_length}" else ""} \
      ${if defined(host_min_identity) then "--host-min-identity ${host_min_identity}" else ""} \
      ${if defined(filter_duplicates) then "--filter-duplicates ${filter_duplicates}" else ""} \
      ${if defined(skip_pre_bwa_repartition) then "--skip-pre-bwa-repartition ${skip_pre_bwa_repartition}" else ""} \
      ${if defined(filter_reads_per_partition) then "--filter-reads-per-partition ${filter_reads_per_partition}" else ""}

    if [ ! -f "${paired_bam_output_path}" ]; then
    	echo "File ${paired_bam_output_path} not found, creating empty BAM"
    	printf "@HD     VN:1.5  SO:queryname\n@RG     ID:A    SM:${sample_name}\n" | samtools view -hb - > ${paired_bam_output_path}
    fi

    if [ ! -f "${unpaired_bam_output_path}" ]; then
    	echo "File ${unpaired_bam_output_path} not found, creating empty BAM"
    	printf "@HD     VN:1.5  SO:queryname\n@RG     ID:A    SM:${sample_name}\n" | samtools view -hb - > ${unpaired_bam_output_path}
    fi
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File paired_bam_out = "${paired_bam_output_path}"
    File unpaired_bam_out = "${unpaired_bam_output_path}"
    File filter_metrics = "${sample_name}.pathseq.filter_metrics"
  }
}

task PathSeqAlign {

  # Inputs for this task
  String sample_name
  File input_paired_bam
  File input_unpaired_bam

  File microbe_bwa_image
  File microbe_dict

  Int? microbe_min_seed_length
  Int? max_alternate_hits
  Int? bwa_score_threshold

  String paired_bam_output_path = "${sample_name}.microbe_aligned.paired.bam"
  String unpaired_bam_output_path = "${sample_name}.microbe_aligned.unpaired.bam"

  File? gatk4_jar_override

  # Default to WARNING which will avoid excessive Spark logging
  String verbosity = "WARNING"

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 140
  Int preemptible_attempts = 3
  Float additional_disk_gb = 50
  Int cpu = 8
  Boolean use_ssd = false

  Int disk_size = ceil((2.5*size(input_paired_bam, "GB")) + (2.5*size(input_unpaired_bam, "GB")) + size(microbe_bwa_image, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil((machine_mem - size(microbe_bwa_image, "MB")) * 0.8)

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqBwaSpark \
      --paired-input ${input_paired_bam} \
      --unpaired-input ${input_unpaired_bam} \
      --paired-output ${paired_bam_output_path} \
      --unpaired-output ${unpaired_bam_output_path} \
      --microbe-bwa-image ${microbe_bwa_image} \
      --microbe-dict ${microbe_dict} \
      --verbosity ${verbosity} \
      ${if defined(microbe_min_seed_length) then "--microbe-min-seed-length ${microbe_min_seed_length}" else ""} \
      ${if defined(max_alternate_hits) then "--max-alternate-hits ${max_alternate_hits}" else ""} \
      ${if defined(bwa_score_threshold) then "--bwa-score-threshold ${bwa_score_threshold}" else ""}
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File paired_bam_out = "${paired_bam_output_path}"
    File unpaired_bam_out = "${unpaired_bam_output_path}"
  }
}

task PathSeqScore {

  # Inputs for this task
  String sample_name
  File input_paired_bam
  File input_unpaired_bam

  File taxonomy_file

  Boolean? divide_by_genome_length
  Float? min_score_identity
  Float? identity_margin
  Boolean? not_normalized_by_kingdom
  Int? score_reads_per_partition_estimate

  String bam_output_path = "${sample_name}.pathseq.bam"
  String scores_output_path = "${sample_name}.pathseq.tsv"
  String score_metrics_output_path = "${sample_name}.pathseq.score_metrics"

  File? gatk4_jar_override

  # Runtime parameters
  String gatk_docker
  Int mem_gb = 8
  Int preemptible_attempts = 3
  Float additional_disk_gb = 10
  Int cpu = 2
  Boolean use_ssd = false

  Int disk_size = ceil(size(input_paired_bam, "GB") + size(input_unpaired_bam, "GB") + additional_disk_gb)

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = mem_gb * 1000
  Int command_mem = ceil(machine_mem * 0.8)

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
    gatk --java-options "-Xmx${command_mem}m" \
      PathSeqScoreSpark \
      --paired-input ${input_paired_bam} \
      --unpaired-input ${input_unpaired_bam} \
      --output ${bam_output_path} \
      --scores-output ${scores_output_path} \
      --score-metrics ${score_metrics_output_path} \
      --taxonomy-file ${taxonomy_file} \
      ${if defined(min_score_identity) then "--min-score-identity ${min_score_identity}" else ""} \
      ${if defined(identity_margin) then "--identity-margin ${identity_margin}" else ""} \
      ${if defined(divide_by_genome_length) then "--divide-by-genome-length ${divide_by_genome_length}" else ""} \
      ${if defined(not_normalized_by_kingdom) then "--not-normalized-by-kingdom ${not_normalized_by_kingdom}" else ""} \
      ${if defined(score_reads_per_partition_estimate) then "--score-reads-per-partition-estimate ${score_reads_per_partition_estimate}" else ""}
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_size + if use_ssd then " SSD" else " HDD"
    preemptible: preemptible_attempts
    cpu: cpu
  }
  output {
    File bam_out = "${sample_name}.pathseq.bam"
    File scores = "${sample_name}.pathseq.tsv"
    File score_metrics = "${sample_name}.pathseq.score_metrics"
  }
}

# Extracts filter metrics as fraction of total reads
task ProcessFilterMetrics {
  File metrics_file
  String docker
  Int preemptible_tries = 3

  String dollar = "$"

  command <<<
    set -euo pipefail
    if [ ! -s ${metrics_file} ]; then
      echo "" | awk -v OFS="\t" '{print 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}'
    else
      sed -n '8p' ${metrics_file} \
        | awk -F"\t" -v OFS="\t" '{if (${dollar}1 > 0) {print ${dollar}2/${dollar}1, ${dollar}3/${dollar}1, ${dollar}4/${dollar}1, ${dollar}5/${dollar}1, ${dollar}6/${dollar}1, ${dollar}7/${dollar}1, ${dollar}8/${dollar}1, ${dollar}9/${dollar}1, ${dollar}0/${dollar}1, ${dollar}11/${dollar}1} else {print 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}'
    fi > metrics.txt

    cut -f1 metrics.txt > frac_after_prealigned_filter.txt
    cut -f2 metrics.txt> frac_after_qual_cpx_filter.txt
    cut -f3 metrics.txt > frac_after_host_filter.txt
    cut -f4 metrics.txt > frac_after_dedup.txt
    cut -f5 metrics.txt > frac_final_paired.txt
    cut -f6 metrics.txt > frac_final_unpaired.txt
    cut -f7 metrics.txt > frac_final_total.txt
    cut -f8 metrics.txt > frac_qual_cpx_filtered.txt
    cut -f9 metrics.txt > frac_host_filtered.txt
    cut -f10 metrics.txt > frac_dup_filtered.txt
  >>>
  output {
    File metrics_file_out = metrics_file
    Float frac_after_prealigned_filter = read_float("frac_after_prealigned_filter.txt")
    Float frac_after_qual_cpx_filter = read_float("frac_after_qual_cpx_filter.txt")
    Float frac_after_host_filter = read_float("frac_after_host_filter.txt")
    Float frac_after_dedup = read_float("frac_after_dedup.txt")
    Float frac_final_paired = read_float("frac_final_paired.txt")
    Float frac_final_unpaired = read_float("frac_final_unpaired.txt")
    Float frac_final_total = read_float("frac_final_total.txt")
    Float frac_qual_cpx_filtered = read_float("frac_qual_cpx_filtered.txt")
    Float frac_host_filtered = read_float("frac_host_filtered.txt")
    Float frac_dup_filtered = read_float("frac_dup_filtered.txt")
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker
    memory: "2 GiB"
    cpu: "1"
    disks: "local-disk 10 HDD"
  }
}

# Extracts score metrics
task ProcessScoreMetrics {
  File metrics_file
  String docker
  Int preemptible_tries = 3

  command <<<
    set -euo pipefail
    sed -n '8p' ${metrics_file} > metrics.txt
    cut -f1 metrics.txt > non_host_mapped_reads.txt
    cut -f2 metrics.txt > non_host_unmapped_reads.txt
  >>>
  output {
    Int non_host_mapped_reads = read_int("non_host_mapped_reads.txt")
    Int non_host_unmapped_reads = read_int("non_host_unmapped_reads.txt")
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker
    memory: "2 GiB"
    cpu: "1"
    disks: "local-disk 10 HDD"
  }
}
