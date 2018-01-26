########################################################################################################################
## PathSeq Pipeline WDL
########################################################################################################################
##
## Runs the PathSeq pipeline
##
## For further info see the GATK Documentation for the PathSeqPipelineSpark tool:
##   https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_spark_pathseq_PathSeqPipelineSpark.php
##
########################################################################################################################
##
## Input requirements :
## - Sequencing data in BAM format
## - Host and microbe references files available in the GATK Resource Bundle (available on FTP):
##     https://software.broadinstitute.org/gatk/download/bundle
##
## - Input BAM files must comply with the following requirements:
## - - file must pass validation by ValidateSamFile
## - - all reads must have an RG tag
## - - one or more read groups all belong to a single sample (SM)
##
## Output:
## - BAM file containing microbe-mapped reads and reads of unknown sequence
## - Tab-separated value (.tsv) file of taxonomic abundance scores
## - Picard-style metrics files for the filter and scoring phases of the pipeline
##
########################################################################################################################

task PathseqPipeline {

  # Inputs for this task
  String sample_name
  File input_bam

  File kmer_file
  File filter_bwa_image
  File microbe_bwa_image
  File microbe_fasta
  File microbe_fasta_dict
  File taxonomy_file

  Boolean is_host_aligned
  Boolean? skip_quality_filters
  Boolean? filter_duplicates
  Boolean? skip_pre_bwa_repartition
  Boolean? divide_by_genome_length
  Int? filter_bwa_seed_length
  Int? host_min_identity
  Int? min_clipped_read_length
  Int? bam_partition_size
  Float? min_score_identity
  Float? identity_margin

  String bam_output_path = "${sample_name}.pathseq.bam"
  String scores_output_path = "${sample_name}.pathseq.tsv"
  String filter_metrics_output_path = "${sample_name}.pathseq.filter_metrics"
  String score_metrics_output_path = "${sample_name}.pathseq.score_metrics"

  File? gatk4_jar_override

  # Runtime parameters
  Int? mem_gb
  String gatk_docker
  Int? preemptible_attempts
  Int? disk_space_gb
  Int? cpu
  Boolean use_ssd = true

  # You may have to change the following two parameter values depending on the task requirements
  Int default_ram_mb = 208000
  # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
  Int default_disk_space_gb = 400
  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
  Int command_mem = machine_mem - 4000

  Float default_min_score_identity = 0.9
  Float default_identity_margin = 0.02

  command <<<
    set -e
    export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
      gatk --java-options "-Xmx${command_mem}m" \
      PathSeqPipelineSpark \
      --input ${input_bam} \
      --output ${bam_output_path} \
      --scores-output ${scores_output_path} \
      --filter-metrics ${filter_metrics_output_path} \
      --score-metrics ${score_metrics_output_path} \
      --kmer-file ${kmer_file} \
      --filter-bwa-image ${filter_bwa_image} \
      --microbe-bwa-image ${microbe_bwa_image} \
      --microbe-fasta ${microbe_fasta} \
      --taxonomy-file ${taxonomy_file} \
      --bam-partition-size ${select_first([bam_partition_size, 4000000])} \
      --is-host-aligned ${is_host_aligned} \
      --skip-quality-filters ${select_first([skip_quality_filters, false])} \
      --min-clipped-read-length ${select_first([min_clipped_read_length, 60])} \
      --filter-bwa-seed-length ${select_first([filter_bwa_seed_length, 19])} \
      --host-min-identity ${select_first([host_min_identity, 30])} \
      --filter-duplicates ${select_first([filter_duplicates, true])} \
      --skip-pre-bwa-repartition ${select_first([skip_pre_bwa_repartition, false])} \
      --min-score-identity ${select_first([min_score_identity, default_min_score_identity])} \
      --identity-margin ${select_first([identity_margin, default_identity_margin])} \
      --divide-by-genome-length ${select_first([divide_by_genome_length, true])}
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem + " MB"
    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
    cpu: select_first([cpu, 32])
  }
  output {
    File bam_output_pathFile = "${sample_name}.pathseq.bam"
    File outputScoresFile = "${sample_name}.pathseq.tsv"
    File outputFilterMetricsFile = "${sample_name}.pathseq.filter_metrics"
    File outputScoreMetricsFile = "${sample_name}.pathseq.score_metrics"
  }
}

workflow PathSeqPipelineWorkflow {

  String sample_name
  File input_bam

  File kmer_file
  File filter_bwa_image
  File microbe_bwa_image
  File microbe_fasta
  File microbe_fasta_dict
  File taxonomy_file

  Boolean is_host_aligned
  Boolean? filter_duplicates
  Boolean? skip_pre_bwa_repartition
  Boolean? divide_by_genome_length
  Int? filter_bwa_seed_length
  Int? host_min_identity
  Int? min_clipped_read_length
  Int? bam_partition_size
  Float? min_score_identity
  Float? identity_margin

  File? gatk4_jar_override

  # Runtime parameters
  Int? mem_gb
  String gatk_docker
  Int? preemptible_attempts
  Int? cpu

  # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
  Int? increase_disk_size

  # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
  # Cromwell error from asking for 0 disk when the input is less than 1GB.
  # Also Spark requires some temporary storage.
  Int additional_disk = select_first([increase_disk_size, 20])

  # Disk sizes for Downsample and PathSeq tasks
  Float disk_space_gb = size(input_bam, "GB") + size(kmer_file, "GB") + size(filter_bwa_image, "GB") + size(microbe_bwa_image, "GB") + size(microbe_fasta, "GB") + additional_disk

  call PathseqPipeline {
    input:
      sample_name=sample_name,
      input_bam=input_bam,
      kmer_file=kmer_file,
      filter_bwa_image=filter_bwa_image,
      microbe_bwa_image=microbe_bwa_image,
      microbe_fasta=microbe_fasta,
      microbe_fasta_dict=microbe_fasta_dict,
      taxonomy_file=taxonomy_file,
      is_host_aligned=is_host_aligned,
      filter_duplicates=filter_duplicates,
      skip_pre_bwa_repartition=skip_pre_bwa_repartition,
      divide_by_genome_length=divide_by_genome_length,
      min_clipped_read_length=min_clipped_read_length,
      bam_partition_size=bam_partition_size,
      min_score_identity=min_score_identity,
      identity_margin=identity_margin,
      host_min_identity=host_min_identity,
      filter_bwa_seed_length=filter_bwa_seed_length,
      gatk4_jar_override=gatk4_jar_override,
      mem_gb=mem_gb,
      gatk_docker=gatk_docker,
      preemptible_attempts=preemptible_attempts,
      disk_space_gb=disk_space_gb,
      cpu=cpu
  }
  output {
    File pathseqBam = PathseqPipeline.bam_output_pathFile
    File pathseqScores = PathseqPipeline.outputScoresFile
    File pathseqFilterMetrics = PathseqPipeline.outputFilterMetricsFile
    File pathseqScoreMetrics = PathseqPipeline.outputScoreMetricsFile
  }
}