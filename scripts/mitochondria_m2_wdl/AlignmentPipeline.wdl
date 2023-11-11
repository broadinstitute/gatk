version 1.0

workflow AlignmentPipeline {

  meta {
    description: "Uses BWA to align unmapped bam and marks duplicates."
  }

  input {
    File input_bam
    String basename = basename(input_bam, ".bam")

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa

    #Optional runtime arguments
    Int? preemptible_tries
    Int disk_size = ceil(size(input_bam, "GB") * 4 + size(mt_fasta, "GB") + size(mt_fasta_index, "GB") + size(mt_amb, "GB") + size(mt_ann, "GB") + size(mt_bwt, "GB") + size(mt_pac, "GB") + size(mt_sa, "GB")) + 20
    Int mem = 6
  }

  parameter_meta {
    input_bam: "Input is an unaligned subset bam of NuMT and chrM reads and their mates. All reads must be paired."
    mt_aligned_bam: "Output is aligned duplicate marked coordinate sorted bam."
  }

  call GetBwaVersion

  call AlignAndMarkDuplicates {
    input:
      input_bam = input_bam,
      bwa_version = GetBwaVersion.version,
      output_bam_basename = basename + ".realigned",
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      ref_amb = mt_amb,
      ref_ann = mt_ann,
      ref_bwt = mt_bwt,
      ref_pac = mt_pac,
      ref_sa = mt_sa,
      preemptible_tries = preemptible_tries,
      disk_size = disk_size,
      mem = mem
  }

  output {
    File mt_aligned_bam = AlignAndMarkDuplicates.output_bam
    File mt_aligned_bai = AlignAndMarkDuplicates.output_bam_index
    File duplicate_metrics = AlignAndMarkDuplicates.duplicate_metrics
  }
}

task AlignAndMarkDuplicates {
  input {
    File input_bam
    String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 2 -Y $bash_ref_fasta"
    String bwa_version
    String output_bam_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    String? read_name_regex

    # runtime
    Int? preemptible_tries
    Int disk_size
    Int mem
  }

  String basename = basename(input_bam, ".bam")
  String metrics_filename = basename + ".metrics"

  meta {
    description: "Aligns with BWA and MergeBamAlignment, then Marks Duplicates. Outputs a coordinate sorted bam."
  }
  parameter_meta {
    input_bam: "Unmapped bam"
    bwa_version: "BWA version to be added to header of aligned bam"
  }
  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{ref_fasta}
    java -Xms5000m -jar /usr/gitc/picard.jar \
      SamToFastq \
      INPUT=~{input_bam} \
      FASTQ=/dev/stdout \
      INTERLEAVE=true \
      NON_PF=true | \
    /usr/gitc/~{bwa_commandline} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    java -Xms3000m -jar /usr/gitc/picard.jar \
      MergeBamAlignment \
      VALIDATION_STRINGENCY=SILENT \
      EXPECTED_ORIENTATIONS=FR \
      ATTRIBUTES_TO_RETAIN=X0 \
      ATTRIBUTES_TO_REMOVE=NM \
      ATTRIBUTES_TO_REMOVE=MD \
      ALIGNED_BAM=/dev/stdin \
      UNMAPPED_BAM=~{input_bam} \
      OUTPUT=mba.bam \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      PAIRED_RUN=true \
      SORT_ORDER="unsorted" \
      IS_BISULFITE_SEQUENCE=false \
      ALIGNED_READS_ONLY=false \
      CLIP_ADAPTERS=false \
      MAX_RECORDS_IN_RAM=2000000 \
      ADD_MATE_CIGAR=true \
      MAX_INSERTIONS_OR_DELETIONS=-1 \
      PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
      PROGRAM_RECORD_ID="bwamem" \
      PROGRAM_GROUP_VERSION="~{bwa_version}" \
      PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
      PROGRAM_GROUP_NAME="bwamem" \
      UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
      ALIGNER_PROPER_PAIR_FLAGS=true \
      UNMAP_CONTAMINANT_READS=true \
      ADD_PG_TAG_TO_READS=false

    java -Xms4000m -jar /usr/gitc/picard.jar \
      MarkDuplicates \
      INPUT=mba.bam \
      OUTPUT=md.bam \
      METRICS_FILE=~{metrics_filename} \
      VALIDATION_STRINGENCY=SILENT \
      ~{"READ_NAME_REGEX=" + read_name_regex} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      ASSUME_SORT_ORDER="queryname" \
      CLEAR_DT="false" \
      ADD_PG_TAG_TO_READS=false

    java -Xms4000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=md.bam \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      MAX_RECORDS_IN_RAM=300000
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: mem + " GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

task GetBwaVersion {
  meta {
    description: "Gets version of BWA"
    volatile: true
  }
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    memory: "1 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }
  output {
    String version = read_string(stdout())
  }
}
