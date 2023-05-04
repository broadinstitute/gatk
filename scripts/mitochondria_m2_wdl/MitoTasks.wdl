version 1.0

task SubsetBamToChrM {
  input {
    File input_bam
    File input_bai
    File? mt_interval_list
    String? contig_name
    String basename = basename(basename(input_bam, ".cram"), ".bam")
    String? requester_pays_project
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

    # runtime
    File? gatk_override
    String? gatk_docker_override
    Int? preemptible_tries
  }
  Float ref_size = if defined(ref_fasta) then size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB") else 0
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  meta {
    description: "Subsets a whole genome bam to just Mitochondria reads"
  }
  parameter_meta {
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
    input_bam: {
      localization_optional: true
    }
    input_bai: {
      localization_optional: true
    }
  }
  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk PrintReads \
      ~{"-R " + ref_fasta} \
      ~{"-L " + mt_interval_list} \
      ~{"-L " + contig_name} \
      --read-filter MateOnSameContigOrNoMappedMateReadFilter \
      --read-filter MateUnmappedAndUnmappedReadFilter \
      ~{"--gcs-project-for-requester-pays " + requester_pays_project} \
      -I ~{input_bam} \
      --read-index ~{input_bai} \
      -O ~{basename}.bam
  >>>
  runtime {
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File output_bam = "~{basename}.bam"
    File output_bai = "~{basename}.bai"
  }
}

task CollectWgsMetrics {
  input {
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_index
    File? mt_interval_list
    Int? read_length
    Int? coverage_cap

    Int? preemptible_tries
  }

  Int read_length_for_optimization = select_first([read_length, 151])
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  meta {
    description: "Collect coverage metrics"
  }
  parameter_meta {
    read_length: "Read length used for optimization only. If this is too small CollectWgsMetrics might fail. Default is 151."
  }

  command <<<
    set -e

    java -Xms2000m -jar /usr/gitc/picard.jar \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=metrics.txt \
      USE_FAST_ALGORITHM=true \
      ~{"INTERVALS=" + mt_interval_list} \
      READ_LENGTH=~{read_length_for_optimization} \
      ~{"COVERAGE_CAP=" + coverage_cap} \
      INCLUDE_BQ_HISTOGRAM=true \
      THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.txt

    R --vanilla <<CODE
      df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
      write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
      write.table(df[,"MEDIAN_COVERAGE"], "median_coverage.txt", quote=F, col.names=F, row.names=F)
    CODE
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
  }
  output {
    File metrics = "metrics.txt"
    File theoretical_sensitivity = "theoretical_sensitivity.txt"
    Int mean_coverage = read_int("mean_coverage.txt")
    Float median_coverage = read_float("median_coverage.txt")
  }
}

task M2 {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_bam
    File input_bai
    File? mt_interval_list
    Int? max_reads_per_alignment_start
    String? m2_extra_args
    Boolean? make_bamout
    Boolean compress
    File? gatk_override
    # runtime
    String? gatk_docker_override
    Int? mem
    Int? preemptible_tries
  }

  Int max_reads_per_alignment_start_arg = select_first([max_reads_per_alignment_start, 75])
  String output_vcf = "raw" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  meta {
    description: "Mutect2 for calling Snps and Indels"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        ~{"-L " + mt_interval_list} \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O ~{output_vcf} \
        ~{true='--bam-output bamout.bam' false='' make_bamout} \
        ~{m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start ~{max_reads_per_alignment_start_arg} \
        --max-mnp-distance 0
  >>>
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: machine_mem + " MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File raw_vcf = "~{output_vcf}"
      File raw_vcf_index = "~{output_vcf_index}"
      File stats = "~{output_vcf}.stats"
      File output_bamOut = "bamout.bam"
  }
}

task SplitMultiAllelicSites {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_vcf
    String base_name
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
  }

  String output_vcf = base_name + ".final.split.vcf"
  String output_vcf_index = output_vcf + ".idx"

  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    gatk LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{input_vcf} \
      -O ~{output_vcf} \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac
  }
  output {
    File split_vcf = "~{output_vcf}"
    File split_vcf_index = "~{output_vcf}"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  } 
}

task SplitMultiAllelicsAndRemoveNonPassSites {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File filtered_vcf
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
  }

  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    gatk LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

      gatk SelectVariants \
        -V split.vcf \
        -O splitAndPassOnly.vcf \
        --exclude-filtered
  
  }
  output {
    File vcf_for_haplochecker = "splitAndPassOnly.vcf"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  } 
}

task GetContamination {
  input {
    File input_vcf
    # runtime
    Int? preemptible_tries
  }

  Int disk_size = ceil(size(input_vcf, "GB")) + 20

  meta {
    description: "Uses new Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_vcf: "Filtered and split multi-allelic sites VCF for mitochondria"
  }
  command <<<
  set -e
  PARENT_DIR="$(dirname "~{input_vcf}")"
  java -jar /usr/mtdnaserver/haplocheckCLI.jar "${PARENT_DIR}"

  sed 's/\"//g' output > output-noquotes

  grep "SampleID" output-noquotes > headers
  FORMAT_ERROR="Bad contamination file format"
  if [ `awk '{print $2}' headers` != "Contamination" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $6}' headers` != "HgMajor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $8}' headers` != "HgMinor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $14}' headers` != "MeanHetLevelMajor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $15}' headers` != "MeanHetLevelMinor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi

  grep -v "SampleID" output-noquotes > output-data
  awk -F "\t" '{print $2}' output-data > contamination.txt
  awk -F "\t" '{print $6}' output-data > major_hg.txt
  awk -F "\t" '{print $8}' output-data > minor_hg.txt
  awk -F "\t" '{print $14}' output-data > mean_het_major.txt
  awk -F "\t" '{print $15}' output-data > mean_het_minor.txt
  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-dsde-methods/haplochecker:haplochecker-0124"
  }
  output {
    File contamination_file = "output-noquotes"
    String hasContamination = read_string("contamination.txt") 
    String major_hg = read_string("major_hg.txt")
    String minor_hg = read_string("minor_hg.txt")
    Float major_level = read_float("mean_het_major.txt")
    Float minor_level = read_float("mean_het_minor.txt")
  }
}

task Filter {
  input {
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File raw_vcf
    File raw_vcf_index
    File raw_vcf_stats
    Boolean compress
    String base_name

    String? m2_extra_filtering_args
    Int max_alt_allele_count = 4
    Float? vaf_filter_threshold
    Float? f_score_beta

    Boolean run_contamination
    String? hasContamination
    Float? contamination_major
    Float? contamination_minor
    Float? verifyBamID 

    File? artifact_prone_sites
    File? artifact_prone_sites_index

    File? gatk_override
    String? gatk_docker_override

  # runtime
    Int? preemptible_tries
  }

  String output_vcf = base_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(raw_vcf, "GB") + ref_size) + 20
  Float hc_contamination = if run_contamination && hasContamination == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
  Float max_contamination = if defined(verifyBamID) && verifyBamID > hc_contamination then verifyBamID else hc_contamination

  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
      f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx2500m" FilterMutectCalls -V ~{raw_vcf} \
        -R ~{ref_fasta} \
        -O filtered.vcf \
        --stats ~{raw_vcf_stats} \
        ~{m2_extra_filtering_args} \
        --max-alt-allele-count ~{max_alt_allele_count} \
        --mitochondria-mode \
        ~{"--min-allele-fraction " + vaf_filter_threshold} \
        ~{"--f-score-beta " + f_score_beta} \
        ~{"--contamination-estimate " + max_contamination}

      gatk VariantFiltration -V filtered.vcf \
        -O ~{output_vcf} \
        --apply-allele-specific-filters \
        ~{"--mask " + artifact_prone_sites} \
        --mask-name "artifact_prone_site"

  >>>
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: "4 MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File filtered_vcf = "~{output_vcf}"
      File filtered_vcf_index = "~{output_vcf_index}"
      Float contamination = "~{hc_contamination}"
  }
}

