version 1.0

workflow MitochondriaPipelineSlim {

  meta {
    description: "Takes in a WGS bam or cram and outputs VCF of SNP/Indel calls on the mitochondria.  Note this pipeline does not perform any realignment and just uses read-pairs that map to chrM."
    allowNestedInputs: true
  }

  input {
    File wgs_aligned_input_bam_or_cram
    File wgs_aligned_input_bam_or_cram_index

    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File mt_interval_list
    File? artifact_prone_sites
    File? artifact_prone_sites_index


    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    String? requester_pays_project
    File? gatk_override
    String? gatk_docker_override
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Float? verifyBamID
    Boolean compress_output_vcf = false

    #Optional runtime arguments
    Int? preemptible_tries
  }

  parameter_meta {
    wgs_aligned_input_bam_or_cram: "Full WGS bam or cram"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
    vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
    f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    mt_interval_list: "Picard style interval list file, with header and single interval representing chrM, eg chrM 1 16569 + ."
  }

  # pull out the base name of either the cram or bam file
  String base_name = basename(basename(wgs_aligned_input_bam_or_cram, ".bam"),".cram")

 # SubsetBamToChrM is called ONLY to reduce computational costs, since passing a smaller input file 
 # to CollectWgsMetrics and CallM2 reduce costs
 call SubsetBamToChrM {
    input:
      input_bam = wgs_aligned_input_bam_or_cram,
      input_bai = wgs_aligned_input_bam_or_cram_index,
      mt_interval_list = mt_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      requester_pays_project = requester_pays_project,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      preemptible_tries = preemptible_tries
  }

  call CollectWgsMetrics {
    input:
      input_bam = SubsetBamToChrM.output_bam,
      input_bai = SubsetBamToChrM.output_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      mt_interval_list = mt_interval_list,
      read_length = max_read_length,
      coverage_cap = 100000,
      preemptible_tries = preemptible_tries
  }

  Int? M2_mem = if CollectWgsMetrics.mean_coverage > 25000 then 14 else 7

  call M2 as CallMt {
    input:
      input_bam = SubsetBamToChrM.output_bam,
      input_bai = SubsetBamToChrM.output_bai,
      mt_interval_list = mt_interval_list,
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_args = select_first([m2_extra_args, ""]),
      mem = M2_mem,
      preemptible_tries = preemptible_tries
  }

  call Filter as InitialFilter {
    input:
      raw_vcf = CallMt.raw_vcf,
      raw_vcf_index = CallMt.raw_vcf_idx,
      raw_vcf_stats = CallMt.stats,
      base_name = base_name,
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = 0,
      f_score_beta = f_score_beta,
      artifact_prone_sites = artifact_prone_sites,
      artifact_prone_sites_index = artifact_prone_sites_index,
      run_contamination = false,
      preemptible_tries = preemptible_tries,
      run_contamination = false,
      hasContamination = false,
      contamination_major = 0,
      contamination_minor = 0
  }

 
  call SplitMultiAllelicsAndRemoveNonPassSites {
    input:
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      filtered_vcf = InitialFilter.filtered_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      preemptible_tries = preemptible_tries
  }

  call GetContamination {
    input:
      input_vcf = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker,
      preemptible_tries = preemptible_tries
  }

  call Filter as FilterContamination {
    input:
      raw_vcf = InitialFilter.filtered_vcf,
      raw_vcf_index = InitialFilter.filtered_vcf_idx,
      raw_vcf_stats = CallMt.stats,
      run_contamination = true,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
      verifyBamID = verifyBamID,
      base_name = base_name,
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      artifact_prone_sites = artifact_prone_sites,
      artifact_prone_sites_index = artifact_prone_sites_index,
      preemptible_tries = preemptible_tries
 }

  # This is a temporary task to handle "joint calling" until Mutect2 can produce a GVCF.
  # This proivdes coverage at each base so low coverage sites can be considered ./. rather than 0/0.
#  call CoverageAtEveryBase {
#    input:
#      input_bam = SubsetBamToChrM.output_bam,
#      input_bai = SubsetBamToChrM.output_bai,
#      mt_interval_list = mt_interval_list,
#      ref_fasta = ref_fasta,
#      ref_fasta_index = ref_fasta_index,
#      preemptible_tries = preemptible_tries
#  }
  
  call SplitMultiAllelicSites {
    input:
      input_vcf = FilterContamination.filtered_vcf,
      base_name = base_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      preemptible_tries = preemptible_tries
  }

  output {
    File out_vcf = FilterContamination.filtered_vcf
    File out_vcf_index = FilterContamination.filtered_vcf_idx
    File split_vcf = SplitMultiAllelicSites.split_vcf
    File split_vcf_index = SplitMultiAllelicSites.split_vcf_index
    File coverage_metrics =  CollectWgsMetrics.metrics
    File contamination_metrics = GetContamination.contamination_file
    #File base_level_coverage_metrics = CoverageAtEveryBase.table
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    String major_haplogroup = GetContamination.major_hg
    Float contamination = FilterContamination.contamination
  }
}

task SubsetBamToChrM {
  input {
    File input_bam
    File input_bai
    File mt_interval_list
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
      -L ~{mt_interval_list} \
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

#task CoverageAtEveryBase {
#  input {
#    File input_bam
#    File input_bai
#    File mt_interval_list
#    File ref_fasta
#    File ref_fasta_index
#    Int? preemptible_tries
#  }
#  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
#  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20
#
#  meta {
#    description: "Remove this hack once there's a GVCF solution."
#  }
#  command <<<
#    set -e
#    
#    java -jar /usr/gitc/picard.jar CollectHsMetrics \
#      I=~{input_bam} \
#      R=~{ref_fasta} \
#      PER_BASE_COVERAGE=per_base_coverage.tsv \
#      O=region.metrics \
#      TI=mt_interval_list \
#      BI=mt_interval_list \
#      COVMAX=20000 \
#      SAMPLE_SIZE=1
#  >>>
#
#  runtime {
#    disks: "local-disk " + disk_size + " HDD"
#    memory: "1200 MB"
#    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
#    preemptible: select_first([preemptible_tries, 5])
#  }
#  output {
#    File table = "per_base_coverage.tsv"
#  }
#}

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

task CollectWgsMetrics {
  input {
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_index
    File mt_interval_list
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
      INTERVALS=~{mt_interval_list} \
      READ_LENGTH=~{read_length_for_optimization} \
      ~{"COVERAGE_CAP=" + coverage_cap} \
      INCLUDE_BQ_HISTOGRAM=true \
      THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.txt

    R --vanilla <<CODE
      df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
      write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
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
  }
}

task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    File mt_interval_list
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
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
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
        -L ~{mt_interval_list} \
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
      File raw_vcf_idx = "~{output_vcf_index}"
      File stats = "~{output_vcf}.stats"
      File output_bamOut = "bamout.bam"
  }
}

task Filter {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File raw_vcf
    File raw_vcf_index
    File raw_vcf_stats
    Boolean compress
    String base_name

    String? m2_extra_filtering_args
    Int max_alt_allele_count
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
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
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
      File filtered_vcf_idx = "~{output_vcf_index}"
      Float contamination = "~{hc_contamination}"
  }
}

task SplitMultiAllelicsAndRemoveNonPassSites {
  input {
    File ref_fasta
    File ref_fai
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
