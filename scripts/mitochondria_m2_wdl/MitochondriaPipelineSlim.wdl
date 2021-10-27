version 1.0

import "MitoTasks.wdl" as MitoTasks
#import "https://personal.broadinstitute.org/scalvo/WDL/MitoTasks.wdl" as MitoTasks
#import "https://api.firecloud.org/ga4gh/v1/tools/mitochondria:MitoTasks/versions/1/plain-WDL/descriptor" as MitoTasks

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
    String? m2_extra_filtering_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Float? verifyBamID
    Boolean compress_output_vcf = false
    Boolean output_coverage_at_every_base = false

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
 call MitoTasks.SubsetBamToChrM as SubsetBamToChrM {
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

  call MitoTasks.CollectWgsMetrics as CollectWgsMetrics {
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

  call MitoTasks.M2 as CallMt {
    input:
      input_bam = SubsetBamToChrM.output_bam,
      input_bai = SubsetBamToChrM.output_bai,
      mt_interval_list = mt_interval_list,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_args = select_first([m2_extra_args, ""]),
      mem = M2_mem,
      preemptible_tries = preemptible_tries
  }

 # Filtering takes place in two passes
 # The first pass is used to create a list of PASS sites that is
 # passed to Haplochecker to estimate mitochondrial contamination
 # and then a second filter is run using this estimated contamination.
 call MitoTasks.Filter as InitialFilter {
    input:
      raw_vcf = CallMt.raw_vcf,
      raw_vcf_index = CallMt.raw_vcf_index,
      raw_vcf_stats = CallMt.stats,
      base_name = base_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = select_first([m2_extra_filtering_args, ""]),
      vaf_filter_threshold = 0,
      artifact_prone_sites = artifact_prone_sites,
      artifact_prone_sites_index = artifact_prone_sites_index,
      f_score_beta = f_score_beta,
      run_contamination = false,
      preemptible_tries = preemptible_tries
  }

 # This is needed to get input for Haplochecker
 call MitoTasks.SplitMultiAllelicsAndRemoveNonPassSites as SplitMultiAllelicsAndRemoveNonPassSites {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      filtered_vcf = InitialFilter.filtered_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      preemptible_tries = preemptible_tries
  }

  # Call Haplochecker on only PASS sites
  call MitoTasks.GetContamination as GetContamination {
    input:
      input_vcf = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker,
      preemptible_tries = preemptible_tries
  }

  # Now second filter using this estimated contamination
  call MitoTasks.Filter as FilterContamination {
    input:
      raw_vcf = InitialFilter.filtered_vcf,
      raw_vcf_index = InitialFilter.filtered_vcf_index,
      raw_vcf_stats = CallMt.stats,
      run_contamination = true,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
      verifyBamID = verifyBamID,
      base_name = base_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = select_first([m2_extra_filtering_args, ""]),
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      artifact_prone_sites = artifact_prone_sites,
      artifact_prone_sites_index = artifact_prone_sites_index,
      preemptible_tries = preemptible_tries
 }

  call MitoTasks.SplitMultiAllelicSites as SplitMultiAllelicSites {
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

  # This is a temporary task to handle "joint calling" until Mutect2 can produce a GVCF.
  # This proivdes coverage at each base so low coverage sites can be considered ./. rather than 0/0.
  if(output_coverage_at_every_base) {
   call CoverageAtEveryBase {
    input:
      input_bam = SubsetBamToChrM.output_bam,
      input_bai = SubsetBamToChrM.output_bai,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      mt_interval_list = mt_interval_list
   }
  }

  output {
    File out_vcf = FilterContamination.filtered_vcf
    File out_vcf_index = FilterContamination.filtered_vcf_index
    File split_vcf = SplitMultiAllelicSites.split_vcf
    File split_vcf_index = SplitMultiAllelicSites.split_vcf_index
    File coverage_metrics =  CollectWgsMetrics.metrics
    File contamination_metrics = GetContamination.contamination_file
    String major_haplogroup = GetContamination.major_hg
    Float contamination = FilterContamination.contamination
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    File? base_level_coverage_metrics = CoverageAtEveryBase.table
  }
}

task CoverageAtEveryBase {
  input {
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File mt_interval_list

    Int? preemptible_tries
  }
  Int disk_size = ceil(size(input_bam, "GB") + size(ref_fasta, "GB") * 2) + 20

  meta {
    description: "Remove this hack once there's a GVCF solution."
  }

  command <<<
    set -e

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=~{input_bam} \
      R=~{ref_fasta} \
      PER_BASE_COVERAGE=per_base_coverage.tsv \
      O=mtdna.metrics \
      TI=~{mt_interval_list} \
      BI=~{mt_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

  >>>
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "1200 MB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File table = "per_base_coverage.tsv"
  }
}


