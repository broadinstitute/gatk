version 1.0

import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates
#import "https://api.firecloud.org/ga4gh/v1/tools/mitochondria:AlignmentPipeline/versions/1/plain-WDL/descriptor" as AlignAndMarkDuplicates

import "MitoTasks.wdl" as MitoTasks
#import "https://api.firecloud.org/ga4gh/v1/tools/mitochondria:MitoTasks/versions/1/plain-WDL/descriptor" as MitoTasks
#import "https://personal.broadinstitute.org/scalvo/WDL/MitoTasks.wdl" as MitoTasks

workflow AlignAndCall {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File unmapped_bam
    Float? autosomal_coverage
    String base_name

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa
    File blacklisted_sites
    File blacklisted_sites_index

    #Shifted reference is used for calling the control region (edge of mitochondria reference).
    #This solves the problem that BWA doesn't support alignment to circular contigs.
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa

    File shift_back_chain

    File? gatk_override
    String? gatk_docker_override
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf
    Float? verifyBamID
    Int? max_low_het_sites

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    #Optional runtime arguments
    Int? preemptible_tries
  }

  parameter_meta {
    unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      preemptible_tries = preemptible_tries
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt {
    input:
      input_bam = unmapped_bam,
      mt_dict = mt_shifted_dict,
      mt_fasta = mt_shifted_fasta,
      mt_fasta_index = mt_shifted_fasta_index,
      mt_amb = mt_shifted_amb,
      mt_ann = mt_shifted_ann,
      mt_bwt = mt_shifted_bwt,
      mt_pac = mt_shifted_pac,
      mt_sa = mt_shifted_sa,
      preemptible_tries = preemptible_tries
  }

  call MitoTasks.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bai = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      read_length = max_read_length,
      coverage_cap = 100000,
      preemptible_tries = preemptible_tries
  }

  Int? M2_mem = if CollectWgsMetrics.mean_coverage > 25000 then 14 else 7

  call MitoTasks.M2 as CallMt {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bai = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:576-16024 ",
      mem = M2_mem,
      preemptible_tries = preemptible_tries
  }

  call MitoTasks.M2 as CallShiftedMt {
    input:
      input_bam = AlignToShiftedMt.mt_aligned_bam,
      input_bai = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta = mt_shifted_fasta,
      ref_fasta_index = mt_shifted_fasta_index,
      ref_dict = mt_shifted_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:8025-9144 ",
      mem = M2_mem,
      preemptible_tries = preemptible_tries
  }

  call LiftoverAndCombineVcfs {
    input:
      shifted_vcf = CallShiftedMt.raw_vcf,
      vcf = CallMt.raw_vcf,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shift_back_chain = shift_back_chain,
      preemptible_tries = preemptible_tries
  }

  call MergeStats {
    input:
      shifted_stats = CallShiftedMt.stats,
      non_shifted_stats = CallMt.stats,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      preemptible_tries = preemptible_tries
  }

  call MitoTasks.Filter as InitialFilter {
    input:
      raw_vcf = LiftoverAndCombineVcfs.merged_vcf,
      raw_vcf_index = LiftoverAndCombineVcfs.merged_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      vaf_filter_threshold = 0,
      artifact_prone_sites = blacklisted_sites,
      artifact_prone_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      run_contamination = false,
      preemptible_tries = preemptible_tries
  }

 
  call MitoTasks.SplitMultiAllelicsAndRemoveNonPassSites as SplitMultiAllelicsAndRemoveNonPassSites {
    input:
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      filtered_vcf = InitialFilter.filtered_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override
  }

  call MitoTasks.GetContamination as GetContamination {
    input:
      input_vcf = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker,
      preemptible_tries = preemptible_tries
  }

  call MitoTasks.Filter as FilterContamination {
    input:
      raw_vcf = InitialFilter.filtered_vcf,
      raw_vcf_index = InitialFilter.filtered_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      run_contamination = true,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
      verifyBamID = verifyBamID,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      artifact_prone_sites = blacklisted_sites,
      artifact_prone_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      preemptible_tries = preemptible_tries
 }

  if ( defined(autosomal_coverage) ) {
    call FilterNuMTs {
      input:
        filtered_vcf = FilterContamination.filtered_vcf,
        ref_fasta = mt_fasta,
        ref_fai = mt_fasta_index,
        ref_dict = mt_dict,
        autosomal_coverage = autosomal_coverage,
        gatk_override = gatk_override,
        gatk_docker_override = gatk_docker_override,
        compress = compress_output_vcf,
        preemptible_tries = preemptible_tries
    }
  }

  File low_het_vcf = select_first([FilterNuMTs.numt_filtered_vcf, FilterContamination.filtered_vcf])

  call FilterLowHetSites {
    input:
      filtered_vcf = low_het_vcf,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      max_low_het_sites = max_low_het_sites,
      gatk_override = gatk_override,
      gatk_docker_override = gatk_docker_override,
      compress = compress_output_vcf,
      base_name = base_name,
      preemptible_tries = preemptible_tries
  }


  output {
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File out_vcf = FilterLowHetSites.final_filtered_vcf
    File out_vcf_index = FilterLowHetSites.final_filtered_vcf_idx
    File input_vcf_for_haplochecker = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
    File duplicate_metrics = AlignToMt.duplicate_metrics
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    File contamination_metrics = GetContamination.contamination_file
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    Float median_coverage = CollectWgsMetrics.median_coverage
    String major_haplogroup = GetContamination.major_hg
    Float contamination = FilterContamination.contamination
  }
}


task LiftoverAndCombineVcfs {
  input {
    File shifted_vcf
    File vcf
    String basename = basename(shifted_vcf, ".vcf")

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File shift_back_chain

    # runtime
    Int? preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(shifted_vcf, "GB") + ref_size) + 20

  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    vcf: "VCF of the rest of chrM on original reference"
    ref_fasta: "Original (not shifted) chrM reference"
    shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    java -jar /usr/gitc/picard.jar LiftoverVcf \
      I=~{shifted_vcf} \
      O=~{basename}.shifted_back.vcf \
      R=~{ref_fasta} \
      CHAIN=~{shift_back_chain} \
      REJECT=~{basename}.rejected.vcf

    java -jar /usr/gitc/picard.jar MergeVcfs \
      I=~{basename}.shifted_back.vcf \
      I=~{vcf} \
      O=~{basename}.merged.vcf
    >>>
    runtime {
      disks: "local-disk " + disk_size + " HDD"
      memory: "1200 MB"
      docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
      preemptible: select_first([preemptible_tries, 5])
    }
    output{
        # rejected_vcf should always be empty
        File rejected_vcf = "~{basename}.rejected.vcf"
        File merged_vcf = "~{basename}.merged.vcf"
        File merged_vcf_index = "~{basename}.merged.vcf.idx"
    }
}

task MergeStats {
  input {
    File shifted_stats
    File non_shifted_stats
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
  }

  command{
    set -e

    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk MergeMutectStats --stats ~{shifted_stats} --stats ~{non_shifted_stats} -O raw.combined.stats
  }
  output {
    File stats = "raw.combined.stats"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  }
}

task FilterNuMTs {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    Float? autosomal_coverage
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
    Boolean compress
  }
  
  String basename = basename(filtered_vcf, ".vcf")
  String output_vcf = basename + ".numt" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  
  parameter_meta {
    autosomal_coverage: "Median coverage of the autosomes for filtering potential polymorphic NuMT variants"
  }

  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    gatk NuMTFilterTool \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O ~{output_vcf} \
      --autosomal-coverage ~{autosomal_coverage}
  
  }
  output {
    File numt_filtered_vcf = "~{output_vcf}"
    File numt_filtered_vcf_idx = "~{output_vcf_index}"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  } 
}

task FilterLowHetSites {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    String base_name
    Int? max_low_het_sites
    Int? preemptible_tries
    File? gatk_override
    String? gatk_docker_override
    Boolean compress
  }

  String output_vcf = base_name + ".final" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Int max_sites = select_first([max_low_het_sites, 3])
  
  command {
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
    gatk MTLowHeteroplasmyFilterTool \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O ~{output_vcf} \
      --max-allowed-low-hets ~{max_sites}
  
  }
  output {
    File final_filtered_vcf = "~{output_vcf}"
    File final_filtered_vcf_idx = "~{output_vcf_index}"
  }
  runtime {
      docker: select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.1.7.0"])
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  } 
}
