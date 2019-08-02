version 1.0

import "AlignAndCall.wdl" as AlignAndCall

workflow MitochondriaPipeline {

  meta {
    description: "Takes in fully aligned hg38 bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File wgs_aligned_input_bam_or_cram
    String contig_name = "chrM"
    Float? autosomal_coverage

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length

    # Full reference is only requred if starting with a CRAM (BAM doesn't need these files)
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

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
    File blacklisted_sites_shifted
    File blacklisted_sites_shifted_index

    File shift_back_chain

    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list

    File? gatk_override
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf = false

    #Optional runtime arguments
    Int? preemptible_tries
  }

  parameter_meta {
    wgs_aligned_input_bam_or_cram: "Full WGS hg38 bam or cram"
    autosomal_coverage: "Median coverage of full input bam"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
    vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
    f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    contig_name: "Name of mitochondria contig in reference that wgs_aligned_input_bam_or_cram is aligned to"
  }

  call SubsetBamToChrM {
    input:
      input_bam = wgs_aligned_input_bam_or_cram,
      contig_name = contig_name,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      gatk_override = gatk_override,
      preemptible_tries = preemptible_tries
  }

  call RevertSam {
    input:
      input_bam = SubsetBamToChrM.output_bam,
      preemptible_tries = preemptible_tries
  }

  call AlignAndCall.AlignAndCall as AlignAndCall {
    input:
      unmapped_bam = RevertSam.unmapped_bam,
      autosomal_coverage = autosomal_coverage,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      mt_shifted_dict = mt_shifted_dict,
      mt_shifted_fasta = mt_shifted_fasta,
      mt_shifted_fasta_index = mt_shifted_fasta_index,
      mt_shifted_amb = mt_shifted_amb,
      mt_shifted_ann = mt_shifted_ann,
      mt_shifted_bwt = mt_shifted_bwt,
      mt_shifted_pac = mt_shifted_pac,
      mt_shifted_sa = mt_shifted_sa,
      blacklisted_sites_shifted = blacklisted_sites_shifted,
      blacklisted_sites_shifted_index = blacklisted_sites_shifted_index,
      shift_back_chain = shift_back_chain,
      gatk_override = gatk_override,
      m2_extra_args = m2_extra_args,
      m2_filter_extra_args = m2_filter_extra_args,
      vaf_filter_threshold = vaf_filter_threshold,
      f_score_beta = f_score_beta,
      compress_output_vcf = compress_output_vcf,
      max_read_length = max_read_length,
      preemptible_tries = preemptible_tries
  }

  # This is a temporary task to handle "joint calling" until Mutect2 can produce a GVCF.
  # This proivdes coverage at each base so low coverage sites can be considered ./. rather than 0/0.
  call CoverageAtEveryBase {
    input:
      input_bam_regular_ref = AlignAndCall.mt_aligned_bam,
      input_bam_regular_ref_index = AlignAndCall.mt_aligned_bai,
      input_bam_shifted_ref = AlignAndCall.mt_aligned_shifted_bam,
      input_bam_shifted_ref_index = AlignAndCall.mt_aligned_shifted_bai,
      shift_back_chain = shift_back_chain,
      control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
      non_control_region_interval_list = non_control_region_interval_list,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shifted_ref_fasta = mt_shifted_fasta,
      shifted_ref_fasta_index = mt_shifted_fasta_index,
      shifted_ref_dict = mt_shifted_dict
  }

  output {
    File subset_bam = SubsetBamToChrM.output_bam
    File subset_bai = SubsetBamToChrM.output_bai
    File mt_aligned_bam = AlignAndCall.mt_aligned_bam
    File mt_aligned_bai = AlignAndCall.mt_aligned_bai
    File out_vcf = AlignAndCall.out_vcf
    File out_vcf_index = AlignAndCall.out_vcf_index
    File duplicate_metrics = AlignAndCall.duplicate_metrics
    File coverage_metrics = AlignAndCall.coverage_metrics
    File theoretical_sensitivity_metrics = AlignAndCall.theoretical_sensitivity_metrics
    File contamination_metrics = AlignAndCall.contamination_metrics
    File base_level_coverage_metrics = CoverageAtEveryBase.table
    Int mean_coverage = AlignAndCall.mean_coverage
    String major_haplogroup = AlignAndCall.major_haplogroup
    Float contamination = AlignAndCall.contamination
  }
}

task SubsetBamToChrM {
  input {
    File input_bam
    String contig_name
    String basename = basename(basename(input_bam, ".cram"), ".bam")
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

    File? gatk_override

    # runtime
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
  }
  command <<<
    set -e
    export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

    gatk PrintReads \
      ~{"-R " + ref_fasta} \
      -L ~{contig_name} \
      --read-filter MateOnSameContigOrNoMappedMateReadFilter \
      --read-filter MateUnmappedAndUnmappedReadFilter \
      -I ~{input_bam} \
      -O ~{basename}.bam
  >>>
  runtime {
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File output_bam = "~{basename}.bam"
    File output_bai = "~{basename}.bai"
  }
}

task RevertSam {
  input {
    File input_bam
    String basename = basename(input_bam, ".bam")

    # runtime
    Int? preemptible_tries
  }
  Int disk_size = ceil(size(input_bam, "GB") * 2.5) + 20

  meta {
    description: "Removes alignment information while retaining recalibrated base qualities and original alignment tags"
  }
  parameter_meta {
    input_bam: "aligned bam"
  }
  command {
    java -Xmx1000m -jar /usr/gitc/picard.jar \
    RevertSam \
    INPUT=~{input_bam} \
    OUTPUT_BY_READGROUP=false \
    OUTPUT=~{basename}.bam \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=false
  }
  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1552931386"
    preemptible: select_first([preemptible_tries, 5])
  }
  output {
    File unmapped_bam = "~{basename}.bam"
  }
}

task CoverageAtEveryBase {
  input {
    File input_bam_regular_ref
    File input_bam_regular_ref_index
    File input_bam_shifted_ref
    File input_bam_shifted_ref_index
    File shift_back_chain
    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File shifted_ref_fasta
    File shifted_ref_fasta_index
    File shifted_ref_dict

    Int? preemptible_tries
  }
  Int disk_size = ceil(size(input_bam_regular_ref, "GB") + size(input_bam_shifted_ref, "GB") + size(ref_fasta, "GB") * 2) + 20

  meta {
    description: "Remove this hack once there's a GVCF solution."
  }
  command <<<
    set -e

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=~{input_bam_regular_ref} \
      R=~{ref_fasta} \
      PER_BASE_COVERAGE=non_control_region.tsv \
      O=non_control_region.metrics \
      TI=~{non_control_region_interval_list} \
      BI=~{non_control_region_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    java -jar /usr/gitc/picard.jar CollectHsMetrics \
      I=~{input_bam_shifted_ref} \
      R=~{shifted_ref_fasta} \
      PER_BASE_COVERAGE=control_region_shifted.tsv \
      O=control_region_shifted.metrics \
      TI=~{control_region_shifted_reference_interval_list} \
      BI=~{control_region_shifted_reference_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    R --vanilla <<CODE
      shift_back = function(x) {
        if (x < 8570) {
          return(x + 8000)
        } else {
          return (x - 8569)
        }
      }

      control_region_shifted = read.table("control_region_shifted.tsv", header=T)
      shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
      control_region_shifted[,"pos"] = shifted_back

      beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
      end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

      non_control_region = read.table("non_control_region.tsv", header=T)
      combined_table = rbind(beginning, non_control_region, end)
      write.table(combined_table, "per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")

    CODE
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
