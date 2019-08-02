version 1.0

import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates

workflow AlignAndCall {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File unmapped_bam
    Float? autosomal_coverage

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

    File? gatk_override
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta
    Boolean compress_output_vcf

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

  call CollectWgsMetrics {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      read_length = max_read_length,
      coverage_cap = 100000,
      preemptible_tries = preemptible_tries
  }

  call GetContamination {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      preemptible_tries = preemptible_tries
  }

  Int? M2_mem = if CollectWgsMetrics.mean_coverage > 25000 then 14 else 7

  call M2 as CallMt {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bai = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:576-16024 ",
      mem = M2_mem,
      preemptible_tries = preemptible_tries
  }

  call M2 as CallShiftedMt {
    input:
      input_bam = AlignToShiftedMt.mt_aligned_bam,
      input_bai = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta = mt_shifted_fasta,
      ref_fai = mt_shifted_fasta_index,
      ref_dict = mt_shifted_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
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
      preemptible_tries = preemptible_tries
  }

  call Filter {
    input:
      raw_vcf = LiftoverAndCombineVcfs.final_vcf,
      raw_vcf_index = LiftoverAndCombineVcfs.final_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      gatk_override = gatk_override,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      contamination = GetContamination.minor_level,
      autosomal_coverage = autosomal_coverage,
      vaf_filter_threshold = vaf_filter_threshold,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      preemptible_tries = preemptible_tries
  }

  output {
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File out_vcf = Filter.filtered_vcf
    File out_vcf_index = Filter.filtered_vcf_idx
    File duplicate_metrics = AlignToMt.duplicate_metrics
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    File contamination_metrics = GetContamination.contamination_file
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    String major_haplogroup = GetContamination.major_hg
    Float contamination = GetContamination.minor_level
  }
}

task GetContamination {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    Int qual = 20
    Int map_qual = 30
    Float vaf = 0.01

    # runtime
    Int? preemptible_tries
  }

  String basename = basename(input_bam, ".bam")
  Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB")
  Int disk_size = ceil(size(input_bam, "GB") + ref_size) + 20

  meta {
    description: "Uses Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_bam: "Bam aligned to chrM"
    ref_fasta: "chrM reference"
  }
  command {
  set -e

  java -jar /usr/mtdnaserver/mitolib.jar haplochecker \
    --in ~{input_bam} \
    --ref ~{ref_fasta} \
    --out haplochecker_out \
    --QUAL ~{qual} \
    --MAPQ ~{map_qual} \
    --VAF ~{vaf}

python3 <<CODE

import csv

with open("haplochecker_out/~{basename}.contamination.txt") as output:
    reader = csv.DictReader(output, delimiter='\t')
    for row in reader:
        print(row["MajorHG"], file=open("major_hg.txt", 'w'))
        print(row["MajorLevel"], file=open("major_level.txt", 'w'))
        print(row["MinorHG"], file=open("minor_hg.txt", 'w'))
        print(row["MinorLevel"], file=open("minor_level.txt", 'w'))
CODE
  }
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "gatkworkflows/mtdnaserver:1.2"
  }
  output {
    File contamination_file = "haplochecker_out/~{basename}.contamination.txt"
    String major_hg = read_string("major_hg.txt")
    Float major_level = read_float("major_level.txt")
    String minor_hg = read_string("minor_hg.txt")
    Float minor_level = read_float("minor_level.txt")
  }
}

task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
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
      O=~{basename}.final.vcf
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
        File final_vcf = "~{basename}.final.vcf"
        File final_vcf_index = "~{basename}.final.vcf.idx"
    }
}

task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    String? m2_extra_args
    Boolean? make_bamout
    Boolean compress
    File? gga_vcf
    File? gga_vcf_idx
    File? gatk_override
    # runtime
    Int? mem
    Int? preemptible_tries
  }

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
    gga_vcf: "VCF for genotype given alleles mode"
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx~{command_mem}m" Mutect2 \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        ~{"--genotyping-mode GENOTYPE_GIVEN_ALLELES --alleles " + gga_vcf} \
        -O ~{output_vcf} \
        ~{true='--bam-output bamout.bam' false='' make_bamout} \
        ~{m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0
  >>>
  runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
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
    Float? vaf_cutoff

    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float contamination
    Float? autosomal_coverage
    Float? vaf_filter_threshold
    Float? f_score_beta

    File blacklisted_sites
    File blacklisted_sites_index

    File? gatk_override

  # runtime
    Int? preemptible_tries
  }

  String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Float ref_size = size(ref_fasta, "GB") + size(ref_fai, "GB")
  Int disk_size = ceil(size(raw_vcf, "GB") + ref_size) + 20

  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      autosomal_coverage: "Median coverage of the autosomes for filtering potential polymorphic NuMT variants"
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
        ~{"--autosomal-coverage " + autosomal_coverage} \
        ~{"--min-allele-fraction " + vaf_filter_threshold} \
        ~{"--f-score-beta " + f_score_beta} \
        --contamination-estimate ~{contamination}

      gatk VariantFiltration -V filtered.vcf \
        -O ~{output_vcf} \
        --mask ~{blacklisted_sites} \
        --mask-name "blacklisted_site"

  >>>
  runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
      memory: "4 MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File filtered_vcf = "~{output_vcf}"
      File filtered_vcf_idx = "~{output_vcf_index}"
  }
}

task MergeStats {
  input {
    File shifted_stats
    File non_shifted_stats
    Int? preemptible_tries
    File? gatk_override
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
      docker: "us.gcr.io/broad-gatk/gatk:4.1.1.0"
      memory: "3 MB"
      disks: "local-disk 20 HDD"
      preemptible: select_first([preemptible_tries, 5])
  }
}
