workflow MitochondriaCalling {

  meta {
    description: "SNP/Indel calling for mitochondria using Mutect2. Takes aligned mitochondria bam and outputs filtered VCF plus coverage metrics."
  }
  parameter_meta {
    input_bam: "Bam of reads realigned to the mitochondria"
    blacklisted_sites: "List of sites that will always be filtered"
    vcf: "VCF of SNPs and INDELs called and filtered by Mutect2"
  }
  File input_bam
  File input_bam_index
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  Float lod_cutoff
  File? gatk_override
  String? m2_extra_args
  File blacklisted_sites
  File blacklisted_sites_index
  Int mean_coverage
  Float? autosomal_coverage
  Float contamination
  #max_read_length is only used for an optimization. If it's too small CollectWgsMetrics might fail, but the results are not affected by this number.
  Int? max_read_length

  String bn = basename(input_bam, ".bam")

  #Optional runtime arguments
  Int? preemptible_tries

  Int? M2_mem = if mean_coverage > 25000 then 14 else 7

  call M2AndFilter {
    input:
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      input_bam = input_bam,
      input_bai = input_bam_index,
      compress = false,
      max_alt_allele_count = 4,
      lod = lod_cutoff,
      contamination = contamination,
      gatk_override = gatk_override,
      m2_extra_args = m2_extra_args,
      mem = M2_mem,
      autosomal_coverage = autosomal_coverage,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      preemptible_tries = preemptible_tries
  }

  output {
    File vcf = M2AndFilter.filtered_vcf
    File vcf_index = M2AndFilter.filtered_vcf_idx
  }
}

task M2AndFilter {
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
  Float? autosomal_coverage

  String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"

  String? m2_extra_filtering_args
  Int max_alt_allele_count
  Float lod
  Float contamination

  File blacklisted_sites
  File blacklisted_sites_index

  File? gatk_override

  # runtime
  Int? mem
  Int? preemptible_tries
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
    prune: "Pruning factor should be approximately (mean coverage / 500)"
    autosomal_coverage: "Median coverage of the autosomes for annotating potential polymorphic NuMT variants"
  }
  command <<<
      set -e

      export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      gatk --java-options "-Xmx${command_mem}m" Mutect2 \
        -R ${ref_fasta} \
        -I ${input_bam} \
        ${"--genotyping-mode GENOTYPE_GIVEN_ALLELES --alleles " + gga_vcf} \
        -O raw.vcf \
        ${true='--bam-output bamout.bam' false='' make_bamout} \
        ${m2_extra_args} \
        ${"--median-autosomal-coverage " + autosomal_coverage} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0

      gatk --java-options "-Xmx2500m" FilterMutectCalls -V raw.vcf \
        -O filtered.vcf \
        ${m2_extra_filtering_args} \
        --max-alt-allele-count ${max_alt_allele_count} \
        --tumor-lod ${lod} \
        --mitochondria-mode \
        --contamination-estimate ${contamination}

      gatk VariantFiltration -V filtered.vcf \
        -O ${output_vcf} \
        --mask ${blacklisted_sites} \
        --mask-name "blacklisted_site"
  >>>
  runtime {
      docker: "us.gcr.io/broad-gatk/gatk:4.1.0.0"
      memory: machine_mem + " MB"
      disks: "local-disk " + disk_size + " HDD"
      preemptible: select_first([preemptible_tries, 5])
      cpu: 2
  }
  output {
      File filtered_vcf = "${output_vcf}"
      File filtered_vcf_idx = "${output_vcf_index}"
      File output_bamOut = "bamout.bam"
  }
}
