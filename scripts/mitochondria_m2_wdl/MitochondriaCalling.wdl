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
  Float? autosomal_coverage
  Float contamination
  #max_read_length is only used for an optimization. If it's too small CollectWgsMetrics might fail, but the results are not affected by this number.
  Int? max_read_length

  String bn = basename(input_bam, ".bam")

  #Optional runtime arguments
  Int? preemptible_tries

  call CollectWgsMetrics {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      read_length = max_read_length,
      coverage_cap = 100000,
      preemptible_tries = preemptible_tries
  }

  Int? M2_mem = if CollectWgsMetrics.mean_coverage > 25000 then 14 else 7

  call M2AndFilter {
    input:
      ref_fasta = ref_fasta,
      ref_fai = ref_fasta_index,
      ref_dict = ref_dict,
      input_bam = input_bam,
      input_bai = input_bam_index,
      compress = false,
      # The minimum pruning value for the assembly graph in M2 defaults to 1. This only makes sense at low depths so we
      # need this number to scale with mean coverage. Ideally this knob should be removed once dynamic pruning has been
      # added to M2.
      prune = (CollectWgsMetrics.mean_coverage / 500) + 1,
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
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    Int mean_coverage = CollectWgsMetrics.mean_coverage
    File vcf = M2AndFilter.filtered_vcf
    File vcf_index = M2AndFilter.filtered_vcf_idx
  }
}

task CollectWgsMetrics {
  File input_bam
  File input_bam_index
  File ref_fasta
  File ref_fasta_index
  Int? read_length
  Int read_length_for_optimization = select_first([read_length, 151])
  Int? coverage_cap

  Int? preemptible_tries
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
      INPUT=${input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=${ref_fasta} \
      OUTPUT=metrics.txt \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=${read_length_for_optimization} \
      ${"COVERAGE_CAP=" + coverage_cap} \
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
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
  }
  output {
    File metrics = "metrics.txt"
    File theoretical_sensitivity = "theoretical_sensitivity.txt"
    Int mean_coverage = read_int("mean_coverage.txt")
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
  Int prune
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
        --min-pruning ${prune} \
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
      docker: "us.gcr.io/broad-gatk/gatk:4.0.11.0"
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
