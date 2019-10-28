workflow ReadsPipelineSparkWorkflow {
  File ref_fasta
  File ref_fasta_index
  File? ref_fasta_gzi
  File ref_dict
  File input_bam
  String output_vcf_basename
  File known_sites
  File known_sites_index

  String gatk_docker
  String gatk
  Int? num_cpu_override
  Int num_cpu = select_first([num_cpu_override, 16])
  # max 3000GB for local SSDs, see https://cloud.google.com/compute/docs/disks/#localssds
  Int? disk_size_override
  Int disk_size = select_first([disk_size_override, 3000])

  call ReadsPipelineSpark {
    input:
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_fasta_gzi = ref_fasta_gzi,
      ref_dict = ref_dict,
      input_bam = input_bam,
      output_vcf_basename = output_vcf_basename,
      known_sites = known_sites,
      known_sites_index = known_sites_index,
      gatk_docker = gatk_docker,
      gatk = gatk,
      num_cpu = num_cpu,
      disk_size = disk_size
  }

  output {
    File output_vcf = ReadsPipelineSpark.output_vcf
  }
}

task ReadsPipelineSpark {
  File ref_fasta
  File ref_fasta_index
  File? ref_fasta_gzi
  File ref_dict
  File input_bam
  String output_vcf_basename
  File known_sites
  File known_sites_index

  String gatk_docker
  String gatk
  Int num_cpu
  Int disk_size

  command {
    ${gatk} \
      ReadsPipelineSpark \
      -R ${ref_fasta} \
      -I ${input_bam} \
      -O ${output_vcf_basename}.vcf \
      --known-sites ${known_sites} \
      -pairHMM AVX_LOGLESS_CACHING \
      --conf 'spark.local.dir=./tmp'
  }

  runtime {
    docker: gatk_docker
    cpu: num_cpu
    disks: "local-disk " + disk_size + " LOCAL"
  }

  output {
    File output_vcf = "${output_vcf_basename}.vcf"
  }
}