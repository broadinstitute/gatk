task M2 {
  File gatk4_jar
  File intervals
  File ref_fasta 
  File ref_fasta_index 
  File ref_dict 
  File tumor_bam
  File tumor_bam_index
  String tumor_sample_name
  File normal_bam
  File normal_bam_index
  String normal_sample_name
  File? pon
  File? pon_index
  File dbsnp
  File dbsnp_index
  File cosmic
  File cosmic_index

  command <<<
    commandline="java -Xmx4g -jar ${gatk4_jar} Mutect2 \
    -R ${ref_fasta} \
    -I ${tumor_bam} \
    -I ${normal_bam} \
    -tumor ${tumor_sample_name} \
    -normal ${normal_sample_name} \
    --dbsnp ${dbsnp} \
    --cosmic ${cosmic} \
    -L ${intervals} \
    -O ${tumor_sample_name}-vs-${normal_sample_name}.vcf"

    if [ ! -z "$pon" ]
    then
      commandline="$commandline --normal_panel ${pon}" 
    fi

    $commandline
  >>>

  output {
    File output_vcf = "${tumor_sample_name}-vs-${normal_sample_name}.vcf"
  }
}

task GatherVCFs {
  File picard_jar
  Array[File] input_vcfs
  String output_vcf_name

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command <<<
    java -Xmx2g -jar ${picard_jar} \
    MergeVcfs \
    INPUT=${sep=' INPUT=' input_vcfs} \
    OUTPUT=${output_vcf_name}.vcf
  >>>

  output {
    File output_vcf = "${output_vcf_name}.vcf"
    File output_vcf_index = "${output_vcf_name}.vcf.idx"
  }
}

task Filter {
  File gatk4_jar
  File unfiltered_calls
  String output_vcf_name

  command <<<
  	java -Xmx4g -jar ${gatk4_jar} FilterMutectCalls \
  	-V ${unfiltered_calls} \
    -O ${output_vcf_name}.vcf
  >>>

  output {
    File output_vcf = "${output_vcf_name}.vcf"
  }
}

# Warning: this task does not work in the cloud.
task SplitIntervals {
  File picard_jar
  Int scatterCount
  File intervals

  command <<<
    mkdir intervalFileDir

    java -jar ${picard_jar} IntervalListTools \
    I=${intervals} \
    O=intervalFileDir \
    SCATTER_COUNT=${scatterCount}
  >>>

  output {
    Array[File] intervalFiles = glob("intervalFileDir/temp_*/scattered.interval_list")
  }
}

workflow Mutect2 {
  File gatk4_jar
  File intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File tumor_bam
  File tumor_bam_index
  String tumor_sample_name
  File normal_bam
  File normal_bam_index
  String normal_sample_name
  File? pon
  File? pon_index
  Int scatter_count
  File dbsnp
  File dbsnp_index
  File cosmic
  File cosmic_index
  File picard_jar

  call SplitIntervals {
    input:
      picard_jar = picard_jar,
      scatterCount = scatter_count,
      intervals = intervals
  }

  scatter (subintervals in SplitIntervals.intervalFiles ) {
    call M2 {
      input: 
        gatk4_jar = gatk4_jar,
        intervals = subintervals,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        tumor_bam = tumor_bam,
        tumor_bam_index = tumor_bam_index,
        tumor_sample_name = tumor_sample_name,
        normal_bam = normal_bam,
        normal_bam_index = normal_bam_index,
        normal_sample_name = normal_sample_name,
        pon = pon,
        pon_index = pon_index,
        dbsnp = dbsnp,
        dbsnp_index = dbsnp_index,
        cosmic = cosmic,
        cosmic_index = cosmic_index
    }
  }

  call GatherVCFs {
    input:
      picard_jar = picard_jar,
      input_vcfs = M2.output_vcf,
      output_vcf_name = "${tumor_sample_name}-vs-${normal_sample_name}-unfiltered"
  }

call Filter {
    input:
      gatk4_jar = gatk4_jar,
      unfiltered_calls = GatherVCFs.output_vcf,
      output_vcf_name = "${tumor_sample_name}-vs-${normal_sample_name}-filtered"
  }

  output {
        File unfiltered_vcf = GatherVCFs.output_vcf
        File filtered_vcf = Filter.output_vcf
    }
}