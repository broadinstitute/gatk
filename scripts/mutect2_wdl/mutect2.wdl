#  Run Mutect 2 on a single tumor-normal pair
#
#  Description of inputs
#  gatk4_jar: java jar file containing gatk 4 (protected)
#  intervals: genomic intervals
#  ref_fasta, ref_fasta_index, ref_dict: reference genome, index, and dictionary
#  tumor_bam, tumor_bam_index, and tumor_sample_name: self-explanatory
#  normal_bam, normal_bam_index, and normal_sample_name: self-explanatory
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  dbsnp, dbsnp_index: optional database of known germline variants
#  cosmic, cosmic_index: optional database of known somatic variants
#  is_run_orientation_bias_filter: if true, run the orientation bias filter post-processing step


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
  File? dbsnp
  File? dbsnp_index
  File? cosmic
  File? cosmic_index

  command {
    java -Xmx4g -jar ${gatk4_jar} Mutect2 \
    -R ${ref_fasta} \
    -I ${tumor_bam} \
    -I ${normal_bam} \
    -tumor ${tumor_sample_name} \
    -normal ${normal_sample_name} \
    ${"--dbsnp " + dbsnp} \
    ${"--cosmic " + cosmic} \
    ${"--normal_panel " + pon} \
    -L ${intervals} \
    -O ${tumor_sample_name}-vs-${normal_sample_name}.vcf
  }

  runtime {
    docker: "$__docker__"
    memory: "5 GB"
    disks: "local-disk " + 500 + " HDD"
  }

  output {
    File output_vcf = "${tumor_sample_name}-vs-${normal_sample_name}.vcf"
  }
}

task MergeVCFs {
  File gatk4_jar
  Array[File] input_vcfs
  String output_vcf_name

  # using MergeVcfs instead of GatherVcfs so we can create indices
  # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
  command {
    java -Xmx2g -jar ${gatk4_jar} MergeVcfs -I ${sep=' -I ' input_vcfs} -O ${output_vcf_name}.vcf
  }

  runtime {
    docker: "$__docker__"
    memory: "3 GB"
    disks: "local-disk " + 300 + " HDD"
  }

  output {
    File output_vcf = "${output_vcf_name}.vcf"
    File output_vcf_index = "${output_vcf_name}.vcf.idx"
  }
}

task Filter {
  File gatk4_jar
  File unfiltered_calls
  String output_vcf_name

  command {
  	java -Xmx4g -jar ${gatk4_jar} FilterMutectCalls -V ${unfiltered_calls} -O ${output_vcf_name}.vcf
  }

  runtime {
    docker: "$__docker__"
    memory: "5 GB"
    disks: "local-disk " + 500 + " HDD"
  }

  output {
    File output_vcf = "${output_vcf_name}.vcf"
  }
}

# Warning: this task does not work in the cloud.
task SplitIntervals {
  File gatk4_jar
  Int scatterCount
  File intervals

  command {
    mkdir intervalFileDir
    java -jar ${gatk4_jar} IntervalListTools -I ${intervals} -O intervalFileDir -SCATTER_COUNT ${scatterCount}
  }

  runtime {
    docker: "$__docker__"
    memory: "3 GB"
    disks: "local-disk " + 100 + " HDD"
  }

  output {
    Array[File] intervalFiles = glob("intervalFileDir/temp_*/scattered.intervals")
  }
}

task CollectSequencingArtifactMetrics {
    String gatk4_jar
    File bam_file
    String output_prepend
    File ref_fasta

    command {
            java -jar ${gatk4_jar} CollectSequencingArtifactMetrics \
                -I ${bam_file} -O ${output_prepend}  -R ${ref_fasta} --VALIDATION_STRINGENCY SILENT
    }

    output {
        File pre_adapter_detail_metrics = "${output_prepend}.pre_adapter_detail_metrics"
        File pre_adapter_summary_metrics = "${output_prepend}.pre_adapter_summary_metrics"
        File bait_bias_detail_metrics = "${output_prepend}.bait_bias_detail_metrics"
        File bait_bias_summary_metrics = "${output_prepend}.bait_bias_summary_metrics"
    }
}

task FilterByOrientationBias {
    String gatk4_jar
    String output_prepend
    File m2_vcf
    File pre_adapter_detail_metrics

    command {
            java -jar ${gatk4_jar} FilterByOrientationBias \
                -V ${m2_vcf} -P ${pre_adapter_detail_metrics} --output ${output_prepend}.ob_filtered.vcf
    }

    output {
        File orientation_bias_vcf = "${output_prepend}.ob_filtered.vcf"
        File orientation_bias_vcf_summary = "${output_prepend}.ob_filtered.vcf.summary"
    }
}

workflow Mutect2 {
  # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
  String gatk4_jar
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
  File? dbsnp
  File? dbsnp_index
  File? cosmic
  File? cosmic_index
  Boolean is_run_orientation_bias_filter

  call SplitIntervals {
    input:
      gatk4_jar = gatk4_jar,
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

  call MergeVCFs {
    input:
      gatk4_jar = gatk4_jar,
      input_vcfs = M2.output_vcf,
      output_vcf_name = "${tumor_sample_name}-vs-${normal_sample_name}-unfiltered"
  }

  call Filter {
    input:
      gatk4_jar = gatk4_jar,
      unfiltered_calls = MergeVCFs.output_vcf,
      output_vcf_name = "${tumor_sample_name}-vs-${normal_sample_name}-filtered"
  }

  if(is_run_orientation_bias_filter) {
      call CollectSequencingArtifactMetrics {
        input:
          gatk4_jar=gatk4_jar,
          bam_file=tumor_bam,
          ref_fasta=ref_fasta,
          output_prepend="${tumor_sample_name}",

      }

      call FilterByOrientationBias {
        input:
           gatk4_jar=gatk4_jar,
           output_prepend="${tumor_sample_name}",
           m2_vcf=Filter.output_vcf,
           pre_adapter_detail_metrics=CollectSequencingArtifactMetrics.pre_adapter_detail_metrics,
      }
  }

  output {
        File unfiltered_vcf = MergeVCFs.output_vcf
        File filtered_vcf = Filter.output_vcf
        File? ob_filtered_vcf = FilterByOrientationBias.orientation_bias_vcf
    }
}
