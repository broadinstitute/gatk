import "mutect2_multi_sample.wdl" as m2_multi

workflow Mutect2_Panel {
    # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
	String gatk4_jar
	Int scatter_count
	File normals_list
	File? intervals
	File ref_fasta
	File ref_fasta_index
	File ref_dict
    String m2_docker
    File? gatk4_jar_override
    Int preemptible_attempts
    File picard_jar
    String? m2_extra_args
    String pon_name

    call m2_multi.Mutect2_Multi {
        input:
            gatk4_jar = gatk4_jar,
            scatter_count = scatter_count,
            pair_list = normals_list,
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            is_run_orientation_bias_filter = false,
            is_run_oncotator = false,
            m2_docker = m2_docker,
            oncotator_docker = m2_docker,   #unused dummy value
            gatk4_jar_override = gatk4_jar_override,
            preemptible_attempts = preemptible_attempts,
            artifact_modes = ["G/T"],   #unused dummy value
            picard_jar = picard_jar,
            m2_extra_args = m2_extra_args
    }

    call CreatePanel {
        input:
            gatk4_jar = gatk4_jar,
            input_vcfs = Mutect2_Multi.unfiltered_vcf_files,
            input_vcfs_idx = Mutect2_Multi.unfiltered_vcf_index_files,
            output_vcf_name = pon_name,
            gatk4_jar_override = gatk4_jar_override,
            preemptible_attempts = preemptible_attempts,
            m2_docker = m2_docker
    }

    output {
        File pon = CreatePanel.output_vcf
        File pon_idx = CreatePanel.output_vcf_index
        Array[File] normal_calls = Mutect2_Multi.unfiltered_vcf_files
        Array[File] normal_calls_idx = Mutect2_Multi.unfiltered_vcf_index_files
    }
}

task CreatePanel {
      String gatk4_jar
      Array[File] input_vcfs
      Array[File] input_vcfs_idx
      String output_vcf_name
      File? gatk4_jar_override
      Int preemptible_attempts
      String m2_docker

      command {
            # Use GATK Jar override if specified
            GATK_JAR=${gatk4_jar}
            if [[ "${gatk4_jar_override}" == *.jar ]]; then
                GATK_JAR=${gatk4_jar_override}
            fi

            java -Xmx2g -jar $GATK_JAR CreateSomaticPanelOfNormals -vcfs ${sep=' -vcfs ' input_vcfs} -O ${output_vcf_name}.vcf
      }

      runtime {
            docker: "${m2_docker}"
            memory: "5 GB"
            disks: "local-disk " + 300 + " HDD"
            preemptible: "${preemptible_attempts}"
      }

      output {
            File output_vcf = "${output_vcf_name}.vcf"
            File output_vcf_index = "${output_vcf_name}.vcf.idx"
      }
}