#Inputs are same as mutect2_multi_sample.wdl with the addtion of a truth_list input, which is a tsv of the form
# truth1.vcf    truth1.vcf.idx
# truth2.vcf    truth2.vcf.idx
# . . .
# The rows of this file correspond to the rows of the pair_list input

import "mutect2_multi_sample.wdl" as m2_multi

workflow Mutect2_Multi_Concordance {
    # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
	String gatk4_jar
	Int scatter_count
	File pair_list
	File intervals
	File ref_fasta
	File ref_fasta_index
	File ref_dict
	File? pon
	File? pon_index
	File? dbsnp
	File? dbsnp_index
	File? cosmic
	File? cosmic_index
	File? variants_for_contamination
    File? variants_for_contamination_index
	Boolean is_run_orientation_bias_filter
    String m2_docker
    File? gatk4_jar_override
    Int preemptible_attempts
    Array[String] artifact_modes
    File picard_jar

    File truth_list
    Array[Array[String]] truth = read_tsv(truth_list)

    call m2_multi.Mutect2_Multi {
      input:
            gatk4_jar = gatk4_jar,
        	scatter_count = scatter_count,
        	pair_list = pair_list,
        	intervals = intervals,
        	ref_fasta = ref_fasta,
        	ref_fasta_index = ref_fasta_index,
        	ref_dict = ref_dict,
        	pon = pon,
        	pon_index = pon_index,
        	dbsnp = dbsnp,
        	dbsnp_index = dbsnp_index,
        	cosmic = cosmic,
        	cosmic_index = cosmic_index,
        	variants_for_contamination = variants_for_contamination,
            variants_for_contamination_index = variants_for_contamination_index,
        	is_run_orientation_bias_filter = is_run_orientation_bias_filter,
        	is_run_oncotator = false,
            m2_docker = m2_docker,
            oncotator_docker = "NO_ONCOTATOR",
            gatk4_jar_override = gatk4_jar_override,
            preemptible_attempts = preemptible_attempts,
            artifact_modes = artifact_modes,
            picard_jar = picard_jar
    }

       scatter (n in range(length(truth))) {
            call Concordance {
                input:
                    gatk4_jar = gatk4_jar,
                    gatk4_jar_override = gatk4_jar_override,
                    intervals = intervals,
                    truth_vcf = truth[n][0],
                    truth_vcf_idx = truth[n][1],
                    eval_vcf = Mutect2_Multi.filtered_vcf_files[n],
                    eval_vcf_idx = Mutect2_Multi.filtered_vcf_index_files[n],
                    m2_docker = m2_docker,
                    preemptible_attempts = preemptible_attempts
            }
         }

    output {
      Array[File] tpfn = Concordance.tpfn
      Array[File] tpfn_idx = Concordance.tpfn_idx
      Array[File] tpfp = Concordance.tpfp
      Array[File] tpfp_idx = Concordance.tpfp_idx
      Array[File] summary = Concordance.summary
    }
}

    task Concordance {
      String gatk4_jar
      File? gatk4_jar_override
      File? intervals
      File truth_vcf
      File truth_vcf_idx
      File eval_vcf
      File eval_vcf_idx
      String m2_docker
      Int preemptible_attempts

      command {
        # Use GATK Jar override if specified
          GATK_JAR=${gatk4_jar}
          if [[ "${gatk4_jar_override}" == *.jar ]]; then
              GATK_JAR=${gatk4_jar_override}
          fi

          java -jar $GATK_JAR Concordance ${"-L " + intervals} \
            -truth ${truth_vcf} -eval ${eval_vcf} \
            -tpfn "true_positives_and_false_negatives.vcf" \
            -tpfp "true_positives_and_false_positives.vcf" \
            -summary summary.tsv
      }

        runtime {
            memory: "5 GB"
            docker: "${m2_docker}"
            disks: "local-disk " + 400 + " HDD"
             preemptible: "${preemptible_attempts}"
        }

      output {
            File tpfn = "true_positives_and_false_negatives.vcf"
            File tpfn_idx = "true_positives_and_false_negatives.vcf.idx"
            File tpfp = "true_positives_and_false_positives.vcf"
            File tpfp_idx = "true_positives_and_false_positives.vcf.idx"
            File summary = "summary.tsv"
      }
}