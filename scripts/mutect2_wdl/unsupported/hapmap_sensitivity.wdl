### This is the concordance part of the CRSP sensitivity validation

#Conceptual Overview
# To measure sensitivity, we sequence a pool 5, 10 or 20 normal samples. The pool contains a variety of allele fractions,
# like a real tumor sample. These normal samples were also sequenced individually as part of HapMap, so we have "truth" vcfs
# for them.  We calculate sensitivity by comparing the calls to the truth data.  For a variety of reasons we don't call
# against a matched normal.

#Workflow Steps
# 1.  Restrict a huge HapMap wgs VCF to given lists of samples and intervals.
# 2.  Annotate this VCF using a bam sequenced from a pool derived from the given samples
# 3.  Run Mutect in tumor-only mode on this pooled bam.
# 4.  Compare Mutect calls to the truth data and output a table of true positives and false negatives along with
#     annotations from the truth VCF prepared in steps 1 and 2.

# Here we implement steps 3 and 4.

#Inputs are identical to mutect2_multi_sample_concordance.wdl except for the addition of a python script for generating plots

import "mutect2_multi_sample_concordance.wdl" as concordance


task ConvertToTable {
  String gatk4_jar
  File? gatk4_jar_override
  File input_vcf
  File input_vcf_idx

  command {
  # Use GATK Jar override if specified
            GATK_JAR=${gatk4_jar}
            if [[ "${gatk4_jar_override}" == *.jar ]]; then
                GATK_JAR=${gatk4_jar_override}
            fi

        java -jar $GATK_JAR VariantsToTable -V ${input_vcf} -F STATUS -F BAM_DEPTH -F AF_EXP -F TYPE -O "result.table"
  }

  runtime {
          memory: "5 GB"
      }

  output {
    File table = "result.table"
  }
}

task AnalyzeSensitivity {
    File input_table
    File python_sensitivity_script

    command {
        python ${python_sensitivity_script} ${input_table}
    }

    runtime {
            memory: "5 GB"
        }

    output {
        File snp_table = "SNP_sensitivity.tsv"
        File snp_plot = "SNP_sensitivity.png"
        File indel_table = "Indel_sensitivity.tsv"
        File indel_plot = "Indel_sensitivity.png"
    }
}


workflow HapmapSensitivity {
   # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
  	String gatk4_jar
  	Int scatter_count
  	File pair_list
  	File? intervals
  	File ref_fasta
  	File ref_fasta_index
  	File ref_dict
  	File? pon
  	File? pon_index
  	Boolean is_run_orientation_bias_filter
    String m2_docker
    File? gatk4_jar_override
    Int preemptible_attempts
    Array[String] artifact_modes
    File picard_jar
    File truth_list

    File python_sensitivity_script

  call concordance.Mutect2_Multi_Concordance {
    input:
     # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
    	gatk4_jar = gatk4_jar,
    	scatter_count = scatter_count,
    	pair_list = pair_list,
    	intervals = intervals,
    	ref_fasta = ref_fasta,
    	ref_fasta_index = ref_fasta_index,
    	ref_dict = ref_dict,
    	pon = pon,
    	pon_index = pon_index,
    	is_run_orientation_bias_filter = is_run_orientation_bias_filter,
        m2_docker = m2_docker,
        gatk4_jar_override = gatk4_jar_override,
        preemptible_attempts = preemptible_attempts,
        artifact_modes = artifact_modes,
        picard_jar = picard_jar,
        truth_list = truth_list
  }

  scatter(n in range(length(read_tsv(pair_list)))) {
    call ConvertToTable {
        input:
              gatk4_jar = gatk4_jar,
              gatk4_jar_override = gatk4_jar_override,
              input_vcf = Mutect2_Multi_Concordance.tpfn[n],
              input_vcf_idx = Mutect2_Multi_Concordance.tpfn_idx[n]
      }

      call AnalyzeSensitivity {
           input:
            input_table = ConvertToTable.table,
            python_sensitivity_script = python_sensitivity_script
      }
  }

  output {
    Array[File] snp_table = AnalyzeSensitivity.snp_table
    Array[File] snp_plot = AnalyzeSensitivity.snp_plot
    Array[File] indel_table = AnalyzeSensitivity.indel_table
    Array[File] indel_plot = AnalyzeSensitivity.indel_plot
    Array[File] summary = Mutect2_Multi_Concordance.summary
    Array[File] tpfn = Mutect2_Multi_Concordance.tpfn
    Array[File] tpfn_idx = Mutect2_Multi_Concordance.tpfn_idx
  }
}