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

import "mutect2.wdl" as m2

task Concordance {
  File gatk_jar
  File? intervals
  File truth_vcf
  File truth_vcf_idx
  File eval_vcf
  File eval_vcf_idx

  command {
      java -jar ${gatk_jar} Concordance ${"-L " + intervals} \
        -truth ${truth_vcf} -eval ${eval_vcf} -tpfn "true_positives_and_false_negatives.vcf" \
        -summary summary.tsv
  }

    runtime {
        memory: "5 GB"
    }

  output {
        File output_vcf = "true_positives_and_false_negatives.vcf"
        File output_vcf_idx = "true_positives_and_false_negatives.vcf.idx"
        File summary = "summary.tsv"
  }
}

task ConvertToTable {
  File gatk_jar
  File input_vcf
  File input_vcf_idx

  command {
        java -jar ${gatk_jar} VariantsToTable -V ${input_vcf} -F STATUS -F BAM_DEPTH -F AF_EXP -F TYPE -O "result.table"
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
        File indel_table = "INDEL_sensitivity.tsv"
        File indel_plot = "INDEL_sensitivity.png"
    }
}


workflow HapmapSensitivity {
  #general
  File gatk_jar
  File intervals

  #for running Mutect
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? pon
  File? pon_index
  Int scatter_count
  Boolean is_run_orientation_bias_filter
  Array[String] artifact_modes

  #hapmap
  #the output of hapmap_sensitivity_truth.wdl
  File hapmap_truth_vcf
  File hapmap_truth_vcf_idx
  File pooled_bam
  File pooled_bam_index
  String pooled_sample_name

  File python_sensitivity_script

  call m2.Mutect2 {
    input:
      gatk4_jar = "DUMMY",
      intervals = intervals,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      tumor_bam = pooled_bam,
      tumor_bam_index = pooled_bam_index,
      tumor_sample_name = pooled_sample_name,
      pon = pon,
      pon_index = pon_index,
      scatter_count = scatter_count,
      is_run_orientation_bias_filter = is_run_orientation_bias_filter,
      is_run_oncotator = false,
      m2_docker = "broadinstitute/gatk-protected:1.0.0.0-alpha1.2.4",
      oncotator_docker = "DUMMY",
      gatk4_jar_override = gatk_jar,
      preemptible_attempts = 2,
      artifact_modes = artifact_modes
  }

  call Concordance {
    input:
        gatk_jar = gatk_jar,
        intervals = intervals,
        truth_vcf = hapmap_truth_vcf,
        truth_vcf_idx = hapmap_truth_vcf_idx,
        eval_vcf = Mutect2.filtered_vcf,
        eval_vcf_idx = Mutect2.filtered_vcf_index
  }

  call ConvertToTable {
    input:
          gatk_jar = gatk_jar,
          input_vcf = Concordance.output_vcf,
          input_vcf_idx = Concordance.output_vcf_idx
  }

  call AnalyzeSensitivity {
       input:
        input_table = ConvertToTable.table,
        python_sensitivity_script = python_sensitivity_script
  }

  output {
    File snp_table = AnalyzeSensitivity.snp_table
    File snp_plot = AnalyzeSensitivity.snp_plot
    File indel_table = AnalyzeSensitivity.indel_table
    File indel_plot = AnalyzeSensitivity.indel_plot
  }
}