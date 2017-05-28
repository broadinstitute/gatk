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

# Here we implement steps 3 and 4 for several replicate each of 5-plex, 10-plex, and 20-plex pools.


import "hapmap_sensitivity.wdl" as single_plex

workflow HapmapSensitivityAllPlexes {
   # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
  	String gatk4_jar
  	Int scatter_count

  	File five_plex_bam_list
  	File ten_plex_bam_list
  	File twenty_plex_bam_list

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

    File five_plex_truth_list
    File ten_plex_truth_list
    File twenty_plex_truth_list

    String? m2_args
    String? m2_filtering_args

    File? intervals

    File python_sensitivity_script

  call single_plex.HapmapSensitivity as FivePlex {
      input:
          gatk4_jar = gatk4_jar,
          scatter_count = scatter_count,
          bam_list = five_plex_bam_list,
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
          truth_list = five_plex_truth_list,
          m2_args = m2_args,
          m2_filtering_args = m2_filtering_args,
          prefix = "5plex",
          python_sensitivity_script = python_sensitivity_script,
          intervals = intervals
  }

  call single_plex.HapmapSensitivity as TenPlex {
      input:
          gatk4_jar = gatk4_jar,
          scatter_count = scatter_count,
          bam_list = ten_plex_bam_list,
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
          truth_list = ten_plex_truth_list,
          m2_args = m2_args,
          m2_filtering_args = m2_filtering_args,
          prefix = "10plex",
          python_sensitivity_script = python_sensitivity_script,
          intervals = intervals
  }

  call single_plex.HapmapSensitivity as TwentyPlex {
      input:
          gatk4_jar = gatk4_jar,
          scatter_count = scatter_count,
          bam_list = twenty_plex_bam_list,
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
          truth_list = twenty_plex_truth_list,
          m2_args = m2_args,
          m2_filtering_args = m2_filtering_args,
          prefix = "20plex",
          python_sensitivity_script = python_sensitivity_script,
          intervals = intervals
  }

  output {
      File snp_table_5_plex = FivePlex.snp_table
      File snp_plot_5_plex = FivePlex.snp_plot
      File indel_table_5_plex = FivePlex.indel_table
      File indel_plot_5_plex = FivePlex.indel_plot
      File summary_5_plex = FivePlex.summary
      File jaccard_5_plex = FivePlex.jaccard_table
      Array[File] tpfn_5_plex = FivePlex.tpfn
      Array[File] tpfn_idx_5_plex = FivePlex.tpfn_idx
      Array[File] ftnfn_5_plex = FivePlex.ftnfn
      Array[File] ftnfn_idx_5_plex = FivePlex.ftnfn_idx
      File snp_table_10_plex = TenPlex.snp_table
      File snp_plot_10_plex = TenPlex.snp_plot
      File indel_table_10_plex = TenPlex.indel_table
      File indel_plot_10_plex = TenPlex.indel_plot
      File summary_10_plex = TenPlex.summary
      File jaccard_10_plex = TenPlex.jaccard_table
      Array[File] tpfn_10_plex = TenPlex.tpfn
      Array[File] tpfn_idx_10_plex = TenPlex.tpfn_idx
      Array[File] ftnfn_10_plex = TenPlex.ftnfn
      Array[File] ftnfn_idx_10_plex = TenPlex.ftnfn_idx
      File snp_table_20_plex = TwentyPlex.snp_table
      File snp_plot_20_plex = TwentyPlex.snp_plot
      File indel_table_20_plex = TwentyPlex.indel_table
      File indel_plot_20_plex = TwentyPlex.indel_plot
      File summary_20_plex = TwentyPlex.summary
      File jaccard_20_plex = TwentyPlex.jaccard_table
      Array[File] tpfn_20_plex = TwentyPlex.tpfn
      Array[File] tpfn_idx_20_plex = TwentyPlex.tpfn_idx
      Array[File] ftnfn_20_plex = TwentyPlex.ftnfn
      Array[File] ftnfn_idx_20_plex = TwentyPlex.ftnfn_idx
  }
}