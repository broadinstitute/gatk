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

# Here we implement these steps, except for the subsampling in step 1, for several replicate bams each of 5-plex, 10-plex, and 20-plex mixtures
# e.g. 4 different 10-plex bams, 3 10-plex bams, and 7 20-plex bams.

import "hapmap_sensitivity.wdl" as single_plex

workflow HapmapSensitivityAllPlexes {
    Int max_depth
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
    File gatk
    Array[String] artifact_modes
    File picard_jar

    File five_plex_preprocessed
    File five_plex_preprocessed_idx
    File ten_plex_preprocessed
    File ten_plex_preprocessed_idx
    File twenty_plex_preprocessed
    File twenty_plex_preprocessed_idx

    String? m2_extra_args
    String? m2_extra_filtering_args

    File? intervals

    File python_script

  call single_plex.HapmapSensitivity as FivePlex {
      input:
          max_depth = max_depth,
          scatter_count = scatter_count,
          bam_list = five_plex_bam_list,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          pon = pon,
          pon_index = pon_index,
          is_run_orientation_bias_filter = is_run_orientation_bias_filter,
          gatk = gatk,
          preemptible_attempts = preemptible_attempts,
          artifact_modes = artifact_modes,
          picard_jar = picard_jar,
          preprocessed_hapmap = five_plex_preprocessed,
          preprocessed_hapmap_idx = five_plex_preprocessed_idx,
          m2_extra_args = m2_extra_args,
          m2_extra_filtering_args = m2_extra_filtering_args,
          prefix = "5plex",
          python_script = python_script,
          intervals = intervals
  }

  call single_plex.HapmapSensitivity as TenPlex {
      input:
          max_depth = max_depth,
          scatter_count = scatter_count,
          bam_list = ten_plex_bam_list,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          pon = pon,
          pon_index = pon_index,
          is_run_orientation_bias_filter = is_run_orientation_bias_filter,
          gatk = gatk,
          preemptible_attempts = preemptible_attempts,
          artifact_modes = artifact_modes,
          picard_jar = picard_jar,
          preprocessed_hapmap = ten_plex_preprocessed,
          preprocessed_hapmap_idx = ten_plex_preprocessed_idx,
          m2_extra_args = m2_extra_args,
          m2_extra_filtering_args = m2_extra_filtering_args,
          prefix = "10plex",
          python_script = python_script,
          intervals = intervals
  }

  call single_plex.HapmapSensitivity as TwentyPlex {
      input:
          max_depth = max_depth,
          scatter_count = scatter_count,
          bam_list = twenty_plex_bam_list,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          pon = pon,
          pon_index = pon_index,
          is_run_orientation_bias_filter = is_run_orientation_bias_filter,
          gatk = gatk,
          preemptible_attempts = preemptible_attempts,
          artifact_modes = artifact_modes,
          picard_jar = picard_jar,
          preprocessed_hapmap = twenty_plex_preprocessed,
          preprocessed_hapmap_idx = twenty_plex_preprocessed_idx,
          m2_extra_args = m2_extra_args,
          m2_extra_filtering_args = m2_extra_filtering_args,
          prefix = "20plex",
          python_script = python_script,
          intervals = intervals
  }

  Array[File] all_plex_sensitivity_tables = [FivePlex.raw_table, TenPlex.raw_table, TwentyPlex.raw_table]

  call single_plex.CombineTables as AllPlexTable { input: input_tables = all_plex_sensitivity_tables, prefix = "all_plex" }

  call single_plex.AnalyzeSensitivity as AllPlex {
    input: input_table = AllPlexTable.table, python_script = python_script, prefix = "all_plex"
  }

  output {
      File snp_table_5_plex = FivePlex.snp_table
      File snp_plot_5_plex = FivePlex.snp_plot
      File indel_table_5_plex = FivePlex.indel_table
      File indel_plot_5_plex = FivePlex.indel_plot
      File summary_5_plex = FivePlex.summary
      File indel_jaccard_5_plex = FivePlex.indel_jaccard_table
      File snp_jaccard_5_plex = FivePlex.snp_jaccard_table
      Array[File] tpfn_5_plex = FivePlex.tpfn
      Array[File] tpfn_idx_5_plex = FivePlex.tpfn_idx
      Array[File] ftnfn_5_plex = FivePlex.ftnfn
      Array[File] ftnfn_idx_5_plex = FivePlex.ftnfn_idx
      File snp_table_10_plex = TenPlex.snp_table
      File snp_plot_10_plex = TenPlex.snp_plot
      File indel_table_10_plex = TenPlex.indel_table
      File indel_plot_10_plex = TenPlex.indel_plot
      File summary_10_plex = TenPlex.summary
      File indel_jaccard_10_plex = TenPlex.indel_jaccard_table
      File snp_jaccard_10_plex = TenPlex.snp_jaccard_table
      Array[File] tpfn_10_plex = TenPlex.tpfn
      Array[File] tpfn_idx_10_plex = TenPlex.tpfn_idx
      Array[File] ftnfn_10_plex = TenPlex.ftnfn
      Array[File] ftnfn_idx_10_plex = TenPlex.ftnfn_idx
      File snp_table_20_plex = TwentyPlex.snp_table
      File snp_plot_20_plex = TwentyPlex.snp_plot
      File indel_table_20_plex = TwentyPlex.indel_table
      File indel_plot_20_plex = TwentyPlex.indel_plot
      File summary_20_plex = TwentyPlex.summary
      File indel_jaccard_20_plex = TwentyPlex.indel_jaccard_table
      File snp_jaccard_20_plex = TwentyPlex.snp_jaccard_table
      Array[File] tpfn_20_plex = TwentyPlex.tpfn
      Array[File] tpfn_idx_20_plex = TwentyPlex.tpfn_idx
      Array[File] ftnfn_20_plex = TwentyPlex.ftnfn
      Array[File] ftnfn_idx_20_plex = TwentyPlex.ftnfn_idx

      File snp_table_all_plex = AllPlex.snp_table
      File snp_plot_all_plex = AllPlex.snp_plot
      File indel_table_all_plex = AllPlex.indel_table
      File indel_plot_all_plex = AllPlex.indel_plot
  }
}
