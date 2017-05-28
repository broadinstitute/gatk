### This is the truth data preparation part of the CRSP sensitivity validation

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

# Here we implement steps 1 and 2 for replicate sets of 5-plex, 10-plex, and 20-plex pools, e.g. 7 replicates each of
# 5-plex, 10-plex, and 20-plex

import "hapmap_sensitivity_truth.wdl" as single_plex

workflow HapmapSensitivityTruthAllPlexes {
    String gatk_jar
    File hapmap
    File hapmap_idx

    File five_plex_pooled_bams_list
    File ten_plex_pooled_bams_list
    File twenty_plex_pooled_bams_list

    File five_plex_sample_file
    File ten_plex_sample_file
    File twenty_plex_sample_file

    File? intervals
    Int min_indel_spacing

    File common_vcf
    File common_vcf_idx

    Int max_depth

    call single_plex.HapmapSensitivityTruth as FivePlex {
        input:
          gatk_jar = gatk_jar,
          hapmap = hapmap,
          hapmap_idx = hapmap_idx,
          pooled_bams_list = five_plex_pooled_bams_list,
          sample_file = five_plex_sample_file,
          intervals = intervals,
          min_indel_spacing = min_indel_spacing,
          common_vcf = common_vcf,
          common_vcf_idx = common_vcf_idx,
          max_depth = max_depth
    }

    call single_plex.HapmapSensitivityTruth as TenPlex {
        input:
          gatk_jar = gatk_jar,
          hapmap = hapmap,
          hapmap_idx = hapmap_idx,
          pooled_bams_list = ten_plex_pooled_bams_list,
          sample_file = ten_plex_sample_file,
          intervals = intervals,
          min_indel_spacing = min_indel_spacing,
          common_vcf = common_vcf,
          common_vcf_idx = common_vcf_idx,
          max_depth = max_depth
    }

    call single_plex.HapmapSensitivityTruth as TwentyPlex {
        input:
          gatk_jar = gatk_jar,
          hapmap = hapmap,
          hapmap_idx = hapmap_idx,
          pooled_bams_list = twenty_plex_pooled_bams_list,
          sample_file = twenty_plex_sample_file,
          intervals = intervals,
          min_indel_spacing = min_indel_spacing,
          common_vcf = common_vcf,
          common_vcf_idx = common_vcf_idx,
          max_depth = max_depth
    }

    output {
        Array[File] five_plex_mixing_fractions_table = FivePlex.mixing_fractions_table
        Array[File] five_plex_output_vcf = FivePlex.output_vcf
        Array[File] ten_plex_mixing_fractions_table = TenPlex.mixing_fractions_table
        Array[File] ten_plex_output_vcf = TenPlex.output_vcf
        Array[File] twenty_plex_mixing_fractions_table = TwentyPlex.mixing_fractions_table
        Array[File] twenty_plex_output_vcf = TwentyPlex.output_vcf
    }

}