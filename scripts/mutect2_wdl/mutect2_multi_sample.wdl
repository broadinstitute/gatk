version 1.0

#  Run Mutect 2 on a list of tumors or tumor-normal pairs
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  pon, pon_idx: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  gnomad, gnomad_idx: optional database of known germline variants, obtainable from http://gnomad.broadinstitute.org/downloads
#  variants_for_contamination, variants_for_contamination_idx: vcf of common variants with allele frequencies fo calculating contamination
#  run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  pair_list: a tab-separated table with no header in the following format:
#   TUMOR_1_BAM</TAB>TUMOR_1_bai</TAB>NORMAL_1_BAM</TAB>NORMAL_1_bai
#   TUMOR_2_BAM</TAB>TUMOR_2_bai</TAB>NORMAL_2_BAM</TAB>NORMAL_2_bai
#   . . .
#  Tumor-only input is the same but without the columns for the normal:
#  TUMOR_1_BAM</TAB>TUMOR_1_bai
#  TUMOR_2_BAM</TAB>TUMOR_2_bai
#   . . .

import "mutect2.wdl" as m2

workflow Mutect2_Multi {
  input {
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File pair_list

    File? pon
    File? pon_idx
    File? gnomad
    File? gnomad_idx
    File? variants_for_contamination
    File? variants_for_contamination_idx
    Boolean? run_orientation_bias_mixture_model_filter
    Int scatter_count
    String? m2_extra_args
    String? m2_extra_filtering_args
    Boolean? compress_vcfs
    Boolean? make_bamout

    String? gcs_project_for_requester_pays

    # Oncotator inputs
    String? sequencing_center
    String? sequence_source

    # funcotator inputs
    Boolean? run_funcotator
    String? funco_reference_version
    File? funco_data_sources_tar_gz
    String? funco_transcript_selection_mode
    File? funco_transcript_selection_list
    Array[String]? funco_annotation_defaults
    Array[String]? funco_annotation_overrides


    # runtime
    String gatk_docker
    Int? preemptible_attempts
    File? gatk_override
  }

  Array[Array[String]] pairs = read_tsv(pair_list)

	scatter( row in pairs ) {
	    #      If the condition is true, variables inside the 'if' block retain their values outside the block.
	    #      Otherwise they are treated as null, which in WDL is equivalent to an empty optional
        if(length(row) == 4) {
            File normal_bam = row[2]
            File normal_bai = row[3]
        }

        call m2.Mutect2 {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_reads = row[0],
                tumor_reads_index = row[1],
                normal_reads = normal_bam,
                normal_reads_index = normal_bai,
                pon = pon,
                pon_idx = pon_idx,
                scatter_count = scatter_count,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                variants_for_contamination = variants_for_contamination,
                variants_for_contamination_idx = variants_for_contamination_idx,
                run_orientation_bias_mixture_model_filter = run_orientation_bias_mixture_model_filter,
                m2_extra_args = m2_extra_args,
                m2_extra_filtering_args = m2_extra_filtering_args,
                sequencing_center = sequencing_center,
                sequence_source = sequence_source,
                run_funcotator = run_funcotator,
                funco_reference_version = funco_reference_version,
                funco_data_sources_tar_gz = funco_data_sources_tar_gz,
                funco_transcript_selection_mode = funco_transcript_selection_mode,
                funco_transcript_selection_list = funco_transcript_selection_list,
                funco_annotation_defaults = funco_annotation_defaults,
                funco_annotation_overrides = funco_annotation_overrides,

                make_bamout = make_bamout,
                compress_vcfs = compress_vcfs,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible = preemptible_attempts,
                gcs_project_for_requester_pays = gcs_project_for_requester_pays
        }
    }

    output {
        Array[File] filtered_vcf = Mutect2.filtered_vcf
        Array[File] filtered_vcf_idx = Mutect2.filtered_vcf_idx
        Array[File?] contamination_tables = Mutect2.contamination_table

        Array[File?] m2_bamout = Mutect2.bamout
        Array[File?] m2_bamout_index = Mutect2.bamout_index
    }
}
