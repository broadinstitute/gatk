#  Run Mutect 2 on a list of tumors or tumor-normal pairs
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  gnomad, gnomad_index: optional database of known germline variants, obtainable from http://gnomad.broadinstitute.org/downloads
#  variants_for_contamination, variants_for_contamination_index: vcf of common variants with allele frequencies fo calculating contamination
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
    # Mutect2 inputs
	File? intervals
	File ref_fasta
	File ref_fai
	File ref_dict
	File pair_list
	Array[Array[String]] pairs = read_tsv(pair_list)
	File? pon
	File? pon_index
	File? gnomad
	File? gnomad_index
	File? variants_for_contamination
    File? variants_for_contamination_index
	Boolean? run_orientation_bias_filter
	Int scatter_count
	Array[String]? artifact_modes
	String? m2_extra_args
    String? m2_extra_filtering_args
    Boolean? compress_vcfs
    Boolean? make_bamout

    # Oncotator inputs
    Boolean? run_oncotator
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    String? sequencing_center
    String? sequence_source
    File? default_config_file

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
    String? oncotator_docker
    Int? preemptible_attempts
    File? gatk_override

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
                tumor_bam = row[0],
                tumor_bai = row[1],
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                pon = pon,
                pon_index = pon_index,
                scatter_count = scatter_count,
                gnomad = gnomad,
                gnomad_index = gnomad_index,
                variants_for_contamination = variants_for_contamination,
                variants_for_contamination_index = variants_for_contamination_index,
                run_orientation_bias_filter = run_orientation_bias_filter,
                artifact_modes = artifact_modes,
                m2_extra_args = m2_extra_args,
                m2_extra_filtering_args = m2_extra_filtering_args,
                run_oncotator = run_oncotator,
                onco_ds_tar_gz = onco_ds_tar_gz,
                onco_ds_local_db_dir = onco_ds_local_db_dir,
                sequencing_center = sequencing_center,
                sequence_source = sequence_source,
                default_config_file = default_config_file,
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
                oncotator_docker = oncotator_docker,
                preemptible_attempts = preemptible_attempts
        }
    }

    output {
        Array[File] filtered_vcf = Mutect2.filtered_vcf
        Array[File] filtered_vcf_idx = Mutect2.filtered_vcf_index
        Array[File?] contamination_tables = Mutect2.contamination_table

        Array[File?] oncotated_m2_mafs = Mutect2.oncotated_m2_maf
        Array[File?] m2_bamout = Mutect2.bamout
        Array[File?] m2_bamout_index = Mutect2.bamout_index
    }
}