#  Run Mutect 2 on a list of tumors or tumor-normal pairs
#
#  Description of inputs
#  gatk4_jar: java jar file containing gatk 4 (protected)
#  intervals: genomic intervals
#  ref_fasta, ref_fasta_index, ref_dict: reference genome, index, and dictionary
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  gnomad, gnomad_index: optional database of known germline variants, obtainable from http://gnomad.broadinstitute.org/downloads
#  variants_for_contamination, variants_for_contamination_index: vcf of common variants with allele frequencies fo calculating contamination
#  is_run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  pair_list: a tab-separated table with no header in the following format:
#   TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE</TAB>NORMAL_1_BAM</TAB>NORMAL_1_BAM_INDEX</TAB>NORMAL_1_SAMPLE
#   TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE</TAB>NORMAL_2_BAM</TAB>NORMAL_2_BAM_INDEX</TAB>NORMAL_2_SAMPLE
#   . . .
#  Tumor-only input is the same but without the columns for the normal:
#  TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE
#  TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE
#   . . .

import "mutect2.wdl" as m2


#
# IMPORTANT: This task will not generate useful results for any backends using docker (incl. JES/cloud).
#
task CreateOutputList {
    String output_name
	Array[String] vcfs
    Int preemptible_attempts

	command {
		for vcf in ${sep=" " vcfs}; do
			echo $vcf
			echo $vcf >> ${output_name}.list
		done
	}

	runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        memory: "1 GB"
        preemptible: "${preemptible_attempts}"
	}

	output {
		File vcf_list = "${output_name}.list"
	}
}


workflow Mutect2_Multi {
    # gatk4_jar needs to be a String input to the workflow in order to work in a Docker image
	String gatk4_jar
	Int scatter_count
	File pair_list
	Array[Array[String]] pairs = read_tsv(pair_list)
	File? intervals
	File ref_fasta
	File ref_fasta_index
	File ref_dict
	File? pon
	File? pon_index
	File? gnomad
	File? gnomad_index
	File? variants_for_contamination
    File? variants_for_contamination_index
	Boolean is_run_orientation_bias_filter
	Boolean is_run_oncotator
    String m2_docker
    String oncotator_docker
    File? gatk4_jar_override
    Int preemptible_attempts
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    Array[String] artifact_modes
    File picard_jar

	scatter( row in pairs ) {
	    #      If the condition is true, variables inside the 'if' block retain their values outside the block.
	    #      Otherwise they are treated as null, which in WDL is equivalent to an empty optional
        if(length(row) == 6) {
            File normal_bam = row[3]
            File normal_bam_index = row[4]
            String normal_sample_name = row[5]
        }

            call m2.Mutect2 {
                input:
                    gatk4_jar = gatk4_jar,
                    intervals = intervals,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    ref_dict = ref_dict,
                    tumor_bam = row[0],
                    tumor_bam_index = row[1],
                    tumor_sample_name = row[2],
                    normal_bam = normal_bam,
                    normal_bam_index = normal_bam_index,
                    normal_sample_name = normal_sample_name,
                    pon = pon,
                    pon_index = pon_index,
                    scatter_count = scatter_count,
                    gnomad = gnomad,
                    gnomad_index = gnomad_index,
                    variants_for_contamination = variants_for_contamination,
                    variants_for_contamination_index = variants_for_contamination_index,
                    is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                    is_run_oncotator = is_run_oncotator,
                    oncotator_docker = oncotator_docker,
                    m2_docker = m2_docker,
                    gatk4_jar_override = gatk4_jar_override,
                    preemptible_attempts = preemptible_attempts,
                    onco_ds_tar_gz = onco_ds_tar_gz,
                    onco_ds_local_db_dir = onco_ds_local_db_dir,
                    artifact_modes = artifact_modes,
                    picard_jar = picard_jar
            }
    }


	call CreateOutputList as unfilteredOutputList {
		input:
		    output_name = "unfiltered",
			vcfs = Mutect2.unfiltered_vcf,
			preemptible_attempts = preemptible_attempts
	}

	call CreateOutputList as filteredOutputList {
        input:
    	    output_name = "filtered",
    	    vcfs = Mutect2.filtered_vcf,
    	    preemptible_attempts = preemptible_attempts
    }

    output {
        File unfiltered_vcfs = unfilteredOutputList.vcf_list
        File filtered_vcfs = filteredOutputList.vcf_list
        Array[File] unfiltered_vcf_files = Mutect2.unfiltered_vcf
        Array[File] filtered_vcf_files = Mutect2.filtered_vcf
        Array[File] filtered_vcf_index_files = Mutect2.filtered_vcf_index
    }
}