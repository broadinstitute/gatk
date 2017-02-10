#  Run Mutect 2 on a list of tumors or tumor-normal pairs
#
#  Description of inputs
#  gatk4_jar: java jar file containing gatk 4 (protected)
#  intervals: genomic intervals
#  ref_fasta, ref_fasta_index, ref_dict: reference genome, index, and dictionary
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  dbsnp, dbsnp_index: optional database of known germline variants
#  cosmic, cosmic_index: optional database of known somatic variants
#  is_run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  pair_list: a tab-separated table with no header in the following format:
#   TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE</TAB>NORMAL_1_BAM</TAB>NORMAL_1_BAM_INDEX</TAB>NORMAL_1_SAMPLE</TAB>
#   TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE</TAB>NORMAL_2_BAM</TAB>NORMAL_2_BAM_INDEX</TAB>NORMAL_2_SAMPLE</TAB>
#   . . .
#  Temporarily, while waiting for a Cromwell bug to be resolved, tumor-only input looks like
#  TUMOR_1_BAM</TAB>TUMOR_1_BAM_INDEX</TAB>TUMOR_1_SAMPLE</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>
#  TUMOR_2_BAM</TAB>TUMOR_2_BAM_INDEX</TAB>TUMOR_2_SAMPLE</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>NO_NORMAL</TAB>
#   . . .
#  That is, you actually write out "NO NORMAL" thrice per row

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
        docker: "ubuntu:14.04"
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
	Boolean is_run_orientation_bias_filter
	Boolean is_run_oncotator
    String m2_docker
    String oncotator_docker
    File? gatk4_jar_override
    Int preemptible_attempts
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    Array[String] artifact_modes

	scatter( row in pairs ) {
	    #The non-hack way, but there's a bug
	    #      In WDL, variables inside the block can be used outside the block.
	    #      If the conditional block is run, they retain their values.
	    #      Otherwise, they evaluate to null, which in WDL is equivalent to an empty optional
	    #      If we simply tried to use eg row[3] below it could cause an out-of-bounds exception
        #if(length(pairs[n]) == 6) {
        #    File normal_bam = row[3]
        #    File normal_bam_index = row[4]
        #    String normal_sample_name = row[5]
        #}

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
                    normal_bam = row[3],
                    normal_bam_index = row[4],
                    normal_sample_name = row[5],
                    pon = pon,
                    pon_index = pon_index,
                    scatter_count = scatter_count,
                    dbsnp = dbsnp,
                    dbsnp_index = dbsnp_index,
                    cosmic = cosmic,
                    cosmic_index = cosmic_index,
                    picard_jar = picard_jar,
                    is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                    is_run_oncotator=is_run_oncotator,
                    oncotator_docker=oncotator_docker,
                    m2_docker = m2_docker,
                    gatk4_jar_override = gatk4_jar_override,
                    preemptible_attempts = preemptible_attempts,
                    onco_ds_tar_gz = onco_ds_tar_gz,
                    onco_ds_local_db_dir = onco_ds_local_db_dir,
                    artifact_modes = artifact_modes
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

    if(is_run_orientation_bias_filter) {
        call CreateOutputList as orientationBiasFilteredOutputList {
            input:
                output_name = "ob_filtered",
                vcfs = select_all(Mutect2.ob_filtered_vcf),
                preemptible_attempts = preemptible_attempts
        }
    }

    output {
        File unfiltered_vcfs = unfilteredOutputList.vcf_list
        File filtered_vcfs = filteredOutputList.vcf_list
        File? ob_filtered_vcfs = orientationBiasFilteredOutputList.vcf_list

        Array[File] unfiltered_vcf_files = Mutect2.unfiltered_vcf
        Array[File] filtered_vcf_files = Mutect2.filtered_vcf
        Array[File?] ob_filtered_vcf_files = Mutect2.ob_filtered_vcf
        Array[File?] oncotated_m2_vcf_files = Mutect2.oncotated_m2_vcf
    }
}