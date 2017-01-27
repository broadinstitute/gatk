#  Run Mutect 2 on a list of tumor-normal pairs
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

import "mutect2.wdl" as m2


task CreateOutputList {
    String output_name
	Array[File] vcfs

	command {
		for vcf in ${sep=" " vcfs}; do
			echo $vcf
			echo $vcf >> ${output_name}.list
		done
	}

	output {
		File vcf_list = "${output_name}.list"
	}
}

workflow Mutect2 {
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

	scatter(pair in pairs) {
		call m2.Mutect2 {
			input:
				gatk4_jar = gatk4_jar,
				intervals = intervals,
				ref_fasta = ref_fasta,
				ref_fasta_index = ref_fasta_index,
				ref_dict = ref_dict,
				tumor_bam = pair[0],
				tumor_bam_index = pair[1],
				tumor_sample_name = pair[2],
				normal_bam = pair[3],
				normal_bam_index = pair[4],
				normal_sample_name = pair[5],
				pon = pon,
				pon_index = pon_index,
				scatter_count = scatter_count,
				dbsnp = dbsnp,
				dbsnp_index = dbsnp_index,
				cosmic = cosmic,
				cosmic_index = cosmic_index,
				is_run_orientation_bias_filter = is_run_orientation_bias_filter
		}
	}

	call CreateOutputList as unfilteredOutputList {
		input:
		    output_name = "unfiltered",
			vcfs = Mutect2.unfiltered_vcf
	}

	call CreateOutputList as filteredOutputList {
        input:
    	    output_name = "filtered",
    	    vcfs = Mutect2.filtered_vcf
    }

    if(is_run_orientation_bias_filter) {
        call CreateOutputList as orientationBiasFilteredOutputList {
            input:
                output_name = "ob_filtered",
                vcfs = select_first(Mutect2.ob_filtered_vcf)
        }
    }

    output {
        File unfiltered_vcfs = unfilteredOutputList.vcf_list
        File filtered_vcfs = filteredOutputList.vcf_list
        File? ob_filtered_vcfs = orientationBiasFilteredOutputList.vcf_list
    }
}