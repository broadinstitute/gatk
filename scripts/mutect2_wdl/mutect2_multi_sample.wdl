import "mutect2.wdl" as m2


task CreateOutputList {
    String output_name
	Array[File] vcfs

	command <<<
		for vcf in ${sep=" " vcfs}; do
			echo $vcf
			echo $vcf >> ${output_name}.list
		done
	>>>

	output {
		File vcf_list = '${output_name}.list'
	}
}

workflow Mutect2 {
	File gatk4_jar
	Int scatter_count
	# pair_list file is a tsv file with the following six columns in this order.
	# tumor_bam, tumor_bam_index, tumor_sample_name, normal_bam, normal_bam_index, normal_sample_name
	File pair_list
	Array[Array[String]] pairs = read_tsv(pair_list)
	File intervals
	File ref_fasta
	File ref_fasta_index
	File ref_dict
	File pon
	File pon_index
	File dbsnp
	File dbsnp_index
	File cosmic
	File cosmic_index
	File picard_jar

	scatter(pair in pairs) {
		call m2.Mutect2 {
			input:
				gatk4_jar=gatk4_jar,
				intervals=intervals,
				ref_fasta=ref_fasta,
				ref_fasta_index=ref_fasta_index,
				ref_dict=ref_dict,
				tumor_bam=pair[0],
				tumor_bam_index=pair[1],
				tumor_sample_name=pair[2],
				normal_bam=pair[3],
				normal_bam_index=pair[4],
				normal_sample_name=pair[5],
				pon=pon,
				pon_index=pon_index,
				scatter_count=scatter_count,
				dbsnp=dbsnp,
				dbsnp_index=dbsnp_index,
				cosmic=cosmic,
				cosmic_index=cosmic_index,
				picard_jar=picard_jar
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
}