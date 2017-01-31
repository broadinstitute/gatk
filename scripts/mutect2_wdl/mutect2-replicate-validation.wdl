# call pairs of replicates as a tumor-normal pair, count the number of variants (i.e. false positives) 
# and report false positive rates

import "mutect2.wdl" as m2

task GatherTables {
    # we assume that each table consists of two lines: one header line and one record
	Array[File] tables

	command {
	    # extract the header from one of the files
		head -n 1 ${tables[0]} > summary.txt

		# then append the record from each table
		for table in ${sep=" " tables}; do
			tail -n +2 $table >> summary.txt						
		done
	}

	runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        memory: "1 GB"
        disks: "local-disk " + 100 + " HDD"
    }

	output {
		File summary = "summary.txt"
	}
}

task CountFalsePositives {
	File gatk4_jar
	File filtered_vcf
	File filtered_vcf_index
	File ref_fasta
	File ref_fasta_index
	File ref_dict
	File intervals

	command {
		java -jar ${gatk4_jar} CountFalsePositives \
		-V ${filtered_vcf} \
		-R ${ref_fasta} \
		-L ${intervals} \
		-O false-positives.txt \
	}

    runtime {
        docker: "broadinstitute/genomes-in-the-cloud:2.2.4-1469632282"
        memory: "5 GB"
        disks: "local-disk " + 500 + " HDD"
    }

	output {
		File false_positive_counts = "false-positives.txt"
	}
}

workflow Mutect2ReplicateValidation {
	File gatk4_jar
	Int scatter_count
	# replicate_pair_list file is a tsv file with the following six columns in this order.
	# tumor_bam, tumor_bam_index, tumor_sample_name, normal_bam, normal_bam_index, normal_sample_name
	File replicate_pair_list
	Array[Array[String]] pairs = read_tsv(replicate_pair_list)
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
	Boolean run_orientation_bias_filter

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
				is_run_orientation_bias_filter=run_orientation_bias_filter
		}

		call CountFalsePositives {
			input:
				gatk4_jar=gatk4_jar,
				filtered_vcf=Mutect2.filtered_vcf,
				filtered_vcf_index=Mutect2.filtered_vcf_index,
				ref_fasta=ref_fasta,
				ref_fasta_index=ref_fasta_index,
				ref_dict=ref_dict,
				intervals=intervals
		}
	}

	call GatherTables {
        input:
            tables=CountFalsePositives.false_positive_counts
	}

	output {
		File summary = GatherTables.summary
	}
}
