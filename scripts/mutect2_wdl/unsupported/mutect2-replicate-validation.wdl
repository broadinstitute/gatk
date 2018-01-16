# call pairs of replicates as a tumor-normal pair, count the number of variants (i.e. false positives) 
# and report false positive rates

import "mutect2.wdl" as m2

workflow Mutect2ReplicateValidation {
	File? intervals
	File ref_fasta
	File ref_fai
	File ref_dict
	Int scatter_count
	# replicate_pair_list file is a tsv file with the following four columns in this order.
	# tumor_bam, tumor_bai, normal_bam, normal_bai
	File replicate_pair_list
	Array[Array[String]] pairs = read_tsv(replicate_pair_list)
	File? pon
	File? pon_index
	File? gnomad
	File? gnomad_index
	Boolean is_run_orientation_bias_filter
	Array[String] artifact_modes
	String? m2_extra_args
    String? m2_extra_filtering_args

    File? gatk_override

    String gatk_docker
    Int? preemptible_attempts

	scatter(pair in pairs) {
		call m2.Mutect2 {
			input:
				intervals = intervals,
				ref_fasta = ref_fasta,
				ref_fai = ref_fai,
				ref_dict = ref_dict,
				tumor_bam = pair[0],
				tumor_bai = pair[1],
				normal_bam = pair[2],
				normal_bai = pair[3],
				pon = pon,
				pon_index = pon_index,
				scatter_count = scatter_count,
				gnomad = gnomad,
				gnomad_index = gnomad_index,
                is_run_orientation_bias_filter = is_run_orientation_bias_filter,
                is_run_oncotator = false,
                preemptible_attempts = preemptible_attempts,
                artifact_modes = artifact_modes,
                m2_extra_args = m2_extra_args,
                m2_extra_filtering_args = m2_extra_filtering_args,
                gatk_override = gatk_override,
				gatk_docker = gatk_docker,
				oncotator_docker = gatk_docker
		}

		call CountFalsePositives {
			input:
				intervals = intervals,
				ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
				filtered_vcf = Mutect2.filtered_vcf,
				filtered_vcf_index = Mutect2.filtered_vcf_index,
				gatk_override = gatk_override,
				gatk_docker = gatk_docker
		}
	}

	call GatherTables { input: tables = CountFalsePositives.false_positive_counts }

	output {
		File summary = GatherTables.summary
		Array[File] false_positives_vcfs = Mutect2.filtered_vcf
		Array[File] unfiltered_vcfs = Mutect2.unfiltered_vcf
	}
}

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
	File? intervals
	File ref_fasta
	File ref_fai
	File ref_dict
	File filtered_vcf
	File filtered_vcf_index

	File? gatk_override

	String gatk_docker

	command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

	    gatk --java-options "-Xmx4g" CountFalsePositives \
		    -V ${filtered_vcf} \
		    -R ${ref_fasta} \
		    ${"-L " + intervals} \
		    -O false-positives.txt \
	}

    runtime {
        docker: gatk_docker
        memory: "5 GB"
        disks: "local-disk " + 500 + " HDD"
    }

	output {
		File false_positive_counts = "false-positives.txt"
	}
}
