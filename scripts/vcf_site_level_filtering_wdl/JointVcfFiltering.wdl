version 1.0

# This is a workflow for filtering a joint callset VCF using INFO level annotations (so filtering is at the site level).
# Note that the input VCFs here may be sharded by genomic position which may be helpful for large cohorts. The script
# will output the same number of shards that are input.
# This portion of the filtering pipeline will assign a SCORE INFO field annotation to each site, but does not yet apply
# the filtering threshold to the final VCF.

workflow JointVcfFiltering {
	input {
		Array[File] vcf
		Array[File] vcf_index
		File sites_only_vcf
		File sites_only_vcf_index
		String basename

		#TODO delete these inputs
		File python_script
		File hyperparameters_json

		String gatk_docker
		File? extract_interval_list
		File? score_interval_list

		String snp_annotations
		String indel_annotations
		File? gatk_override

		String snp_training_resource_command_line = "--resource:hapmap,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz --resource:omni,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz --resource:1000G,training=true,calibration=false gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
		String indel_training_resource_command_line = "--resource:mills,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
	}

	parameter_meta {
		vcf: "An array of input VCFs that are one callset sharded by genomic region."
		sites_only_vcf: "The full VCF callset without any genotype or sample level information."
		basename: "Desired output file basename."
	}

	call ExtractVariantAnnotations as ExtractVariantAnnotationsSNPs {
		input:
			input_vcf = sites_only_vcf,
			input_vcf_index = sites_only_vcf_index,
			mode = "SNP",
			annotations = snp_annotations,
			resources = snp_training_resource_command_line,
			basename = basename,
			interval_list = extract_interval_list,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call ExtractVariantAnnotations as ExtractVariantAnnotationsINDELs {
		input:
			input_vcf = sites_only_vcf,
			input_vcf_index = sites_only_vcf_index,
			mode = "INDEL",
			annotations = indel_annotations,
			resources = indel_training_resource_command_line,
			basename = basename,
			interval_list = extract_interval_list,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call TrainVariantAnnotationModel as TrainVariantAnnotationModelSNPs {
		input:
			annots = ExtractVariantAnnotationsSNPs.annots,
			basename = basename,
			mode = "snp",
			python_script = python_script,
			hyperparameters_json = hyperparameters_json,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	call TrainVariantAnnotationModel as TrainVariantAnnotationModelINDELs {
		input:
			annots = ExtractVariantAnnotationsINDELs.annots,
			basename = basename,
			mode = "indel",
			python_script = python_script,
			hyperparameters_json = hyperparameters_json,
			gatk_override = gatk_override,
			gatk_docker = gatk_docker
	}

	scatter(idx in range(length(vcf))) {
		call ScoreVariantAnnotations as ScoreVariantAnnotationsSNPs {
			input:
				vcf = vcf[idx],
				vcf_index = vcf_index[idx],
				basename = basename,
				mode = "SNP",
				scoring_python_script = python_script,
				annotations = snp_annotations,
				extracted_training_vcf = ExtractVariantAnnotationsSNPs.extracted_training_vcf,
				extracted_training_vcf_index = ExtractVariantAnnotationsSNPs.extracted_training_vcf_index,
				interval_list = score_interval_list,
				model_files = TrainVariantAnnotationModelSNPs.outputs,
				resources = "-resource:hapmap,training=false,calibration=true,prior=15 gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz -resource:omni,training=false,calibration=true,prior=12 gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
				gatk_override = gatk_override,
				gatk_docker = gatk_docker
		}

		call ScoreVariantAnnotations as ScoreVariantAnnotationsINDELs {
			input:
				vcf = ScoreVariantAnnotationsSNPs.output_vcf,
				vcf_index = ScoreVariantAnnotationsSNPs.output_vcf_index,
				basename = basename,
				mode = "INDEL",
				scoring_python_script = python_script,
				annotations = indel_annotations,
				extracted_training_vcf = ExtractVariantAnnotationsINDELs.extracted_training_vcf,
				extracted_training_vcf_index = ExtractVariantAnnotationsINDELs.extracted_training_vcf_index,
				interval_list = score_interval_list,
				model_files = TrainVariantAnnotationModelINDELs.outputs,
				resources = "--resource:mills,training=false,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
				gatk_override = gatk_override,
				gatk_docker = gatk_docker
		}
	}

	output {
		Array[File] variant_filtered_vcf = ScoreVariantAnnotationsINDELs.output_vcf
		Array[File] variant_filtered_vcf_index = ScoreVariantAnnotationsINDELs.output_vcf_index
	}

}

task ExtractVariantAnnotations {
	input {
		String gatk_docker
		File? gatk_override
		File input_vcf
		File input_vcf_index
		String basename
		String mode
		String annotations
		String resources
		File? interval_list

		Int memory_mb = 14000
		Int command_mem = memory_mb - 1000
	}
	Int disk_size = ceil(size(input_vcf, "GB") + 50)
	command {
		set -e
		export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

		gatk --java-options "-Xmx${command_mem}m" \
			ExtractVariantAnnotations \
			-V ~{input_vcf} \
			-O ~{basename}.~{mode} \
			~{annotations} \
			~{"-L " + interval_list} \
			-mode ~{mode} \
			~{resources}
	}
	output {
		File annots = "~{basename}.~{mode}.annot.hdf5"
		File extracted_training_vcf = "~{basename}.~{mode}.vcf.gz"
		File extracted_training_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " LOCAL"
		memory: memory_mb + " MiB"
	}
}

task TrainVariantAnnotationModel {
	input {
		String gatk_docker
		File? gatk_override
		File annots
		String basename
		String mode
		File python_script
		File hyperparameters_json

		Int memory_mb = 14000
		Int command_mem = memory_mb - 1000
	}
	Int disk_size = ceil(size(annots, "GB") + 100)
	command {
		set -e

		conda install -y --name gatk dill

		export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

		mode=$(echo "~{mode}" | awk '{print toupper($0)}')

		gatk --java-options "-Xmx${command_mem}m" \
			TrainVariantAnnotationsModel \
			--annotations-hdf5 ~{annots} \
			-O ~{basename} \
			--python-script ~{python_script} \
			--hyperparameters-json ~{hyperparameters_json} \
			-mode $mode

	}
	output {
		Array[File] outputs = glob("~{basename}.~{mode}.*")
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " LOCAL"
		memory: memory_mb + " MiB"
	}
}

task ScoreVariantAnnotations {
	input {
		String gatk_docker
		File? gatk_override
		File vcf
		File vcf_index
		String basename
		String mode
		File scoring_python_script
		String annotations
		String resources
		File extracted_training_vcf
		File extracted_training_vcf_index
		File? interval_list
		Array[File] model_files

		Int memory_mb = 16000
		Int command_mem = memory_mb - 1000
	}
	Int disk_size = ceil(size(vcf, "GB") *2 + 50)

	command {
		set -e

		ln -s ~{sep=" . && ln -s " model_files} .

		conda install -y --name gatk dill

		export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

		gatk --java-options "-Xmx${command_mem}m" \
			ScoreVariantAnnotations \
			~{"-L " + interval_list} \
			-V ~{vcf} \
			-O ~{basename}.~{mode} \
			--python-script ~{scoring_python_script} \
			--model-prefix ~{basename} \
			~{annotations} \
			-mode ~{mode} \
			--resource:extracted,extracted=true ~{extracted_training_vcf} \
			~{resources}
	}
	output {
		File scores = "~{basename}.~{mode}.scores.hdf5"
		File annots = "~{basename}.~{mode}.annot.hdf5"
		File output_vcf = "~{basename}.~{mode}.vcf.gz"
		File output_vcf_index = "~{basename}.~{mode}.vcf.gz.tbi"
	}
	runtime {
		docker: gatk_docker
		disks: "local-disk " + disk_size + " LOCAL"
		memory: memory_mb + " MiB"
	}
}

