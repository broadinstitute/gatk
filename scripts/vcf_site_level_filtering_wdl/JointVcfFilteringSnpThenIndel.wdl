version 1.0

# This is a workflow for filtering a joint callset VCF using INFO level annotations (so filtering is at the site level).
# Note that the input VCFs here may be sharded by genomic position which may be helpful for large cohorts. The script
# will output the same number of shards that are input.
# This portion of the filtering pipeline will assign a SCORE INFO field annotation to each site, but does not yet apply
# the filtering threshold to the final VCF.

import "JointVcfFiltering.wdl" as JointVcfFiltering

workflow JointVcfFilteringSnpThenIndel {
	input {
		Array[File] input_vcfs
        Array[File] input_vcf_idxs
		File sites_only_vcf
		File sites_only_vcf_idx
		String output_prefix

        # SNP arguments TODO replace with struct
		Array[String] snp_annotations
        String snp_resource_args = "--resource:hapmap,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz --resource:omni,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz --resource:1000G,training=true gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        String? snp_model_backend
        File? snp_python_script
        File? snp_hyperparameters_json
        String? snp_extract_extra_args
        String? snp_train_extra_args
        String? snp_score_extra_args

        # INDEL arguments TODO replace with struct
        Array[String] indel_annotations
        String indel_resource_args = "--resource:mills,training=true,calibration=true gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        String? indel_model_backend
        File? indel_python_script
        File? indel_hyperparameters_json
        String? indel_extract_extra_args
        String? indel_train_extra_args
        String? indel_score_extra_args

        String gatk_docker
		File? gatk_override
	}

	parameter_meta {
		input_vcfs: "An array of input VCFs that are one callset sharded by genomic region."
		sites_only_vcf: "The full VCF callset without any genotype or sample level information."
		output_prefix: "Desired output file prefix."
	}

	call JointVcfFiltering.JointVcfFiltering as JointVcfFilteringSnp {
		input:
			input_vcfs = input_vcfs,
            input_vcf_idxs = input_vcf_idxs,
            sites_only_vcf = sites_only_vcf,
            sites_only_vcf_idx = sites_only_vcf_idx,
            output_prefix = output_prefix,
            annotations = snp_annotations,
            resource_args = snp_resource_args,
            model_backend = snp_model_backend,
            python_script = snp_python_script,
            hyperparameters_json = snp_hyperparameters_json,
            extract_extra_args = "--mode SNP " + select_first([snp_extract_extra_args, ""]),
            train_extra_args = "--mode SNP " + select_first([snp_train_extra_args, ""]),
            score_extra_args = "--mode SNP " + select_first([snp_score_extra_args, ""]),
            gatk_docker = gatk_docker,
            gatk_override = gatk_override
	}

	call JointVcfFiltering.JointVcfFiltering as JointVcfFilteringIndel {
    		input:
    			input_vcfs = JointVcfFilteringSnp.scored_vcfs,
                input_vcf_idxs = JointVcfFilteringSnp.scored_vcf_idxs,
                sites_only_vcf = sites_only_vcf,
                sites_only_vcf_idx = sites_only_vcf_idx,
                output_prefix = output_prefix,
                annotations = indel_annotations,
                resource_args = indel_resource_args,
                model_backend = indel_model_backend,
                python_script = indel_python_script,
                hyperparameters_json = indel_hyperparameters_json,
                extract_extra_args = "--mode INDEL " + select_first([indel_extract_extra_args, ""]),
                train_extra_args = "--mode INDEL " + select_first([indel_train_extra_args, ""]),
                score_extra_args = "--mode INDEL " + select_first([indel_score_extra_args, ""]),
                gatk_docker = gatk_docker,
                gatk_override = gatk_override
    	}

	output {
		Array[File] scored_vcfs = JointVcfFilteringIndel.scored_vcfs
		Array[File] scored_vcf_idxs = JointVcfFilteringIndel.scored_vcf_idxs
	}
}

