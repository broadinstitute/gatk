#
# This simple, *unsupported* WDL takes in a VCF from M2 and a tumor bam file.
#  It produces a new VCF with the filtering results.

workflow test_orientation_bias_filter {
    # Format: tsv
    # entity_id vcf   tumor_bam_file
    #
    #  For example:
    # SAMPLE1<tab>SAMPLE1.m2.vcf<tab>SAMPLE1.bam
    # SAMPLE2<tab>SAMPLE2.m2.vcf<tab>SAMPLE2.bam
    # ....
    File input_table
    Array[Array[String]] m2_vcfs = read_tsv(input_table)
    File db_snp
    String gatk_jar
    File ref_fasta

    scatter (row in m2_vcfs) {
        call CollectSequencingArtifactMetrics {
            input:
                entity_id=row[0],
                bam_file=row[2],
                gatk_jar=gatk_jar,
                ref_fasta=ref_fasta,
                output_location_prepend=row[0]
        }
        call FilterByOrientationBias {
            input:
                entity_id=row[0],
                gatk_jar=gatk_jar,
                m2_vcf=row[1],
                pre_adapter_detail_metrics=CollectSequencingArtifactMetrics.pre_adapter_detail_metrics
        }
    }
}

task CollectSequencingArtifactMetrics {
    String entity_id
    File bam_file
    String output_location_prepend
    String gatk_jar
    # /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta
    File ref_fasta

    command {
        java -jar ${gatk_jar} CollectSequencingArtifactMetrics -I ${bam_file} -O ${output_location_prepend}  -R ${ref_fasta} --VALIDATION_STRINGENCY SILENT
    }

    output {
        File pre_adapter_detail_metrics = "${output_location_prepend}.pre_adapter_detail_metrics"
        File pre_adapter_summary_metrics = "${output_location_prepend}.pre_adapter_summary_metrics"
        File bait_bias_detail_metrics = "${output_location_prepend}.bait_bias_detail_metrics"
        File bait_bias_summary_metrics = "${output_location_prepend}.bait_bias_summary_metrics"
    }
}

task FilterByOrientationBias {
    String entity_id
    String gatk_jar
    File m2_vcf
    File pre_adapter_detail_metrics

    command {
        java -jar ${gatk_jar} FilterByOrientationBias -V ${m2_vcf} -P ${pre_adapter_detail_metrics} --output ${entity_id}.ob_filtered.vcf
    }

    output {
        File orientation_bias_vcf = "${entity_id}.ob_filtered.vcf"
        File orientation_bias_vcf_summary = "${entity_id}.ob_filtered.vcf.summary"
    }
}

task MakeSummaryFileList {
    Array[File] files
    String output_file

    command {
            python <<CODE
	    	   import shutil
		   import os
	    	   str_file_list = """${sep="\n" files}"""
 	    	   fp = fopen(${output_file}, 'w')
 	    	   fp.write(str_file_list)
 	    	   fp.close()
            CODE
    }

    output {
        File summary_file_list = "${output_file}"
    }
}