# Tasks common to both the CNV somatic panel and case workflows.
#
#############

# Pad targets in the target file by the specified amount (this was found to improve sensitivity and specificity)
task PadTargets {
    File? targets
    Int? padding
    File gatk_jar
    Int? mem

    # Determine output filename
    String filename = select_first([targets, ""])
    String base_filename = sub(sub(sub(filename, "gs://", ""), "[/]*.*/", ""), "\\.tsv$", "")

    command {
        echo ${filename}; \
        java -Xmx${default=1 mem}g -jar ${gatk_jar} PadTargets \
            --targets ${targets} \
            --padding ${default=250 padding} \
            --output ${base_filename}.padded.tsv
    }

    output {
        File padded_targets = "${base_filename}.padded.tsv"
    }
}

# Collect proportional coverage
task CollectCoverage {
    File? padded_targets
    File bam
    File bam_idx
    String? transform
    Boolean? keep_non_autosomes
    Boolean? disable_all_read_filters
    Boolean? disable_sequence_dictionary_validation
    Boolean? keep_duplicate_reads
    Int? wgs_bin_size
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    Int? mem

    # If no padded target file is input, then do WGS workflow
    Boolean is_wgs = !defined(padded_targets)

    # Sample name is derived from the bam filename
    String base_filename = sub(sub(sub(bam, "gs://", ""), "[/]*.*/", ""), "\\.bam$", "")
 
    # Output file name depending on type of coverage
    String cov_output_name = if (is_wgs && (select_first([transform, ""])  == "RAW")) then "${base_filename}.coverage.tsv.raw_cov" else "${base_filename}.coverage.tsv"
  
    command <<<
        if [ ${is_wgs} = true ]
            then
                java -Xmx${default=4 mem}g -jar ${gatk_jar} SparkGenomeReadCounts \
                    --input ${bam} \
                    --reference ${ref_fasta} \
                    --binsize ${default=10000 wgs_bin_size} \
                    --keepXYMT ${default="false" keep_non_autosomes} \
                    --disableToolDefaultReadFilters ${default="false" disable_all_read_filters} \
                    --disableSequenceDictionaryValidation ${default="true" disable_sequence_dictionary_validation} \
                    $(if [ ${default="true" keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                    --outputFile ${base_filename}.coverage.tsv
            else
                java -Xmx${default=4 mem}g -jar ${gatk_jar} CalculateTargetCoverage \
                    --input ${bam} \
                    --reference ${ref_fasta} \
                    --targets ${padded_targets} \
                    --groupBy SAMPLE \
                    --transform ${default="PCOV" transform} \
                    --targetInformationColumns FULL \
                    --interval_set_rule UNION \
                    --interval_padding 0 \
                    --secondsBetweenProgressUpdates 10.0 \
                    --disableToolDefaultReadFilters ${default="false" disable_all_read_filters} \
                    --disableSequenceDictionaryValidation ${default="true" disable_sequence_dictionary_validation} \
                    $(if [ ${default="true" keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                    --output ${base_filename}.coverage.tsv
        fi
    >>>

    output {
        String entity_id = base_filename
        File coverage = cov_output_name 
    }
}

# Create a target file with GC annotations
task AnnotateTargets {
    String entity_id
    File targets
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File gatk_jar
    Int? mem

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} AnnotateTargets \
            --targets ${targets} \
            --reference ${ref_fasta} \
            --output ${entity_id}.annotated.tsv
    }

    output {
        File annotated_targets = "${entity_id}.annotated.tsv"
    }
}

# Correct coverage profile(s) for sample-specific GC bias
task CorrectGCBias {
    String entity_id
    File coverage   # This can be either single-sample or multi-sample
    File annotated_targets
    File gatk_jar
    Int? mem

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} CorrectGCBias \
          --input ${coverage} \
          --targets ${annotated_targets} \
          --output ${entity_id}.gc_corrected.tsv
    }

    output {
        File corrected_coverage = "${entity_id}.gc_corrected.tsv"
    }
}
