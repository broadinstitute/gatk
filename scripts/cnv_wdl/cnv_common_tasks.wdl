# Tasks common to both the CNV somatic panel and case workflows.
#
#############

# Pad targets in the target file by the specified amount (this was found to improve sensitivity and specificity)
task PadTargets {
    File targets
    Int? padding
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Determine output filename
    String filename = select_first([targets, ""])
    String base_filename = basename(filename, ".tsv")

    command {
        java -Xmx${default="1" mem}g -jar ${gatk_jar} PadTargets \
            --targets ${targets} \
            --padding ${default="250" padding} \
            --output ${base_filename}.padded.tsv
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 2]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File padded_targets = "${base_filename}.padded.tsv"
    }
}

task PreprocessIntervals {
    File? intervals
    File ref_fasta_dict
    Int? padding
    Int? bin_length
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Determine output filename
    String filename = select_first([intervals, "wgs"])
    String base_filename = basename(filename, ".interval_list")

    command {
        java -Xmx${default="2" mem}g -jar ${gatk_jar} PreprocessIntervals \
            ${"-L " + intervals} \
            -sequenceDictionary ${ref_fasta_dict} \
            --padding ${default="250" padding} \
            --binLength ${default="1000" bin_length} \
            --interval_merging_rule OVERLAPPING_ONLY \
            --output ${base_filename}.preprocessed.interval_list
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 2]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File preprocessed_intervals = "${base_filename}.preprocessed.interval_list"
    }
}

# Create a target file with GC annotations
task AnnotateTargets {
    String entity_id
    File intervals
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default="4" mem}g -jar ${gatk_jar} AnnotateTargets \
            --targets ${intervals} \
            --reference ${ref_fasta} \
            --interval_merging_rule OVERLAPPING_ONLY \
            --output ${entity_id}.annotated.tsv
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File annotated_intervals = "${entity_id}.annotated.tsv"
    }
}

task AnnotateIntervals {
    File intervals
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default="4" mem}g -jar ${gatk_jar} AnnotateIntervals \
            -L ${intervals} \
            --reference ${ref_fasta} \
            --interval_merging_rule OVERLAPPING_ONLY \
            --output annotated_intervals.tsv
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File annotated_intervals = "annotated_intervals.tsv"
    }
}

# Collect read counts for germline workflow (TSV output in target format)
task CollectReadCounts {
    File? padded_targets
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Int? wgs_bin_length
    Boolean? keep_non_autosomes
    Boolean? disable_all_read_filters
    Boolean? disable_sequence_dictionary_validation
    Boolean? keep_duplicate_reads
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # If no padded target file is input, then do WGS workflow
    Boolean is_wgs = !defined(padded_targets)

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")

    String read_counts_tsv_filename = "${base_filename}.readCounts.tsv"
    String read_counts_hdf5_filename = if is_wgs then "${base_filename}.readCounts.hdf5" else ""
    String intervals_filename = if is_wgs then "${base_filename}.readCounts.intervals.tsv" else select_first([padded_targets, ""])

    command <<<
        if [ ${is_wgs} = true ]
            then
                java -Xmx${default="8" mem}g -jar ${gatk_jar} SparkGenomeReadCounts \
                    --input ${bam} \
                    --reference ${ref_fasta} \
                    --binLength ${default="1000" wgs_bin_length} \
                    --keepXYMT ${default="false" keep_non_autosomes} \
                    --disableToolDefaultReadFilters ${default="false" disable_all_read_filters} \
                    --disableSequenceDictionaryValidation ${default="true" disable_sequence_dictionary_validation} \
                    $(if [ ${default="true" keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                    --output ${read_counts_tsv_filename} \
                    --writeHdf5
            else
                java -Xmx${default="4" mem}g -jar ${gatk_jar} CalculateTargetCoverage \
                    --input ${bam} \
                    --reference ${ref_fasta} \
                    --targets ${padded_targets} \
                    --groupBy SAMPLE \
                    --transform RAW \
                    --targetInformationColumns FULL \
                    --interval_set_rule UNION \
                    --interval_merging_rule OVERLAPPING_ONLY \
                    --interval_padding 0 \
                    --secondsBetweenProgressUpdates 10.0 \
                    --disableToolDefaultReadFilters ${default="false" disable_all_read_filters} \
                    --disableSequenceDictionaryValidation ${default="true" disable_sequence_dictionary_validation} \
                    $(if [ ${default="true" keep_duplicate_reads} = true ]; then echo " --disableReadFilter NotDuplicateReadFilter "; else echo ""; fi) \
                    --output ${read_counts_tsv_filename}
        fi
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        String entity_id = base_filename
        File read_counts = read_counts_tsv_filename
        File read_counts_hdf5 = read_counts_hdf5_filename   #"" if is_wgs = false
        File intervals = intervals_filename                 #padded_targets if is_wgs = false
    }
}

# Collect counts for ModelSegments workflow
task CollectCounts {
    File intervals
    File bam
    File bam_idx
    String? output_format
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")
    String counts_filename = if !defined(output_format) then "${base_filename}.counts.hdf5" else "${base_filename}.counts.tsv"

    command {
        java -Xmx${default="8" mem}g -jar ${gatk_jar} CollectFragmentCounts \
            --input ${bam} \
            -L ${intervals} \
            --outputFormat ${default="HDF5" output_format} \
            --interval_merging_rule OVERLAPPING_ONLY \
            --output ${counts_filename}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 8]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        String entity_id = base_filename
        File counts = counts_filename
    }
}

# Collect allelic counts
task CollectAllelicCounts {
    File common_sites
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Int? minimum_base_quality
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")

    String allelic_counts_filename = "${base_filename}.allelicCounts.tsv"

    command {
        java -Xmx${default="8" mem}g -jar ${gatk_jar} CollectAllelicCounts \
            -L ${common_sites} \
            --input ${bam} \
            --reference ${ref_fasta} \
            --minimumBaseQuality ${default="20" minimum_base_quality} \
            --output ${allelic_counts_filename}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        String entity_id = base_filename
        File allelic_counts = allelic_counts_filename
    }
}

# Correct coverage profile(s) for sample-specific GC bias
task CorrectGCBias {
    String entity_id
    File coverage   # This can be either single-sample or multi-sample
    File annotated_intervals
    String gatk_jar

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command {
        java -Xmx${default=4 mem}g -jar ${gatk_jar} CorrectGCBias \
          --input ${coverage} \
          --targets ${annotated_intervals} \
          --output ${entity_id}.gc_corrected.tsv
    }

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(coverage, "GB"))+50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File corrected_coverage = "${entity_id}.gc_corrected.tsv"
    }
}