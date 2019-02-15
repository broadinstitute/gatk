task PreprocessIntervals {
    File? intervals
    File? blacklist_intervals
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Int? padding
    Int? bin_length
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # Determine output filename
    String filename = select_first([intervals, "wgs"])
    String base_filename = basename(filename, ".interval_list")

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" PreprocessIntervals \
            ${"-L " + intervals} \
            ${"-XL " + blacklist_intervals} \
            --sequence-dictionary ${ref_fasta_dict} \
            --reference ${ref_fasta} \
            --padding ${default="250" padding} \
            --bin-length ${default="1000" bin_length} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${base_filename}.preprocessed.interval_list
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File preprocessed_intervals = "${base_filename}.preprocessed.interval_list"
    }
}

task AnnotateIntervals {
    File intervals
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    File? mappability_track_bed
    File? mappability_track_bed_idx
    File? segmental_duplication_track_bed
    File? segmental_duplication_track_bed_idx
    Int? feature_query_lookahead
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500
    
    # Determine output filename
    String filename = select_first([intervals, "wgs.preprocessed"])
    String base_filename = basename(filename, ".interval_list")

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" AnnotateIntervals \
            -L ${intervals} \
            --reference ${ref_fasta} \
            ${"--mappability-track " + mappability_track_bed} \
            ${"--segmental-duplication-track " + segmental_duplication_track_bed} \
            --feature-query-lookahead ${default=1000000 feature_query_lookahead} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${base_filename}.annotated.tsv
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File annotated_intervals = "${base_filename}.annotated.tsv"
    }
}

task FilterIntervals {
    File intervals
    File? blacklist_intervals
    File? annotated_intervals
    Array[File]? read_count_files
    Float? minimum_gc_content
    Float? maximum_gc_content
    Float? minimum_mappability
    Float? maximum_mappability
    Float? minimum_segmental_duplication_content
    Float? maximum_segmental_duplication_content
    Int? low_count_filter_count_threshold
    Float? low_count_filter_percentage_of_samples
    Float? extreme_count_filter_minimum_percentile
    Float? extreme_count_filter_maximum_percentile
    Float? extreme_count_filter_percentage_of_samples
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 500
    
    # Determine output filename
    String filename = select_first([intervals, "wgs.preprocessed"])
    String base_filename = basename(filename, ".interval_list")

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" FilterIntervals \
            -L ${intervals} \
            ${"-XL " + blacklist_intervals} \
            ${"--annotated-intervals " + annotated_intervals} \
            ${if defined(read_count_files) then "--input " else ""} ${sep=" --input " read_count_files} \
            --minimum-gc-content ${default="0.1" minimum_gc_content} \
            --maximum-gc-content ${default="0.9" maximum_gc_content} \
            --minimum-mappability ${default="0.9" minimum_mappability} \
            --maximum-mappability ${default="1.0" maximum_mappability} \
            --minimum-segmental-duplication-content ${default="0.0" minimum_segmental_duplication_content} \
            --maximum-segmental-duplication-content ${default="0.5" maximum_segmental_duplication_content} \
            --low-count-filter-count-threshold ${default="5" low_count_filter_count_threshold} \
            --low-count-filter-percentage-of-samples ${default="90.0" low_count_filter_percentage_of_samples} \
            --extreme-count-filter-minimum-percentile ${default="1.0" extreme_count_filter_minimum_percentile} \
            --extreme-count-filter-maximum-percentile ${default="99.0" extreme_count_filter_maximum_percentile} \
            --extreme-count-filter-percentage-of-samples ${default="90.0" extreme_count_filter_percentage_of_samples} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${base_filename}.filtered.interval_list
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File filtered_intervals = "${base_filename}.filtered.interval_list"
    }
}

task CollectCounts {
    File intervals
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    String? format
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")
    String counts_filename = if !defined(format) then "${base_filename}.counts.hdf5" else "${base_filename}.counts.tsv"

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" CollectReadCounts \
            -L ${intervals} \
            --input ${bam} \
            --reference ${ref_fasta} \
            --format ${default="HDF5" format} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${counts_filename}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        String entity_id = base_filename
        File counts = counts_filename
    }
}

task CollectAllelicCounts {
    File common_sites
    File bam
    File bam_idx
    File ref_fasta
    File ref_fasta_fai
    File ref_fasta_dict
    Int? minimum_base_quality
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 13]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")

    String allelic_counts_filename = "${base_filename}.allelicCounts.tsv"

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" CollectAllelicCounts \
            -L ${common_sites} \
            --input ${bam} \
            --reference ${ref_fasta} \
            --minimum-base-quality ${default="20" minimum_base_quality} \
            --output ${allelic_counts_filename}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        String entity_id = base_filename
        File allelic_counts = allelic_counts_filename
    }
}

task ScatterIntervals {
    File interval_list
    Int num_intervals_per_scatter
    String? output_dir
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = machine_mem_mb - 500

    # If optional output_dir not specified, use "out";
    String output_dir_ = select_first([output_dir, "out"])

    String base_filename = basename(interval_list, ".interval_list")

    command <<<
        set -e
        mkdir ${output_dir_}
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}
        
        {
            >&2 echo "Attempting to run IntervalListTools..."
            gatk --java-options "-Xmx${command_mem_mb}m" IntervalListTools \
                --INPUT ${interval_list} \
                --SUBDIVISION_MODE INTERVAL_COUNT \
                --SCATTER_CONTENT ${num_intervals_per_scatter} \
                --OUTPUT ${output_dir_} &&
            # output files are named output_dir_/temp_0001_of_N/scattered.interval_list, etc. (N = number of scatters);
            # we rename them as output_dir_/base_filename.scattered.0000.interval_list, etc.
            ls ${output_dir_}/*/scattered.interval_list | \
                cat -n | \
                while read n filename; do mv $filename ${output_dir_}/${base_filename}.scattered.$(printf "%04d" $n).interval_list; done
            rm -rf ${output_dir_}/temp_*_of_*
        } || {
            # if only a single shard is required, then we can just rename the original interval list
            >&2 echo "IntervalListTools failed because only a single shard is required. Copying original interval list..."
            cp ${interval_list} ${output_dir_}/${base_filename}.scattered.1.interval_list
        }
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] scattered_interval_lists = glob("${output_dir_}/${base_filename}.scattered.*.interval_list")
    }
}

task PostprocessGermlineCNVCalls {
    String entity_id
    Array[File] gcnv_calls_tars
    Array[File] gcnv_model_tars
    Array[File] calling_configs
    Array[File] denoising_configs
    Array[File] gcnvkernel_version
    Array[File] sharded_interval_lists
    File contig_ploidy_calls_tar
    Array[String]? allosomal_contigs
    Int ref_copy_number_autosomal_contigs
    Int sample_index
    File? gatk4_jar_override

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = machine_mem_mb - 1000

    String genotyped_intervals_vcf_filename = "genotyped-intervals-${entity_id}.vcf.gz"
    String genotyped_segments_vcf_filename = "genotyped-segments-${entity_id}.vcf.gz"

    Array[String] allosomal_contigs_args = if defined(allosomal_contigs) then prefix("--allosomal-contig ", select_first([allosomal_contigs])) else []

    String dollar = "$" #WDL workaround for using array[@], see https://github.com/broadinstitute/cromwell/issues/1819

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        sharded_interval_lists_array=(${sep=" " sharded_interval_lists})

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        # also copy over shard config and interval files
        gcnv_calls_tar_array=(${sep=" " gcnv_calls_tars})
        calling_configs_array=(${sep=" " calling_configs})
        denoising_configs_array=(${sep=" " denoising_configs})
        gcnvkernel_version_array=(${sep=" " gcnvkernel_version})
        sharded_interval_lists_array=(${sep=" " sharded_interval_lists})
        calls_args=""
        for index in ${dollar}{!gcnv_calls_tar_array[@]}; do
            gcnv_calls_tar=${dollar}{gcnv_calls_tar_array[$index]}
            mkdir -p CALLS_$index/SAMPLE_${sample_index}
            tar xzf $gcnv_calls_tar -C CALLS_$index/SAMPLE_${sample_index}
            cp ${dollar}{calling_configs_array[$index]} CALLS_$index/
            cp ${dollar}{denoising_configs_array[$index]} CALLS_$index/
            cp ${dollar}{gcnvkernel_version_array[$index]} CALLS_$index/
            cp ${dollar}{sharded_interval_lists_array[$index]} CALLS_$index/
            calls_args="$calls_args --calls-shard-path CALLS_$index"
        done

        # untar models to MODEL_0, MODEL_1, etc directories and build the command line
        gcnv_model_tar_array=(${sep=" " gcnv_model_tars})
        model_args=""
        for index in ${dollar}{!gcnv_model_tar_array[@]}; do
            gcnv_model_tar=${dollar}{gcnv_model_tar_array[$index]}
            mkdir MODEL_$index
            tar xzf $gcnv_model_tar -C MODEL_$index
            model_args="$model_args --model-shard-path MODEL_$index"
        done

        mkdir contig-ploidy-calls
        tar xzf ${contig_ploidy_calls_tar} -C contig-ploidy-calls

        gatk --java-options "-Xmx${command_mem_mb}m" PostprocessGermlineCNVCalls \
            $calls_args \
            $model_args \
            ${sep=" " allosomal_contigs_args} \
            --autosomal-ref-copy-number ${ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls contig-ploidy-calls \
            --sample-index ${sample_index} \
            --output-genotyped-intervals ${genotyped_intervals_vcf_filename} \
            --output-genotyped-segments ${genotyped_segments_vcf_filename}

        rm -rf CALLS_*
        rm -rf MODEL_*
        rm -rf contig-ploidy-calls
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File genotyped_intervals_vcf = genotyped_intervals_vcf_filename
        File genotyped_segments_vcf = genotyped_segments_vcf_filename
    }
}
