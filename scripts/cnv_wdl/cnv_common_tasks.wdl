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
    File? mappability_track
    File? segmental_duplication_track
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

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx${command_mem_mb}m" AnnotateIntervals \
            -L ${intervals} \
            --reference ${ref_fasta} \
            ${"--mappability-track " + mappability_track} \
            ${"--segmental-duplication-track " + segmental_duplication_track} \
            --feature-query-lookahead ${default=1000000 feature_query_lookahead} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output annotated_intervals.tsv
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File annotated_intervals = "annotated_intervals.tsv"
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

    # Runtime parameters
    String gatk_docker
    Int? mem_gb
    Int? disk_space_gb
    Boolean use_ssd = false
    Int? cpu
    Int? preemptible_attempts

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000

    String base_filename = basename(interval_list, ".interval_list")

    command <<<
        set -e

        grep @ ${interval_list} > header.txt
        grep -v @ ${interval_list} > all_intervals.txt
        split -l ${num_intervals_per_scatter} --numeric-suffixes all_intervals.txt ${base_filename}.scattered.
        for i in ${base_filename}.scattered.*; do cat header.txt $i > $i.interval_list; done
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] scattered_interval_lists = glob("${base_filename}.scattered.*.interval_list")
    }
}

task PostprocessGermlineCNVCalls {
    String entity_id
    Array[File] gcnv_calls_tars
    Array[File] gcnv_model_tars
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
    Boolean allosomal_contigs_specified = defined(allosomal_contigs) && length(select_first([allosomal_contigs, []])) > 0

    String dollar = "$" #WDL workaround for using array[@], see https://github.com/broadinstitute/cromwell/issues/1819

    command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        gcnv_calls_tar_array=(${sep=" " gcnv_calls_tars})
        calls_args=""
        for index in ${dollar}{!gcnv_calls_tar_array[@]}; do
            gcnv_calls_tar=${dollar}{gcnv_calls_tar_array[$index]}
            mkdir CALLS_$index
            tar xzf $gcnv_calls_tar -C CALLS_$index
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

        mkdir extracted-contig-ploidy-calls
        tar xzf ${contig_ploidy_calls_tar} -C extracted-contig-ploidy-calls

        allosomal_contigs_args="--allosomal-contig ${sep=" --allosomal-contig " allosomal_contigs}"

        gatk --java-options "-Xmx${command_mem_mb}m" PostprocessGermlineCNVCalls \
            $calls_args \
            $model_args \
            ${true="$allosomal_contigs_args" false="" allosomal_contigs_specified} \
            --autosomal-ref-copy-number ${ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls extracted-contig-ploidy-calls \
            --sample-index ${sample_index} \
            --output-genotyped-intervals ${genotyped_intervals_vcf_filename} \
            --output-genotyped-segments ${genotyped_segments_vcf_filename}
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