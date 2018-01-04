task PreprocessIntervals {
    File? intervals
    File ref_fasta_dict
    Int? padding
    Int? bin_length
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Determine output filename
    String filename = select_first([intervals, "wgs"])
    String base_filename = basename(filename, ".interval_list")

    command <<<
        set -e
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${default="2" mem}g -jar $GATK_JAR PreprocessIntervals \
            ${"-L " + intervals} \
            --sequence-dictionary ${ref_fasta_dict} \
            --padding ${default="250" padding} \
            --binLength ${default="1000" bin_length} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${base_filename}.preprocessed.interval_list
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 2]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + " HDD"
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
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    command <<<
        set -e
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${default="4" mem}g -jar $GATK_JAR AnnotateIntervals \
            -L ${intervals} \
            --reference ${ref_fasta} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output annotated_intervals.tsv
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 5]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fasta, "GB")) + 50]) + " HDD"
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
    String? format
    File? gatk4_jar_override

    # Runtime parameters
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")
    String counts_filename = if !defined(format) then "${base_filename}.counts.hdf5" else "${base_filename}.counts.tsv"

    command <<<
        set -e
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${default="8" mem}g -jar $GATK_JAR CollectFragmentCounts \
            --input ${bam} \
            -L ${intervals} \
            --format ${default="HDF5" format} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${counts_filename}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: select_first([mem, 8]) + " GB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
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
    Int? mem
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 13000
    Int command_mem = machine_mem - 1000

    # Sample name is derived from the bam filename
    String base_filename = basename(bam, ".bam")

    String allelic_counts_filename = "${base_filename}.allelicCounts.tsv"

    command <<<
        set -e
        GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

        java -Xmx${command_mem}m -jar $GATK_JAR CollectAllelicCounts \
            -L ${common_sites} \
            --input ${bam} \
            --reference ${ref_fasta} \
            --minimum-base-quality ${default="20" minimum_base_quality} \
            --output ${allelic_counts_filename}
    >>>

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(bam, "GB")) + 50]) + " HDD"
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        String entity_id = base_filename
        File allelic_counts = allelic_counts_filename
    }
}