task CNNScoreVariants {
    File input_vcf
    File input_vcf_index
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    String output_prefix
    File? bam_file
    File? bam_file_index
    File? architecture_json
    File? architecture_hd5
    Int? inference_batch_size
    Int? transfer_batch_size
    Int? intra_op_threads
    Int? inter_op_threads
    String? tensor_type

    File interval_list
    File? gatk_override

    # Runtime parameters
    Int? mem_gb
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 6000
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem / 2

command <<<

        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
        CNNScoreVariants \
        ${"-I " + bam_file} \
        ${"--read-index " + bam_file_index} \
        -R ${reference_fasta} \
        -V ${input_vcf} \
        -O ${output_prefix}_cnn_annotated.vcf.gz \
        -L ${interval_list} \
        ${"--architecture " + architecture_json} \
        ${"--tensor-type " + tensor_type} \
        ${"--inference-batch-size " + inference_batch_size} \
        ${"--transfer-batch-size " + transfer_batch_size} \
        ${"--intra-op-threads " + intra_op_threads} \
        ${"--inter-op-threads " + inter_op_threads}

>>>

  runtime {
    docker: "${gatk_docker}"
    memory: machine_mem + " MB"
    disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
    cpu: select_first([cpu, 1])
    zones: "us-central1-b" # Needs to be a zone that guarantees CPUs with AVX see (https://cloud.google.com/compute/docs/regions-zones/)
    bootDiskSizeGb: "16"
  }

  output {
    Array[File] log = glob("gatkStreamingProcessJournal*")
    File cnn_annotated_vcf = "${output_prefix}_cnn_annotated.vcf.gz"
    File cnn_annotated_vcf_index = "${output_prefix}_cnn_annotated.vcf.gz.tbi"
  }
}

task RunHC4 {
    File input_bam
    File input_bam_index
    File reference_fasta
    File reference_dict
    File reference_fasta_index
    String output_prefix
    File interval_list
    String extra_args
    File? gatk_override

    # Runtime parameters
    Int? mem_gb
    String gatk_docker
    Int? preemptible_attempts
    Int disk_space_gb
    Int? cpu

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 8000

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
        HaplotypeCaller \
        -R ${reference_fasta} \
        -I ${input_bam} \
        --read-index ${input_bam_index} \
        -O ${output_prefix}_hc4.vcf.gz \
        -L ${interval_list} \
        -bamout ${output_prefix}_bamout.bam \
        ${extra_args}
    }

    output {
        File bamout = "${output_prefix}_bamout.bam"
        File bamout_index = "${output_prefix}_bamout.bai"
        File raw_vcf = "${output_prefix}_hc4.vcf.gz"
        File raw_vcf_index = "${output_prefix}_hc4.vcf.gz.tbi"
    }
    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        # Note that the space before SSD and HDD should be included.
        disks: "local-disk " + sub(disk_space_gb, "\\..*", "") + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: "16"
    }
}


task FilterVariantTranches {
    File input_vcf
    File input_vcf_index
    Array[File] resources
    Array[File] resources_index
    String output_prefix
    String snp_tranches
    String indel_tranches
    String info_key
    String? extra_args
    File? gatk_override

    # Runtime parameters
    Int? mem_gb
    String gatk_docker
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    String output_vcf = "${output_prefix}_cnn_filtered.vcf.gz"

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
    Int default_disk_space_gb = 200

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

command <<<
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
        FilterVariantTranches \
        -V ${input_vcf} \
        --output ${output_vcf} \
        -resource ${sep=" -resource " resources} \
        -info-key ${info_key} \
        ${snp_tranches} \
        ${indel_tranches} \
        ${extra_args}
>>>

  runtime {
    docker: "${gatk_docker}"
    memory: machine_mem + " MB"
    # Note that the space before HDD and HDD should be included.
    disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
    cpu: select_first([cpu, 1])
    bootDiskSizeGb: "16"
  }

  output {
    File cnn_filtered_vcf = "${output_vcf}"
    File cnn_filtered_vcf_index = "${output_vcf}.tbi"
  }
}

task SplitIntervals {
    # inputs
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    Int scatter_count
    String? split_intervals_extra_args

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
            SplitIntervals \
            -R ${ref_fasta} \
            ${"-L " + intervals} \
            -scatter ${scatter_count} \
            -O ./ \
            ${split_intervals_extra_args}
    }

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: "16"
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task MergeVCFs {
    # inputs
    Array[File] input_vcfs
    String output_prefix

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu

    String output_vcf = "${output_prefix}_cnn_scored.vcf.gz"

    Int default_disk_space_gb = 100
    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx${command_mem}m" MergeVcfs \
        -I ${sep=' -I ' input_vcfs} -O "${output_vcf}"
    }

    runtime {
        docker: "${gatk_docker}"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
        bootDiskSizeGb: "16"
    }

    output {
        File merged_vcf = "${output_vcf}"
        File merged_vcf_index = "${output_vcf}.tbi"
    }
}

task CramToBam {
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File cram_file
    String output_prefix

    # Runtime parameters
    Int? mem_gb
    Int? preemptible_attempts
    Int disk_space_gb
    Int? cpu

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000


    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

command <<<
  ls -ltr ${cram_file} ${reference_fasta} &&
  echo "ls (1): complete" &&
  samtools view -h -T ${reference_fasta} ${cram_file} |
  samtools view -b -o ${output_prefix}.bam - &&
  echo "samtools view: complete" &&
  ls -ltr . &&
  echo "ls (2): complete" &&
  samtools index -b ${output_prefix}.bam &&
  echo "samtools index: complete" &&
  ls -ltr . &&
  echo "ls (3): complete" &&
  mv ${output_prefix}.bam.bai ${output_prefix}.bai &&
  echo "mv: complete" &&
  ls -ltr . &&
  echo "ls (4): complete"
  >>>
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.1.1"
    memory: machine_mem + " MB"

    # Note that the space before SSD and HDD should be included.
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
    cpu: select_first([cpu, 1])
  }

  output {
    File output_bam = "${output_prefix}.bam"
    File output_bam_index = "${output_prefix}.bai"
  }
}

task SamtoolsMergeBAMs {
    Array[File] input_bams
    String output_prefix
    Int disk_space_gb
    command {
        samtools merge ${output_prefix}_bamout.bam ${sep=' ' input_bams}
        samtools index ${output_prefix}_bamout.bam ${output_prefix}_bamout.bai
    }

    output {
        File bamout = "${output_prefix}_bamout.bam"
        File bamout_index = "${output_prefix}_bamout.bai"
    }

  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.1.1"
    memory: "16 GB"
    disks: "local-disk " + disk_space_gb + " HDD"
  }
}
