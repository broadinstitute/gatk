# Run the hap.py VCF evaluation over input vcfs given a validated truth vcf and confidence region
workflow HappyWorkflow {
    File vcf_fofn           # Text file with list of VCFs to evaluate with hap.py

    File reference_fasta
    File reference_dict
    File reference_fasta_index

    File truth_vcf
    File truth_vcf_index
    File truth_bed

    File rscript

    Int? preemptible_attempts
    Int? disk_space
    Int? mem_gb
    Int? cpu

    call RunHappy {
        input:
            vcf_fofn = vcf_fofn,
            truth_vcf = truth_vcf,
            truth_vcf_index = truth_vcf_index,
            truth_bed = truth_bed,
            reference_fasta = reference_fasta,
            reference_dict = reference_dict,
            reference_fasta_index = reference_fasta_index,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_space = disk_space,
            preemptible_attempts = preemptible_attempts
    }

    call RunHappyPlots{
        input:
            happy_outputs = RunHappy.happy_outputs,
            rscript = rscript,
            cpu = cpu,
            mem_gb = mem_gb,
            disk_space = disk_space,
            preemptible_attempts = preemptible_attempts
    }

    output {
        RunHappy.*
        RunHappyPlots.*
    }
}

task RunHappy {
    File vcf_fofn
    Array[File] vcf_files=read_lines(vcf_fofn)

    File reference_fasta
    File reference_dict
    File reference_fasta_index

    File truth_vcf
    File truth_vcf_index
    File truth_bed

    # Runtime parameters
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    command {
        for vcf_file in ${sep=" " vcf_files}; do
            vname=$(basename "$vcf_file")
            /opt/hap.py/bin/hap.py \
            ${truth_vcf} \
            "$vcf_file" \
            -f ${truth_bed} \
            -r ${reference_fasta} \
            -o ./happy_"$vname"
        done
    }

    output {
        Array[File] happy_outputs = glob("./happy_*")
    }

    runtime {
        docker: "pkrusche/hap.py"

        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }
}

task RunHappyPlots {
    Array[File] happy_outputs
    File rscript

    # Runtime parameters
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    command {
        for file in ${sep=" " happy_outputs}; do
            mv "$file" ./
        done
        find `pwd`

        Rscript ${rscript}
    }

    output {
        Array[File] plots = glob("*png")
    }

    runtime {
        docker: "rocker/tidyverse"

        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }
}


