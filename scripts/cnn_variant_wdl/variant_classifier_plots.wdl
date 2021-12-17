# Run VCF evaluation over input vcf given a validated truth vcf and confidence region
workflow VariantClassifierPlots {
    File call_vcf                   # VCF to be evaluated
    File call_vcf_index             # Index of VCF to be evaluated
    String? call_sample
    String score_key

    File? truth_vcf                 # Optional truth VCF. If provided, plot colors show true positives and
    File? truth_vcf_index           # true negatives in green with false positives in red and false negatives in yellow.
    String? truth_sample            # Otherwise, plot colors show filtered variants in red and passing variant in green.

    File? intervals

    File rscript

    String gatk_docker
    File? gatk_override

    Int? preemptible_attempts
    Int? disk_space
    Int? mem_gb
    Int? cpu

    if(defined(truth_vcf)){
        call MakeTables {
            input:
                call_vcf = call_vcf,
                call_vcf_index = call_vcf_index,
                call_sample = call_sample,
                score_key = score_key,
                truth_vcf = truth_vcf,
                truth_vcf_index = truth_vcf_index,
                truth_sample = truth_sample,
                intervals = intervals,
                gatk_docker = gatk_docker,
                cpu = cpu,
                mem_gb = mem_gb,
                disk_space = disk_space,
                preemptible_attempts = preemptible_attempts
        }

        call MakePlots{
            input:
                rscript = rscript,
                call_table = MakeTables.call_table,
                truth_table = MakeTables.truth_table,
                score_key = score_key,
                cpu = cpu,
                mem_gb = mem_gb,
                disk_space = disk_space,
                preemptible_attempts = preemptible_attempts
        }

        output {
            MakeTables.*
            MakePlots.*
        }
    }

    if(!defined(truth_vcf)){
        call MakeTableNoTruth {
            input:
                call_vcf = call_vcf,
                call_vcf_index = call_vcf_index,
                call_sample = call_sample,
                score_key = score_key,
                gatk_docker = gatk_docker,
                cpu = cpu,
                mem_gb = mem_gb,
                disk_space = disk_space,
                preemptible_attempts = preemptible_attempts
        }

        call MakePlots as MakePlotsNoTruth {
            input:
                rscript = rscript,
                call_table = MakeTableNoTruth.call_table,
                score_key = score_key,
                cpu = cpu,
                mem_gb = mem_gb,
                disk_space = disk_space,
                preemptible_attempts = preemptible_attempts
        }

        output {
            MakeTableNoTruth.*
            MakePlotsNoTruth.*
        }

    }

}

task MakeTables {
    File call_vcf
    File call_vcf_index
    String? call_sample
    String score_key

    File? truth_vcf
    File? truth_vcf_index
    String? truth_sample

    File? intervals

    # Runtime parameters
    String gatk_docker
    File? gatk_override
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu

    String call_table_name = basename(call_vcf) + ".table"
    String sd_fix_vcf = "call_sd_fix.vcf"

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
            UpdateVcfSequenceDictionary \
            --INPUT=${call_vcf} \
            --OUTPUT=${sd_fix_vcf} \
            -SD=${truth_vcf}

        gatk IndexFeatureFile -F ${sd_fix_vcf}

        gatk --java-options "-Xmx${command_mem}m" \
            GenotypeConcordance \
            --CALL_VCF=${sd_fix_vcf} \
            ${"--CALL_SAMPLE=" + call_sample} \
            --TRUTH_VCF=${truth_vcf} \
            ${"--TRUTH_SAMPLE=" + truth_sample} \
            ${"--INTERVALS=" + intervals} \
            --OUTPUT_VCF=true \
            --IGNORE_FILTER_STATUS \
            -O=concordance

        gatk --java-options "-Xmx${command_mem}m" \
            VariantsToTable \
            -V ${sd_fix_vcf} \
            -F CHROM -F POS -F REF -F ALT -F FILTER -F ${score_key} \
            -F EVENTLENGTH -F AC -F MULTI-ALLELIC -F TRANSITION -F TYPE \
            --show-filtered \
            -O ${call_table_name}

        gatk --java-options "-Xmx${command_mem}m" \
            VariantsToTable \
            -V concordance.genotype_concordance.vcf.gz \
            -F CHROM -F POS -F REF -F ALT -F CONC_ST \
            -O truth.table
    }

    output {
        File call_table = "${call_table_name}"
        File truth_table = "truth.table"
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        # Note that the space before SSD and HDD should be included.
        disks: "local-disk " + default_disk_space_gb + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
        zones: "us-east4-a"
        bootDiskSizeGb: "16"
    }
}

task MakeTableNoTruth {
    File call_vcf
    File call_vcf_index
    String? call_sample
    String score_key

    # Runtime parameters
    String gatk_docker
    File? gatk_override
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu

    String call_table_name = basename(call_vcf) + ".table"

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
            VariantsToTable \
            -V ${call_vcf} \
            -F CHROM -F POS -F REF -F ALT -F FILTER -F ${score_key} \
            -F EVENTLENGTH -F AC -F MULTI-ALLELIC -F TRANSITION -F TYPE \
            --show-filtered \
            -O ${call_table_name}
    }

    output {
        File call_table = "${call_table_name}"
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        # Note that the space before SSD and HDD should be included.
        disks: "local-disk " + default_disk_space_gb + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
        zones: "us-east4-a"
        bootDiskSizeGb: "16"
    }
}

task MakePlots {
    File rscript
    File call_table
    File? truth_table
    String score_key

    # Runtime parameters
    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).
    Int default_disk_space_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    command {
        Rscript ${rscript} ${call_table} ${truth_table} ${score_key}
    }

    output {
        Array[File] plots = glob("*png")
    }

    runtime {
        docker: "rocker/tidyverse:4.1"

        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }
}
