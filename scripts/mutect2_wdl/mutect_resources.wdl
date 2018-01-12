# given a gnomAD vcf, produce:
# 1) a sites-only vcf where the only INFO field is allele frequency (AF)
#   this is used as the gnomad input to mutect2.wdl
# 2) this sites-only vcf further restricted to a given minimum allele frequency
#   this is used as the variants_for_contamination input to mutect2.wdl

workflow MutectResources {
    # inputs
    File gnomad
    File gnomad_idx
    File? intervals
    Float minimum_allele_frequency

    File? gatk_override

    # runtime
    String gatk_docker

    call SelectIntervals {
        input:
            input_vcf = gnomad,
            input_vcf_idx = gnomad_idx,
            intervals = intervals,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    call MakeAlleleFrequencyOnlyVcf {
        input:
            input_vcf = SelectIntervals.output_vcf,
            output_name = "gnomad-for-mutect2",
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    call SelectCommonBiallelicSNPs {
        input:
            input_vcf = MakeAlleleFrequencyOnlyVcf.output_vcf,
            input_vcf_idx = MakeAlleleFrequencyOnlyVcf.output_vcf_idx,
            minimum_allele_frequency = minimum_allele_frequency,
            output_name = "mutect2-contamination-variants",
            gatk_override = gatk_override,
            gatk_docker = gatk_docker
    }

    output {
        File m2_gnomad = MakeAlleleFrequencyOnlyVcf.output_vcf
        File m2_gnomad_idx = MakeAlleleFrequencyOnlyVcf.output_vcf_idx
        File contamination_variants = SelectCommonBiallelicSNPs.output_vcf
        File contamination_variants_idx = SelectCommonBiallelicSNPs.output_vcf_idx
    }
}

task SelectIntervals {
    # inputs
    File input_vcf
    File input_vcf_idx
    File? intervals

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    Int machine_mem = select_first([mem, 3])
    Int command_mem = machine_mem - 1

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx${command_mem}g" SelectVariants -V ${input_vcf} ${"-L " + intervals} -O selected.vcf --lenient
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
    }

    output {
        File output_vcf = "selected.vcf"
        File output_vcf_index = "selected.vcf.idx"
    }
}


# clear ID and QUAL fields and delete all INFO fields other than AF
# note that input must be a plain-text vcf, not a vcf.gz.
# this task re-indexes and compresses the output vcf
task MakeAlleleFrequencyOnlyVcf {
    # inputs
    File input_vcf
    String output_name

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    Int machine_mem = select_first([mem, 3])
    Int command_mem = machine_mem - 1

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        grep '#' ${input_vcf} > header
        grep -v '#' ${input_vcf} | grep -P "\tPASS\t" > body

        # delete any INFO fields before AF (PASS<TAB><other info fields>;AF=___ --> PASS<TAB>AF=____)
        # delete any INFO fields after AF (;<other info fields><end of line> --> nothing)
        sed -e 's/PASS\t.*AF/PASS\tAF/g' -e 's/[;].*$//g' body > simplified_info

        # replace ID (3rd) and QUAL (6th) columns with '.' (empty)
        while read contig pos id ref alt qual filter info; do
          printf "$contig\t$pos\t.\t$ref\t$alt\t.\t$filter\t$info\n"
        done < simplified_info > simplified_body

        cat header simplified_body > simplified.vcf

        gatk --java-options "-Xmx${command_mem}g" IndexFeatureFile -F simplified.vcf
        gatk --java-options "-Xmx${command_mem}g" SelectVariants -V simplified.vcf -O ${output_name}.vcf.gz
        rm header body simplified_info simplified_body simplified.vcf simplified.vcf.idx
    }

    output {
        File output_vcf = "${output_name}.vcf.gz"
        File output_vcf_idx = "${output_name}.vcf.gz.tbi"
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
    }
}

task SelectCommonBiallelicSNPs {
    # inputs
    File input_vcf
    File input_vcf_idx
    Float minimum_allele_frequency
    String output_name

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    Int machine_mem = select_first([mem, 3])
    Int command_mem = machine_mem - 1

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}g" SelectVariants \
            -V ${input_vcf} \
            -select-type SNP -restrict-alleles-to BIALLELIC \
            -select "AF > ${minimum_allele_frequency}" \
            -O ${output_name}.vcf.gz \
            --lenient
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        cpu: select_first([cpu, 1])
    }

    output {
        File output_vcf = "${output_name}.vcf.gz"
        File output_vcf_idx = "${output_name}.vcf.gz.tbi"
    }
}