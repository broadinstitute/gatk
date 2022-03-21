version 1.0

task SNPsVariantRecalibratorCreateModel {

    input {
        String recalibration_filename
        String tranches_filename
        String model_report_filename

        Array[String] recalibration_tranche_values
        Array[String] recalibration_annotation_values

        File sites_only_variant_filtered_vcf
        File sites_only_variant_filtered_vcf_index

        File hapmap_resource_vcf
        File omni_resource_vcf
        File one_thousand_genomes_resource_vcf
        File dbsnp_resource_vcf
        File hapmap_resource_vcf_index
        File omni_resource_vcf_index
        File one_thousand_genomes_resource_vcf_index
        File dbsnp_resource_vcf_index
        Boolean use_allele_specific_annotations
        Int max_gaussians = 6
        Int? machine_mem_gb

        Int disk_size
    }

    Int machine_mem = select_first([machine_mem_gb, 100])
    Int java_mem = machine_mem - 5

    command <<<
        set -euo pipefail

        gatk --java-options -Xms~{java_mem}g \
        VariantRecalibrator \
        -V ~{sites_only_variant_filtered_vcf} \
        -O ~{recalibration_filename} \
        --tranches-file ~{tranches_filename} \
        --trust-all-polymorphic \
        -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
        -an ~{sep=' -an ' recalibration_annotation_values} \
        ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
        -mode SNP \
        --sample-every-Nth-variant 1 \
        --output-model ~{model_report_filename} \
        --max-gaussians ~{max_gaussians} \
        -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
        -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
        -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
    >>>

    runtime {
        memory: "~{machine_mem} GiB"
        cpu: "2"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    }

    output {
        File model_report = "~{model_report_filename}"
    }
}

task GatherTranches {

    input {
        Array[File] tranches
        Array[String] output_tranche_values
        String output_filename
        String mode
        Int disk_size
        File? gatk_override
    }

    parameter_meta {
        tranches: {
            localization_optional: true
        }
    }

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        tranches_fofn=~{write_lines(tranches)}

        # Jose says:
        # Cromwell will fall over if we have it try to localize tens of thousands of files,
        # so we manually localize files using gsutil.
        # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
        # PAPI doesn't do.

        # This is here to deal with the JES bug where commands may be run twice
        rm -rf tranches
        mkdir tranches
        RETRY_LIMIT=5

        count=0
        until cat $tranches_fofn | gsutil -m cp -L cp.log -c -I tranches/; do
        sleep 1
        ((count++)) && ((count >= $RETRY_LIMIT)) && break
        done
        if [ "$count" -ge "$RETRY_LIMIT" ]; then
        echo 'Could not copy all the tranches from the cloud' && exit 1
        fi

        cat $tranches_fofn | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list

        gatk --java-options -Xms6g \
        GatherTranches \
        --input inputs.list \
        --mode ~{mode} \
        -tranche ~{sep=' -tranche ' output_tranche_values} \
        --output ~{output_filename}
    >>>

    runtime {
        memory: "7.5 GiB"
        cpu: "2"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    }

    output {
        File tranches_file = "~{output_filename}"
    }
}

task IndelsVariantRecalibrator {

    input {
        String recalibration_filename
        String tranches_filename
        File? model_report

        Array[String] recalibration_tranche_values
        Array[String] recalibration_annotation_values

        File sites_only_variant_filtered_vcf
        File sites_only_variant_filtered_vcf_index

        File mills_resource_vcf
        File axiomPoly_resource_vcf
        File dbsnp_resource_vcf
        File mills_resource_vcf_index
        File axiomPoly_resource_vcf_index
        File dbsnp_resource_vcf_index
        Boolean use_allele_specific_annotations
        Int max_gaussians = 4

        Int disk_size
        Int? machine_mem_gb
    }

    Int machine_mem = select_first([machine_mem_gb, 30])
    Int java_mem = machine_mem - 5


    command <<<
        set -euo pipefail

        gatk --java-options -Xmx~{java_mem}g \
        VariantRecalibrator \
        -V ~{sites_only_variant_filtered_vcf} \
        -O ~{recalibration_filename} \
        --output-model indels.model \
        --rscript-file indels.Rscript \
        --tranches-file ~{tranches_filename} \
        --trust-all-polymorphic \
        -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
        -an ~{sep=' -an ' recalibration_annotation_values} \
        ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
        -mode INDEL \
        ~{"--input-model " + model_report} \
        --max-gaussians ~{max_gaussians} \
        -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
        -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}
    >>>

    runtime {
        memory: "~{machine_mem} GiB"
        cpu: "2"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    }

    output {
        File recalibration = "~{recalibration_filename}"
        File recalibration_index = "~{recalibration_filename}.idx"
        File tranches = "~{tranches_filename}"
        File model = "indels.model"
        File rscript = "indels.Rscript"
    }
}

task SNPsVariantRecalibrator {

    input {
        String recalibration_filename
        String tranches_filename
        File? model_report

        Array[String] recalibration_tranche_values
        Array[String] recalibration_annotation_values

        File sites_only_variant_filtered_vcf
        File sites_only_variant_filtered_vcf_index

        File hapmap_resource_vcf
        File omni_resource_vcf
        File one_thousand_genomes_resource_vcf
        File dbsnp_resource_vcf
        File hapmap_resource_vcf_index
        File omni_resource_vcf_index
        File one_thousand_genomes_resource_vcf_index
        File dbsnp_resource_vcf_index
        Boolean use_allele_specific_annotations
        Int max_gaussians = 6

        Int disk_size
        Int? machine_mem_gb

    }

    Int auto_mem = ceil(2 * size([sites_only_variant_filtered_vcf,
                                 hapmap_resource_vcf,
                                 omni_resource_vcf,
                                 one_thousand_genomes_resource_vcf,
                                 dbsnp_resource_vcf],
                            "GiB"))
    Int machine_mem = select_first([machine_mem_gb, if auto_mem < 7 then 7 else auto_mem])
    Int java_mem = machine_mem - 5
    String model_report_arg = if defined(model_report) then "--input-model $MODEL_REPORT --output-tranches-for-scatter" else ""

    command <<<
        set -euo pipefail

        MODEL_REPORT=~{model_report}

        gatk --java-options -Xmx~{java_mem}g \
        VariantRecalibrator \
        -V ~{sites_only_variant_filtered_vcf} \
        -O ~{recalibration_filename} \
        --output-model snps.model \
        --rscript-file snps.Rscript \
        --tranches-file ~{tranches_filename} \
        --trust-all-polymorphic \
        -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
        -an ~{sep=' -an ' recalibration_annotation_values} \
        ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
        -mode SNP \
        --sample-every-Nth-variant 1 \
        ~{model_report_arg} \
        --max-gaussians ~{max_gaussians} \
        -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
        -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
        -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
    >>>

    runtime {
        memory: "~{machine_mem} GiB"
        cpu: 2
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
        docker: "us.gcr.io/broad-gatk/gatk:4.1.9.0"
    }

    output {
        File recalibration = "~{recalibration_filename}"
        File recalibration_index = "~{recalibration_filename}.idx"
        File tranches = "~{tranches_filename}"
        File model = "snps.model"
        File rscript = "snps.Rscript"
    }
}
