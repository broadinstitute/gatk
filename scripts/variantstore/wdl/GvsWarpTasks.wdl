version 1.0

task SNPsVariantRecalibratorCreateModel {

    input {
        String recalibration_filename
        String tranches_filename
        Int downsampleFactor
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

        Int disk_size
        File? gatk_override
    }

    command <<<
        set -euo pipefail
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        gatk --java-options -Xms100g \
        VariantRecalibrator \
        -V ~{sites_only_variant_filtered_vcf} \
        -O ~{recalibration_filename} \
        --tranches-file ~{tranches_filename} \
        --trust-all-polymorphic \
        -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
        -an ~{sep=' -an ' recalibration_annotation_values} \
        ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
        -mode SNP \
        --sample-every-Nth-variant ~{downsampleFactor} \
        --output-model ~{model_report_filename} \
        --max-gaussians ~{max_gaussians} \
        -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
        -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
        -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
    >>>

    runtime {
        memory: "104 GiB"
        cpu: "2"
        bootDiskSizeGb: 15
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots:varstore_d8a72b825eab2d979c8877448c0ca948fd9b34c7_change_to_hwe"
    }

    output {
        File model_report = "~{model_report_filename}"
    }
}

task GatherTranches {

    input {
        Array[File] tranches
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
        --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
        GatherTranches \
        --input inputs.list \
        --mode ~{mode} \
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
