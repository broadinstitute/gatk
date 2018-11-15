# Evaluate site-level concordance of an input VCF against a truth VCF.
#
# Evaluates two variant callsets against each other and produces a six-column summary metrics table.
# The summary:
#
#  - stratifies SNP and INDEL calls
#  - tallies true-positive, false-positive and false-negative calls
#  - calculates sensitivity and precision
#
# Description of inputs:
#
#   Required:
#     gatk_docker                    - GATK Docker image in which to run
#     actual_vcf                     - The produced data set.  Variant Context File (VCF) containing the variants.
#     actual_vcf_idx                 - The produced data set index.  Variant Context File (VCF) index for the given actual_vcf.
#     expected_vcf                   - The truth data set.  Variant Context File (VCF) containing the variants.
#     expected_vcf_idx               - The truth data set index.  Variant Context File (VCF) index for the given expected_vcf.
#
#   Optional:
#     intervals                      - File containing intervals over which the comparison should be run.
#     masks                          - File containing intervals to exclude from the comparison.
#     gatk4_jar_override             - Override Jar file containing GATK 4.0.  Use this when overriding the docker JAR or when using a backend without docker.
#     mem                            - Amount of memory to give the runtime environment.
#     disk_space_gb                  - Amount of disk space to give the runtime environment.
#     cpu                            - Number of CPUs to give the runtime environment.
#     boot_disk_size_gb              - Amount of boot disk space to give the runtime environment.
#     use_ssd                        - Whether or not to mandate the use of Solid State Drives in the runtime environment.
#     preemptible_attempts           - Number of times the comparison can be preempted.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
# Adapted from: https://portal.firecloud.org/#methods/davidben/m2-concordance/10/wdl
workflow VCFVariantCallerConcordance {
    ####################################################################################
    # Inputs:
    File actual_vcf
    File actual_vcf_idx
    File expected_vcf
    File expected_vcf_idx

    File? intervals
    File? masks

    ####################################################################################
    # Runtime Inputs:
    String gatk_docker

    File? gatk4_jar_override
    Int? mem
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible_attempts

    ####################################################################################
    # Call sub-workflows:
    call Concordance {
        input:
            eval_vcf = actual_vcf,
            eval_vcf_idx = actual_vcf_idx,
            truth_vcf = expected_vcf,
            truth_vcf_idx = expected_vcf_idx,
            intervals = intervals,
            masks = masks,

            gatk_docker = gatk_docker,
            gatk_override = gatk4_jar_override,
            mem = mem,
            disk_space_gb = disk_space_gb,
            cpu = cpu,
            boot_disk_size_gb = boot_disk_size_gb,
            use_ssd = use_ssd,
            preemptible_attempts = preemptible_attempts
    }

    ####################################################################################
    # Outputs:
    output {
        File fn = Concordance.fn
        File fn_idx = Concordance.fn_idx
        File fp = Concordance.fp
        File fp_idx = Concordance.fp_idx
        File tp = Concordance.tp
        File tp_idx = Concordance.tp_idx
        File ffn = Concordance.ffn
        File ffn_idx = Concordance.ffn_idx
        File summary = Concordance.summary
        File filter_analysis = Concordance.filter_analysis
    }
}

#=======================================================================================

task Concordance {
    ####################################################################################
    # Inputs:
    File? intervals
    File? masks
    File truth_vcf
    File truth_vcf_idx
    File eval_vcf
    File eval_vcf_idx

    ####################################################################################
    # Runtime Inputs:
    String gatk_docker
    File? gatk_override
    Int? mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb
    Boolean? use_ssd

    ####################################################################################
    # Define default values and set up values for running:
    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3000
    # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
    Int default_disk_space_gb = 100

    Int default_boot_disk_size_gb = 15

    Int default_min_preemptable_attempts = 2

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem *1000 else default_ram_mb
    Int command_mem = machine_mem - 1000

    ####################################################################################
    # Do the work:
    command {
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" \
            Concordance \
                ${"-L " + intervals} \
                ${"-XL " + masks} \
                -truth ${truth_vcf} -eval ${eval_vcf} \
                -tpfn "tpfn.vcf" \
                -tpfp "tpfp.vcf" \
                -ftnfn "ftnfn.vcf" \
                -filter-analysis "filter-analysis.txt" \
                -summary "summary.txt"

        grep '#' tpfn.vcf > HEAD
        grep STATUS=FN tpfn.vcf > BODY
        cat HEAD BODY > false_negatives.vcf

        grep '#' tpfp.vcf > HEAD
        grep STATUS=FP tpfp.vcf > BODY
        cat HEAD BODY > false_positives.vcf

        grep '#' tpfp.vcf > HEAD
        grep STATUS=TP tpfp.vcf > BODY
        cat HEAD BODY > true_positives.vcf

        grep '#' ftnfn.vcf > HEAD
        grep STATUS=FFN ftnfn.vcf > BODY
        cat HEAD BODY > filtered_false_negatives.vcf

        for vcf in false_negatives.vcf false_positives.vcf true_positives.vcf filtered_false_negatives.vcf ; do
            #HACK: IndexFeatureFile throws error if vcf is empty, which is possible here especially in the case of false negatives
            gatk --java-options "-Xmx2g" SelectVariants -V $vcf -O tmp.vcf
            mv tmp.vcf $vcf
            mv tmp.vcf.idx $vcf.idx
        done
    }

    ####################################################################################
    # Runtime Environment:
    runtime {
        cpu: select_first([cpu, 1])
        memory: machine_mem + " MB"
        bootDiskSizeGb: select_first([disk_space_gb, default_disk_space_gb])
        disks: "local-disk " + select_first([boot_disk_size_gb, default_boot_disk_size_gb]) + if use_ssd then " SSD" else " HDD"
        docker: "${gatk_docker}"
        preemptible: select_first([preemptible_attempts, 2])
    }

    ####################################################################################
    # Outputs:
    output {
        File fn = "false_negatives.vcf"
        File fn_idx = "false_negatives.vcf.idx"
        File fp = "false_positives.vcf"
        File fp_idx = "false_positives.vcf.idx"
        File tp = "true_positives.vcf"
        File tp_idx = "true_positives.vcf.idx"
        File ffn = "filtered_false_negatives.vcf"
        File ffn_idx = "filtered_false_negatives.vcf.idx"
        File summary = "summary.txt"
        File filter_analysis = "filter-analysis.txt"
    }
}
