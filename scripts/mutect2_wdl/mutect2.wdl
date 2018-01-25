## Copyright Broad Institute, 2017
##
## This WDL workflow runs GATK4 Mutect 2 on a single tumor-normal pair or on a single tumor sample,
## and performs additional filtering and functional annotation tasks.
##
## Main requirements/expectations :
## - One analysis-ready BAM file (and its index) for each sample
##
## Description of inputs:
##
## ** Runtime ** 
## gatk_docker, oncotator_docker: docker images to use for GATK 4 Mutect2 and for Oncotator
## preemptible_attempts: how many preemptions to tolerate before switching to a non-preemptible machine (on Google)
## gatk_override: (optional) local file or Google bucket path to a GATK 4 java jar file to be used instead of the GATK 4 jar
##                in the docker image.  This must be supplied when running in an environment that does not support docker
##                (e.g. SGE cluster on a Broad on-prem VM)
##
## ** Workflow options **
## intervals: genomic intervals (will be used for scatter)
## scatter_count: number of parallel jobs to generate when scattering over intervals
## artifact_modes: types of artifacts to consider in the orientation bias filter (optional)
## m2_extra_args, m2_extra_filtering_args: additional arguments for Mutect2 calling and filtering (optional)
## split_intervals_extra_args: additional arguments for splitting intervals before scattering (optional)
## run_orientation_bias_filter: if true, run the orientation bias filter post-processing step (optional, false by default)
## run_oncotator: if true, annotate the M2 VCFs using oncotator (to produce a TCGA MAF).  Important:  This requires a
##                   docker image and should  not be run in environments where docker is unavailable (e.g. SGE cluster on
##                   a Broad on-prem VM).  Access to docker hub is also required, since the task downloads a public docker image.
##                   (optional, false by default)
##
## ** Primary inputs **
## ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
## tumor_bam, tumor_bam_index: BAM and index for the tumor sample
## normal_bam, normal_bam_index: BAM and index for the normal sample
##
## ** Primary resources ** (optional but strongly recommended)
## pon, pon_index: optional panel of normals in VCF format containing probable technical artifacts (false positves)
## gnomad, gnomad_index: optional database of known germline variants (see http://gnomad.broadinstitute.org/downloads)
## variants_for_contamination, variants_for_contamination_index: VCF of common variants with allele frequencies for calculating contamination
##
## ** Secondary resources ** (for optional tasks)
## onco_ds_tar_gz, default_config_file: Oncotator datasources and config file
## sequencing_center, sequence_source: metadata for Oncotator
##
## Outputs :
## - One VCF file and its index with primary filtering applied; secondary filtering and functional annotation if requested; a bamout.bam
##   file of reassembled reads if requested
##
## Cromwell version support
## - Successfully tested on v29
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## pages at https://hub.docker.com/r/broadinstitute/* for detailed licensing information
## pertaining to the included programs.
workflow Mutect2 {
    # Mutect2 inputs
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File tumor_bam
    File tumor_bai
    File? normal_bam
    File? normal_bai
    File? pon
    File? pon_index
    Int scatter_count
    File? gnomad
    File? gnomad_index
    File? variants_for_contamination
    File? variants_for_contamination_index
    Boolean? run_orientation_bias_filter
    Boolean run_ob_filter = select_first([run_orientation_bias_filter, false])
    Array[String]? artifact_modes
    File? tumor_sequencing_artifact_metrics
    String? m2_extra_args
    String? m2_extra_filtering_args
    String? split_intervals_extra_args
    Boolean? make_bamout
    Boolean generate_bamout = select_first([make_bamout, false])
    Boolean? compress_vcfs
    Boolean compress = select_first([compress_vcfs, false])

    # oncotator inputs
    Boolean? run_oncotator
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    String? sequencing_center
    String? sequence_source
    File? default_config_file

    File? gatk_override

    # runtime
    String gatk_docker
    String basic_bash_docker = "ubuntu:16.04"
    String? oncotator_docker
    Int? preemptible_attempts

    # Do not populate unless you know what you are doing...
    File? auth

    # Use as a last resort to increase the disk given to every task in case of ill behaving data
    Int? emergency_extra_disk

    # Disk sizes used for dynamic sizing
    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fai, "GB"))
    Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bai, "GB"))
    Int gnomad_vcf_size = if defined(gnomad) then ceil(size(gnomad, "GB") + size(gnomad_index, "GB")) else 0
    Int normal_bam_size = if defined(normal_bam) then ceil(size(normal_bam, "GB") + size(normal_bai, "GB")) else 0

    # If no tar is provided, the task downloads one from broads ftp server
    Int onco_tar_size = if defined(onco_ds_tar_gz) then ceil(size(onco_ds_tar_gz, "GB") * 3) else 100
    Int gatk_override_size = if defined(gatk_override) then ceil(size(gatk_override, "GB")) else 0

    # This is added to every task as padding, should increase if systematically you need more disk for every call
    Int disk_pad = 10 + gatk_override_size + select_first([emergency_extra_disk,0])

    # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
    # Large is for Bams/WGS vcfs
    # Small is for metrics/other vcfs
    Float large_input_to_output_multiplier = 2.25
    Float small_input_to_output_multiplier = 2.0

    # logic about output file names -- these are the names *without* .vcf extensions
    String output_basename = basename(tumor_bam, ".bam")
    String unfiltered_name = output_basename + "-unfiltered"
    String filtered_name = output_basename + "-filtered"


    String output_vcf_name = basename(tumor_bam, ".bam") + ".vcf"


    call SplitIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            split_intervals_extra_args = split_intervals_extra_args,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts,
            disk_space = ref_size + ceil(size(intervals, "GB") * small_input_to_output_multiplier) + disk_pad
    }

    Int m2_output_size = tumor_bam_size / scatter_count
    scatter (subintervals in SplitIntervals.interval_files ) {
        call M2 {
            input:
                intervals = subintervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                pon = pon,
                pon_index = pon_index,
                gnomad = gnomad,
                gnomad_index = gnomad_index,
                preemptible_attempts = preemptible_attempts,
                m2_extra_args = m2_extra_args,
                generate_bamout = generate_bamout,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts,
                disk_space = tumor_bam_size + normal_bam_size + ref_size + gnomad_vcf_size + m2_output_size + disk_pad,
                auth = auth
        }

        Float sub_vcf_size = size(M2.output_vcf, "GB")
        Float sub_bamout_size = size(M2.output_bamOut, "GB")
    }

    call SumFloats as SumSubVcfs {
        input:
            sizes = sub_vcf_size,
            preemptible_attempts = preemptible_attempts
    }

    call MergeVCFs {
        input:
            input_vcfs = M2.output_vcf,
            output_name = unfiltered_name,
            compress = compress,
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts,
            disk_space = ceil(SumSubVcfs.total_size * large_input_to_output_multiplier) + disk_pad
    }

    if (generate_bamout) {
        call SumFloats as SumSubBamouts {
            input:
                sizes = sub_bamout_size,
                preemptible_attempts = preemptible_attempts
        }

        call MergeBamOuts {
            input:
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                bam_outs = M2.output_bamOut,
                output_vcf_name = basename(MergeVCFs.merged_vcf, ".vcf"),
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                disk_space = ceil(SumSubBamouts.total_size * large_input_to_output_multiplier) + disk_pad
        }
    }

    if (run_ob_filter && !defined(tumor_sequencing_artifact_metrics)) {
        call CollectSequencingArtifactMetrics {
            input:
                gatk_docker = gatk_docker,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                preemptible_attempts = preemptible_attempts,
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                gatk_override = gatk_override,
                disk_space = tumor_bam_size + ref_size + disk_pad
        }
    }

    if (defined(variants_for_contamination)) {
        call CalculateContamination {
            input:
                gatk_override = gatk_override,
                intervals = intervals,
                preemptible_attempts = preemptible_attempts,
                gatk_docker = gatk_docker,
                tumor_bam = tumor_bam,
                tumor_bai = tumor_bai,
                normal_bam = normal_bam,
                normal_bai = normal_bai,
                variants_for_contamination = variants_for_contamination,
                variants_for_contamination_index = variants_for_contamination_index,
                disk_space = tumor_bam_size + ceil(size(variants_for_contamination, "GB") * small_input_to_output_multiplier) + disk_pad,
                auth = auth
        }
    }

    call Filter {
        input:
            gatk_override = gatk_override,
            gatk_docker = gatk_docker,
            intervals = intervals,
            unfiltered_vcf = MergeVCFs.merged_vcf,
            output_name = filtered_name,
            compress = compress,
            preemptible_attempts = preemptible_attempts,
            contamination_table = CalculateContamination.contamination_table,
            m2_extra_filtering_args = m2_extra_filtering_args,
            disk_space = ceil(size(MergeVCFs.merged_vcf, "GB") * small_input_to_output_multiplier) + disk_pad,
            auth = auth
    }

    if (run_ob_filter) {
        # Get the metrics either from the workflow input or CollectSequencingArtifactMetrics if no workflow input is provided
        File input_artifact_metrics = select_first([tumor_sequencing_artifact_metrics, CollectSequencingArtifactMetrics.pre_adapter_metrics])

        call FilterByOrientationBias {
            input:
                gatk_override = gatk_override,
                input_vcf = Filter.filtered_vcf,
                output_name = filtered_name,
                compress = compress,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts,
                pre_adapter_metrics = input_artifact_metrics,
                artifact_modes = artifact_modes,
                disk_space = ceil(size(Filter.filtered_vcf, "GB") * small_input_to_output_multiplier) + ceil(size(input_artifact_metrics, "GB")) + disk_pad,
                auth = auth
        }
    }


    if (select_first([run_oncotator, false])) {
        File oncotate_vcf_input = select_first([FilterByOrientationBias.filtered_vcf, Filter.filtered_vcf])
        call oncotate_m2 {
            input:
                m2_vcf = oncotate_vcf_input,
                onco_ds_tar_gz = onco_ds_tar_gz,
                onco_ds_local_db_dir = onco_ds_local_db_dir,
                sequencing_center = sequencing_center,
                sequence_source = sequence_source,
                default_config_file = default_config_file,
                case_id = M2.tumor_sample[0],
                control_id = M2.normal_sample[0],
                oncotator_docker = select_first([oncotator_docker, "NO_ONCOTATOR_DOCKER_GIVEN"]),
                preemptible_attempts = preemptible_attempts,
                disk_space = ceil(size(oncotate_vcf_input, "GB") * large_input_to_output_multiplier) + onco_tar_size + disk_pad
        }
    }

    output {
        File unfiltered_vcf = MergeVCFs.merged_vcf
        File unfiltered_vcf_index = MergeVCFs.merged_vcf_index
        File filtered_vcf = select_first([FilterByOrientationBias.filtered_vcf, Filter.filtered_vcf])
        File filtered_vcf_index = select_first([FilterByOrientationBias.filtered_vcf_index, Filter.filtered_vcf_index])
        File? contamination_table = CalculateContamination.contamination_table

        File? oncotated_m2_maf = oncotate_m2.oncotated_m2_maf
        File? preadapter_detail_metrics = CollectSequencingArtifactMetrics.pre_adapter_metrics
        File? bamout = MergeBamOuts.merged_bam_out
        File? bamout_index = MergeBamOuts.merged_bam_out_index
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
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500


    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        mkdir interval-files
        gatk --java-options "-Xmx${command_mem}m" SplitIntervals \
            -R ${ref_fasta} \
            ${"-L " + intervals} \
            -scatter ${scatter_count} \
            -O interval-files \
            ${split_intervals_extra_args}
        cp interval-files/*.intervals .
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        Array[File] interval_files = glob("*.intervals")
    }
}

task M2 {
    # inputs
    File? intervals
    File ref_fasta
    File ref_fai
    File ref_dict
    File tumor_bam
    File tumor_bai
    File? normal_bam
    File? normal_bai
    File? pon
    File? pon_index
    File? gnomad
    File? gnomad_index
    String? m2_extra_args
    Boolean? generate_bamout

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    # Do not populate this unless you know what you are doing...
    File? auth

    command <<<
        set -e

        if [[ "${auth}" == *.json ]]; then
            gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
            GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
            export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
        fi

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        echo "" > normal_name.txt

        gatk --java-options "-Xmx${command_mem}m" GetSampleName -I ${tumor_bam} -O tumor_name.txt
        tumor_command_line="-I ${tumor_bam} -tumor `cat tumor_name.txt`"

        if [[ "${normal_bam}" == *.bam ]]; then
            gatk --java-options "-Xmx${command_mem}m" GetSampleName -I ${normal_bam} -O normal_name.txt
            normal_command_line="-I ${normal_bam} -normal `cat normal_name.txt`"
        fi

        gatk --java-options "-Xmx${command_mem}m" Mutect2 \
            -R ${ref_fasta} \
            $tumor_command_line \
            $normal_command_line \
            ${"--germline-resource " + gnomad} \
            ${"-pon " + pon} \
            ${"-L " + intervals} \
            -O "output.vcf.gz" \
            ${true='--bam-output bamout.bam' false='' generate_bamout} \
            ${m2_extra_args}
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File output_vcf = "output.vcf.gz"
        File output_bamOut = "bamout.bam"
        String tumor_sample = read_string("tumor_name.txt")
        String normal_sample = read_string("normal_name.txt")
    }
}

task MergeVCFs {
    # inputs
    Array[File] input_vcfs
    String output_name
    Boolean compress
    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"


    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 1

    # using MergeVcfs instead of GatherVcfs so we can create indices
    # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.
    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx${command_mem}m" MergeVcfs -I ${sep=' -I ' input_vcfs} -O ${output_vcf}
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File merged_vcf = "${output_vcf}"
        File merged_vcf_index = "${output_vcf_index}"
    }
}

task MergeBamOuts {
    # inputs
    File ref_fasta
    File ref_fai
    File ref_dict
    Array[File]+ bam_outs
    String output_vcf_name

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 1

    command <<<
        # This command block assumes that there is at least one file in bam_outs.
        #  Do not call this task if len(bam_outs) == 0
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx${command_mem}m" GatherBamFiles \
            -I ${sep=" -I " bam_outs} -O ${output_vcf_name}.out.bam -R ${ref_fasta}
        samtools index ${output_vcf_name}.out.bam ${output_vcf_name}.out.bam.bai
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File merged_bam_out = "${output_vcf_name}.out.bam"
        File merged_bam_out_index = "${output_vcf_name}.out.bam.bai"
    }
}

task CollectSequencingArtifactMetrics {
    # inputs
    File ref_fasta
    File ref_fai
    File tumor_bam
    File tumor_bai

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 1

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx${command_mem}m" CollectSequencingArtifactMetrics \
            -I ${tumor_bam} -O "gatk" -R ${ref_fasta} -VALIDATION_STRINGENCY LENIENT
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File pre_adapter_metrics = "gatk.pre_adapter_detail_metrics"
    }
}

task CalculateContamination {
    # inputs
    File? intervals
    File tumor_bam
    File tumor_bai
    File? normal_bam
    File? normal_bai
    File? variants_for_contamination
    File? variants_for_contamination_index

    File? gatk_override

    # runtime
    Int? preemptible_attempts
    String gatk_docker
    Int? disk_space
    Int? mem

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 500

    # Do not populate this unless you know what you are doing...
    File? auth

    command {
        set -e

        if [[ "${auth}" == *.json ]]; then
            gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
            GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
            export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
        fi

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        if [[ "${normal_bam}" == *.bam ]]; then
            gatk --java-options "-Xmx${command_mem}m" GetPileupSummaries -I ${normal_bam} ${"-L " + intervals} -V ${variants_for_contamination} -O normal_pileups.table
            NORMAL_CMD="-matched normal_pileups.table"
        fi

        gatk --java-options "-Xmx${command_mem}m" GetPileupSummaries -I ${tumor_bam} ${"-L " + intervals} -V ${variants_for_contamination} -O pileups.table
        gatk --java-options "-Xmx${command_mem}m" CalculateContamination -I pileups.table -O contamination.table $NORMAL_CMD
    }

    runtime {
        docker: gatk_docker
        memory: command_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
    }

    output {
        File pileups = "pileups.table"
        File contamination_table = "contamination.table"
    }
}

task Filter {
    # inputs
    File? intervals
    File unfiltered_vcf
    String output_name
    Boolean compress
    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
    File? contamination_table
    String? m2_extra_filtering_args

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 500

    # Do not populate this unless you know what you are doing...
    File? auth

    command {
        set -e

        if [[ "${auth}" == *.json ]]; then
            gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
            GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
            export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
        fi

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" FilterMutectCalls -V ${unfiltered_vcf} \
      	    -O ${output_vcf} \
      	    ${"--contamination-table " + contamination_table} \
      	    ${m2_extra_filtering_args}
    }

    runtime {
        docker: gatk_docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File filtered_vcf = "${output_vcf}"
        File filtered_vcf_index = "${output_vcf_index}"
    }
}

task FilterByOrientationBias {
    # input
    File? gatk_override
    File input_vcf
    String output_name
    Boolean compress
    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_index = output_vcf +  if compress then ".tbi" else ".idx"
    File pre_adapter_metrics
    Array[String]? artifact_modes

    # runtime
    Int? preemptible_attempts
    String gatk_docker
    Int? disk_space
    Int? mem
    Int? cpu
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 7000
    Int command_mem = machine_mem - 500

    # Do not populate this unless you know what you are doing...
    File? auth

    command {
        set -e

        if [[ "${auth}" == *.json ]]; then
            gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
            GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
            export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
        fi

        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx${command_mem}m" FilterByOrientationBias \
            -V ${input_vcf} \
            -AM ${sep=" -AM " artifact_modes} \
            -P ${pre_adapter_metrics} \
            -O ${output_vcf}
    }

    runtime {
        docker: gatk_docker
        memory: command_mem + " MB"
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File filtered_vcf = "${output_vcf}"
        File filtered_vcf_index = "${output_vcf_index}"
    }
}

task oncotate_m2 {
    # inputs
    File m2_vcf
    File? onco_ds_tar_gz
    String? onco_ds_local_db_dir
    String? oncotator_exe
    String? sequencing_center
    String? sequence_source
    File? default_config_file
    String case_id
    String? control_id

    # runtime
    String oncotator_docker
    Int? mem
    Int? preemptible_attempts
    Int? disk_space
    Int? cpu
    Boolean use_ssd = false

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem) then mem * 1000 else 3500
    Int command_mem = machine_mem - 500

    command <<<

        # fail if *any* command below (not just the last) doesn't return 0, in particular if wget fails
        set -e

        # local db dir is a directory and has been specified
        if [[ -d "${onco_ds_local_db_dir}" ]]; then
            echo "Using local db-dir: ${onco_ds_local_db_dir}"
            echo "THIS ONLY WORKS WITHOUT DOCKER!"
            ln -s ${onco_ds_local_db_dir} onco_dbdir
        elif [[ "${onco_ds_tar_gz}" == *.tar.gz ]]; then
            echo "Using given tar file: ${onco_ds_tar_gz}"
            mkdir onco_dbdir
            tar zxvf ${onco_ds_tar_gz} -C onco_dbdir --strip-components 1
        else
            echo "Downloading and installing oncotator datasources from Broad FTP site..."
            # Download and untar the db-dir
            wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/oncotator_v1_ds_April052016.tar.gz
            tar zxvf oncotator_v1_ds_April052016.tar.gz
            ln -s oncotator_v1_ds_April052016 onco_dbdir
        fi

        ${default="/root/oncotator_venv/bin/oncotator" oncotator_exe} --db-dir onco_dbdir/ -c $HOME/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt  \
            -v ${m2_vcf} ${case_id}.maf.annotated hg19 -i VCF -o TCGAMAF --skip-no-alt --infer-onps --collapse-number-annotations --log_name oncotator.log \
            -a Center:${default="Unknown" sequencing_center} \
            -a source:${default="Unknown" sequence_source} \
            -a normal_barcode:${control_id} \
            -a tumor_barcode:${case_id} \
            ${"--default_config " + default_config_file}
    >>>

    runtime {
        docker: oncotator_docker
        memory: machine_mem + " GB"
        bootDiskSizeGb: 12
        disks: "local-disk " + select_first([disk_space, 100]) + if use_ssd then " SSD" else " HDD"
        preemptible: select_first([preemptible_attempts, 10])
        cpu: select_first([cpu, 1])
    }

    output {
        File oncotated_m2_maf="${case_id}.maf.annotated"
    }
}

# Calculates sum of a list of floats
task SumFloats {
    Array[Float] sizes

    # Runtime parameters
    Int? preemptible_attempts

    command <<<
        python -c "print ${sep="+" sizes}"
    >>>

    output {
        Float total_size = read_float(stdout())
    }

    runtime {
        docker: "python:2.7"
        disks: "local-disk " + 10 + " HDD"
        preemptible: select_first([preemptible_attempts, 10])
    }
}

