#  Run Mutect 2 on a single tumor-normal pair or on a single tumor sample.
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fasta_index, ref_dict: reference genome, index, and dictionary
#  tumor_bam, tumor_bam_index, and tumor_sample_name: self-explanatory
#  normal_bam, normal_bam_index, and normal_sample_name: self-explanatory
#  pon, pon_index: optional panel of normals and index in vcf format containing known false positves
#  scatter_count: number of parallel jobs when scattering over intervals
#  gnomad, gnomad_index: optional database of known germline variants, obtainable from http://gnomad.broadinstitute.org/downloads
#  variants_for_contamination, variants_for_contamination_index: vcf of common variants with allele frequencies fo calculating contamination
#  is_run_orientation_bias_filter: if true, run the orientation bias filter post-processing step
#  is_run_oncotator: if true, annotate the M2 VCFs using oncotator (to produce a TCGA MAF).  Important:  This requires a docker image and should
#   not be run in environments where docker is unavailable (e.g. SGE cluster on a Broad on-prem VM).  Access to docker
#   hub is also required, since the task will download a public docker image.
#
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
#   this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
#   independent of what is in the docker file.  See the README.md for more info.
#

workflow Mutect2 {
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File tumor_bam
  File tumor_bam_index
  String tumor_sample_name
  File? normal_bam
  File? normal_bam_index
  String? normal_sample_name
  File? pon
  File? pon_index
  Int scatter_count
  File? gnomad
  File? gnomad_index
  File? variants_for_contamination
  File? variants_for_contamination_index
  Boolean is_run_orientation_bias_filter
  Boolean is_run_oncotator
  String m2_docker
  String basic_bash_docker = "ubuntu:16.04"
  String oncotator_docker
  File? gatk4_jar_override
  Int preemptible_attempts
  File? onco_ds_tar_gz
  String? onco_ds_local_db_dir
  Array[String] artifact_modes
  File picard_jar
  String? m2_extra_args
  String? m2_extra_filtering_args
  String? sequencing_center
  String? sequence_source
  File? default_config_file
  File? tumor_sequencing_artifact_metrics
  Boolean is_bamOut = false

  # Do not populate unless you know what you are doing...
  File? auth
  # This value is added to every tasks disk in case the dynamic sizing isnt enough
  Int? emergency_extra_disk

  # Disk sizes used for dynamic sizing
  Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fasta_index, "GB"))
  Int tumor_bam_size = ceil(size(tumor_bam, "GB") + size(tumor_bam_index, "GB"))
  Int gnomad_vcf_size = if defined(gnomad) then ceil(size(gnomad, "GB") + size(gnomad_index, "GB")) else 0
  Int normal_bam_size = if defined(normal_bam) then ceil(size(normal_bam, "GB") + size(normal_bam_index, "GB")) else 0
  # If no tar is provided, the task downloads one from broads ftp server
  Int onco_tar_size = if defined(onco_ds_tar_gz) then ceil(size(onco_ds_tar_gz, "GB") * 3) else 100
  Int gatk4_override_size = if defined(gatk4_jar_override) then ceil(size(gatk4_jar_override, "GB")) else 0

  # This is added to every task as padding, should increase if systematically you need more disk for every call
  Int disk_pad = 10 + ceil(size(picard_jar, "GB")) + gatk4_override_size + select_first([emergency_extra_disk,0])

  # These are multipliers to multipler inputs by to make sure we have enough disk to accommodate for possible output sizes
  # Large is for Bams/WGS vcfs
  # Small is for metrics/other vcfs
  Float large_input_to_output_mulitplier = 2.25
  Float small_input_to_output_mulitplier = 2.0

  call ProcessOptionalArguments {
    input:
      tumor_sample_name = tumor_sample_name,
      normal_bam = normal_bam,
      normal_sample_name = normal_sample_name,
      preemptible_attempts = preemptible_attempts,
      docker = basic_bash_docker
  }

  Int split_interval_size = ref_size + ceil(size(intervals, "GB") * small_input_to_output_mulitplier) + disk_pad

  call SplitIntervals {
    input:
      scatter_count = scatter_count,
      intervals = intervals,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ref_dict = ref_dict,
      gatk4_jar_override = gatk4_jar_override,
      preemptible_attempts = preemptible_attempts,
      m2_docker = m2_docker,
      disk_space_gb = split_interval_size
  }

  # Size M2 differently based on if we are using NIO or not
  Int m2_output_size = tumor_bam_size / scatter_count
  Int m2_nio_disk_size = ((tumor_bam_size + normal_bam_size) / scatter_count) + ref_size + (gnomad_vcf_size / scatter_count) + m2_output_size + disk_pad
  Int m2_no_nio_disk_size = tumor_bam_size + normal_bam_size + ref_size + gnomad_vcf_size + m2_output_size + disk_pad
  Int m2_per_scatter_size = if defined(auth) then m2_nio_disk_size else m2_no_nio_disk_size

  scatter (subintervals in SplitIntervals.interval_files ) {
    call M2 {
      input:
        intervals = subintervals,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        tumor_bam = tumor_bam,
        tumor_bam_index = tumor_bam_index,
        normal_bam = normal_bam,
        normal_bam_index = normal_bam_index,
        pon = pon,
        gnomad = gnomad,
        output_vcf_name = ProcessOptionalArguments.output_name,
        gatk4_jar_override = gatk4_jar_override,
        preemptible_attempts = preemptible_attempts,
        m2_docker = m2_docker,
        m2_extra_args = m2_extra_args,
        is_bamOut = is_bamOut,
        disk_space_gb = m2_per_scatter_size,
        auth = auth
    }

    Float sub_vcf_size = size(M2.output_vcf, "GB")
    Float sub_bamout_size = size(M2.output_bamOut, "GB")
  }

  call SumFloats as SumSubVcfs {
    input:
      sizes = sub_vcf_size,
      preemptible_attempts = preemptible_attempts,
  }

  Int merge_vcf_size = ceil(SumSubVcfs.total_size * large_input_to_output_mulitplier) + disk_pad

  call MergeVCFs {
    input:
      input_vcfs = M2.output_vcf,
      output_vcf_name = ProcessOptionalArguments.output_name,
      gatk4_jar_override = gatk4_jar_override,
      preemptible_attempts = preemptible_attempts,
      disk_space_gb = merge_vcf_size,
      m2_docker = m2_docker
  }

  if (is_bamOut) {
    call SumFloats as SumSubBamouts {
      input:
        sizes = sub_bamout_size,
        preemptible_attempts = preemptible_attempts,
    }

    Int merge_bamout_size = ceil(SumSubBamouts.total_size * large_input_to_output_mulitplier) + disk_pad

    call MergeBamOuts {
      input:
        bam_outs = M2.output_bamOut,
        picard_jar = picard_jar,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        gatk4_jar = gatk4_jar,
        gatk4_jar_override = gatk4_jar_override,
        m2_docker = m2_docker,
        output_vcf_name = ProcessOptionalArguments.output_name,
        disk_space_gb = merge_bamout_size
    }
  }

  if (is_run_orientation_bias_filter && !defined(tumor_sequencing_artifact_metrics)) {
    Int seq_artifact_metrcis_size = tumor_bam_size + ref_size + disk_pad
    call CollectSequencingArtifactMetrics {
      input:
        preemptible_attempts = preemptible_attempts,
        m2_docker = m2_docker,
        tumor_bam = tumor_bam,
        tumor_bam_index = tumor_bam_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        picard_jar = picard_jar,
        disk_space_gb = seq_artifact_metrcis_size
      }
  }

  # If contamination vcf is defined, run CalculateContamination
  if (defined(variants_for_contamination)) {
    Int calculate_contam_size = tumor_bam_size + ceil(size(variants_for_contamination, "GB") * small_input_to_output_mulitplier) + disk_pad
    call CalculateContam {
      input:
        gatk4_jar_override = gatk4_jar_override,
        intervals = intervals,
        preemptible_attempts = preemptible_attempts,
        m2_docker = m2_docker,
        tumor_bam = tumor_bam,
        tumor_bam_index = tumor_bam_index,
        variants_for_contamination = variants_for_contamination,
        variants_for_contamination_index = variants_for_contamination_index,
        disk_space_gb = calculate_contam_size,
        auth = auth
    }
  }

  Int filter_size = ceil(size(MergeVCFs.output_vcf, "GB") * small_input_to_output_mulitplier) + disk_pad

  call Filter {
    input:
      gatk4_jar_override = gatk4_jar_override,
      unfiltered_vcf = MergeVCFs.output_vcf,
      output_vcf_name = ProcessOptionalArguments.output_name,
      m2_docker = m2_docker,
      preemptible_attempts = preemptible_attempts,
      m2_extra_filtering_args = m2_extra_filtering_args,
      contamination_table = CalculateContam.contamination_table,
      disk_space_gb = filter_size,
      auth = auth
  }

  if (is_run_orientation_bias_filter) {
    # Get the metrics either from the workflow input or CollectSequencingArtifactMetrics if no workflow input is provided
    File input_artifact_metrics = select_first([tumor_sequencing_artifact_metrics, CollectSequencingArtifactMetrics.pre_adapter_metrics])
    Int filter_by_orientation_bias_size = ceil(size(Filter.filtered_vcf, "GB") * small_input_to_output_mulitplier) + ceil(size(input_artifact_metrics, "GB")) + disk_pad

    call FilterByOrientationBias {
      input:
        gatk4_jar_override = gatk4_jar_override,
        input_vcf = Filter.filtered_vcf,
        output_vcf_name = ProcessOptionalArguments.output_name,
        m2_docker = m2_docker,
        preemptible_attempts = preemptible_attempts,
        pre_adapter_metrics = input_artifact_metrics,
        artifact_modes = artifact_modes,
        disk_space_gb = filter_by_orientation_bias_size,
        auth = auth
    }
  }

  if (is_run_oncotator) {
    File oncotate_vcf_input = select_first([FilterByOrientationBias.filtered_vcf, Filter.filtered_vcf])
    Int oncotate_size = ceil(size(oncotate_vcf_input, "GB") * large_input_to_output_mulitplier) + onco_tar_size + disk_pad

    call OncotateM2 {
      input:
        m2_vcf = oncotate_vcf_input,
        case_id = tumor_sample_name,
        control_id = normal_sample_name,
        preemptible_attempts = preemptible_attempts,
        oncotator_docker = oncotator_docker,
        onco_ds_tar_gz = onco_ds_tar_gz,
        onco_ds_local_db_dir = onco_ds_local_db_dir,
        sequencing_center = sequencing_center,
        sequence_source = sequence_source,
        default_config_file = default_config_file,
        disk_space_gb = oncotate_size
    }
  }

  output {
    File unfiltered_vcf = MergeVCFs.output_vcf
    File unfiltered_vcf_index = MergeVCFs.output_vcf_index
    File filtered_vcf = select_first([FilterByOrientationBias.filtered_vcf, Filter.filtered_vcf])
    File filtered_vcf_index = select_first([FilterByOrientationBias.filtered_vcf_index, Filter.filtered_vcf_index])
    Array[String] tumor_bam_sample_names = M2.tumor_bam_sample_name
    String tumor_bam_sample_name = tumor_bam_sample_names[0]
    Array[String] normal_bam_sample_names = M2.normal_bam_sample_name
    String normal_bam_sample_name = normal_bam_sample_names[0]

    # select_first() fails if nothing resolves to non-null, so putting in "null" for now.
    File contamination_table = select_first([CalculateContam.contamination_table, "null"])
    File oncotated_m2_maf = select_first([OncotateM2.oncotated_m2_maf, "null"])
    File preadapter_detail_metrics = select_first([CollectSequencingArtifactMetrics.pre_adapter_metrics, "null"])
    File bamout = select_first([MergeBamOuts.merged_bam_out, "null"])
    File bamout_index = select_first([MergeBamOuts.merged_bam_out_index, "null"])
  }
}

task M2 {
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  String tumor_bam
  String tumor_bam_index
  String? normal_bam
  String? normal_bam_index
  String? pon
  String? gnomad
  String output_vcf_name
  File? gatk4_jar_override
  String? m2_extra_args
  Boolean? is_bamOut

  # Runtime parameters
  String m2_docker
  Int? preemptible_attempts
  Int disk_space_gb
  Int? mem

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  # Do not populate this unless you know what you are doing...
  File? auth

  command {
  if [[ "${auth}" == *.json ]]; then
    gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
    GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
    export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
  fi

  # Use GATK Jar override if specified
  GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

  # These two following commands will only be executed if normal_bam is defined, otherwise they will be blank
  ${"java -Xmx" + command_mem + "m -jar $GATK_JAR GetSampleName -I " + normal_bam + " -O normal_name.txt"}
  ${"normal_command_line=\"-I " + normal_bam + " -normal `cat normal_name.txt`\""}

  java -Xmx${command_mem}m -jar $GATK_JAR GetSampleName -I ${tumor_bam} -O tumor_name.txt

  # We need to create a file regardless, even if it stays empty
  touch bamout.bam
  touch normal_name.txt

  java -Xmx${command_mem}m -jar $GATK_JAR Mutect2 \
    -R ${ref_fasta} \
    -I ${tumor_bam} \
    -tumor `cat tumor_name.txt` \
    $normal_command_line \
    ${"--germline_resource " + gnomad} \
    ${"--normal_panel " + pon} \
    ${"-L " + intervals} \
    -O "${output_vcf_name}.vcf" \
    ${true='--bamOutput bamout.bam' false='' is_bamOut} \
    ${m2_extra_args}
  }

  runtime {
    docker: m2_docker
    memory: machine_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    File output_vcf = "${output_vcf_name}.vcf"
    File output_bamOut = "bamout.bam"
    String tumor_bam_sample_name = read_string("tumor_name.txt")
    String normal_bam_sample_name = read_string("normal_name.txt")
  }
}

# HACK: cromwell can't handle the optional normal sample name in output or input --
# string interpolation of optionals only works inside a command block
# thus we use this hack
task ProcessOptionalArguments {
  String tumor_sample_name
  String? normal_bam
  String? normal_sample_name

  # Runtime parameters
  Int? preemptible_attempts
  String docker

  command {
    if [[ "_${normal_bam}" == *.bam ]]; then
      echo "${tumor_sample_name}-vs-${normal_sample_name}" > name.tmp
    else
      echo "${tumor_sample_name}-tumor-only" > name.tmp
    fi
  }

  runtime {
    docker: docker
    memory: "1 GB"
    disks: "local-disk " + 10 + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
      String output_name = read_string("name.tmp")
  }
}

task MergeVCFs {
  Array[File] input_vcfs
  String output_vcf_name
  File? gatk4_jar_override

  # Runtime parameters
  Int? preemptible_attempts
  String m2_docker
  Int disk_space_gb
  Int? mem

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  command {
    java -Xmx${command_mem}m -jar ${default="/root/gatk.jar" gatk4_jar_override} MergeVcfs \
      -I ${sep=' -I ' input_vcfs} \
      -O ${output_vcf_name}.vcf
  }

  runtime {
    docker: m2_docker
    memory: machine_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    File output_vcf = "${output_vcf_name}.vcf"
    File output_vcf_index = "${output_vcf_name}.vcf.idx"
  }
}

task CollectSequencingArtifactMetrics {
  File tumor_bam
  File tumor_bam_index
  File ref_fasta
  File ref_fasta_index
  File picard_jar

  # Runtime parameters
  String m2_docker
  Int? preemptible_attempts
  Int disk_space_gb
  Int? mem

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 7000
  Int command_mem = machine_mem - 500

  command {
    # Sed'ing of the metrics file to make it GATK compatitble has been moved to the FilterByOrientationBias task

    java -Xmx${command_mem}m -jar ${picard_jar} CollectSequencingArtifactMetrics \
      I=${tumor_bam} \
      O="gatk" \
      R=${ref_fasta} \
      VALIDATION_STRINGENCY=LENIENT
  }

  runtime {
    docker: m2_docker
    memory: machine_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    File pre_adapter_metrics = "gatk.pre_adapter_detail_metrics"
  }
}

task CalculateContam {
  File? gatk4_jar_override
  File? intervals
  String tumor_bam
  String tumor_bam_index
  File? variants_for_contamination
  File? variants_for_contamination_index

  # Runtime parameters
  Int? preemptible_attempts
  String m2_docker
  Int disk_space_gb
  Int? mem

  # Do not populate this unless you know what you are doing...
  File? auth

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 7000
  Int command_mem = machine_mem - 500

  command {
    set -e

    if [[ "${auth}" == *.json ]]; then
      gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
      GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
      export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
    fi

    # Use GATK Jar override if specified
    GATK_JAR=${default="/root/gatk.jar" gatk4_jar_override}

    java -Xmx${command_mem}m -jar $GATK_JAR GetPileupSummaries -I ${tumor_bam} ${"-L " + intervals} -V ${variants_for_contamination} -O pileups.table
    java -Xmx${command_mem}m -jar $GATK_JAR CalculateContamination -I pileups.table -O contamination.table
  }

  runtime {
    docker: m2_docker
    memory: command_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    File contamination_table = "contamination.table"
  }
}

task Filter {
  File? gatk4_jar_override
  File unfiltered_vcf
  String output_vcf_name
  String? m2_extra_filtering_args
  File? contamination_table

  # Runtime parameters
  Int? preemptible_attempts
  String m2_docker
  Int disk_space_gb
  Int? mem

  # Do not populate this unless you know what you are doing...
  File? auth

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 7000
  Int command_mem = machine_mem - 500

  command {
    set -e

    if [[ "${auth}" == *.json ]]; then
      gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
      GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
      export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
    fi

    java -Xmx${command_mem}m -jar ${default="/root/gatk.jar" gatk4_jar_override} FilterMutectCalls -V ${unfiltered_vcf} \
      -O ${output_vcf_name}-filtered.vcf ${"-contaminationTable " + contamination_table} \
      ${m2_extra_filtering_args}
  }

  runtime {
    docker: m2_docker
    memory: command_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    File filtered_vcf = "${output_vcf_name}-filtered.vcf"
    File filtered_vcf_index = "${output_vcf_name}-filtered.vcf.idx"
  }
}

task FilterByOrientationBias {
  File? gatk4_jar_override
  File input_vcf
  String output_vcf_name
  File pre_adapter_metrics
  Array[String]? artifact_modes

  # Runtime parameters
  Int? preemptible_attempts
  String m2_docker
  Int disk_space_gb
  Int? mem

  # Do not populate this unless you know what you are doing...
  File? auth

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 7000
  Int command_mem = machine_mem - 500

  command {
    set -e

    if [[ "${auth}" == *.json ]]; then
      gsutil cp ${auth} /root/.config/gcloud/application_default_credentials.json
      GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
      export GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json
    fi

    # Convert to GATK format
    sed -r "s/picard\.analysis\.artifacts\.SequencingArtifactMetrics\\\$PreAdapterDetailMetrics/org\.broadinstitute\.hellbender\.tools\.picard\.analysis\.artifacts\.SequencingArtifactMetrics\$PreAdapterDetailMetrics/g" \
      "${pre_adapter_metrics}" > "gatk.pre_adapter_detail_metrics"

     java -Xmx${command_mem}m -jar ${default="/root/gatk.jar" gatk4_jar_override} FilterByOrientationBias -A ${sep=" -A " artifact_modes} \
       -V ${input_vcf} -P gatk.pre_adapter_detail_metrics --output "${output_vcf_name}-filtered.vcf"
  }

  runtime {
    docker: m2_docker
    memory: command_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    File filtered_vcf = "${output_vcf_name}-filtered.vcf"
    File filtered_vcf_index = "${output_vcf_name}-filtered.vcf.idx"
  }
}

task SplitIntervals {
  Int scatter_count
  File? intervals
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? gatk4_jar_override

  # Runtime parameters
  String m2_docker
  Int disk_space_gb
  Int? preemptible_attempts
  Int? mem

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  command {
    set -e

    mkdir interval-files

    java -Xmx${command_mem}m -jar ${default="/root/gatk.jar" gatk4_jar_override} SplitIntervals \
      -R ${ref_fasta} \
      ${"-L " + intervals} \
      -scatter ${scatter_count} \
      -O interval-files

    cp interval-files/*.intervals .
  }

  runtime {
    docker: m2_docker
    memory: command_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    Array[File] interval_files = glob("*.intervals")
  }
}

task MergeBamOuts {
  String gatk4_jar
  Array[File]+ bam_outs
  File picard_jar
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  File? gatk4_jar_override
  String output_vcf_name

  # Runtime parameters
  Int? mem
  String m2_docker
  Int? preemptible_attempts
  Int disk_space_gb

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 7000
  Int command_mem = machine_mem - 500

  command <<<
    # This command block assumes that there is at least one file in bam_outs.
    #  Do not call this task if len(bam_outs) == 0
    set -e

    java -Xmx${command_mem}m -jar ${picard_jar} GatherBamFiles \
      I=${sep=" I=" bam_outs} \
      O=${output_vcf_name}.out.bam \
      R=${ref_fasta}

    samtools index ${output_vcf_name}.out.bam ${output_vcf_name}.out.bam.bai
  >>>

  runtime {
    docker: m2_docker
    memory: machine_mem + " MB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 10])
  }

  output {
    File merged_bam_out = "${output_vcf_name}.out.bam"
    File merged_bam_out_index = "${output_vcf_name}.out.bam.bai"
  }
}

task OncotateM2 {
  File m2_vcf
  String case_id
  String? control_id
  File? onco_ds_tar_gz
  String? onco_ds_local_db_dir
  String? oncotator_exe
  String? sequencing_center
  String? sequence_source
  File? default_config_file

  # Runtime parameters
  Int? preemptible_attempts
  String oncotator_docker
  Int disk_space_gb
  Int? mem

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1000 else 3500
  Int command_mem = machine_mem - 500

  command <<<
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
      rm ${onco_ds_tar_gz}

    else
      echo "Downloading and installing oncotator datasources from Broad FTP site..."
      # Download and untar the db-dir
      wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/oncotator/oncotator_v1_ds_April052016.tar.gz
      tar zxvf oncotator_v1_ds_April052016.tar.gz
      rm oncotator_v1_ds_April052016.tar.gz
      ln -s oncotator_v1_ds_April052016 onco_dbdir
    fi

    ${default="/root/oncotator_venv/bin/oncotator" oncotator_exe} --db-dir onco_dbdir/ -c $HOME/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt  \
      -v ${m2_vcf} ${case_id}.maf.annotated hg19 -i VCF -o TCGAMAF --skip-no-alt --infer-onps --collapse-number-annotations --log_name oncotator.log \
      -a Center:${default="Unknown" sequencing_center} \
      -a source:${default="Unknown" sequence_source} \
      -a normal_barcode:${default=" " control_id} \
      -a tumor_barcode:${case_id} \
      ${"--default_config " + default_config_file}
    >>>

    runtime {
      docker: oncotator_docker
      memory: machine_mem + " MB"
      disks: "local-disk " + disk_space_gb + " HDD"
      preemptible: select_first([preemptible_attempts, 10])
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
