
task indexFastaFileTask {

  File input_fasta_file

  # runtime
  String gatk_docker
  File? gatk_override
  Int? mem
  Int? preemptible_attempts
  Int? disk_space_gb
  Int? cpu
  Int? boot_disk_size_gb

  # ------------------------------------------------
  # Get machine settings:
  Boolean use_ssd = false

  # You may have to change the following two parameter values depending on the task requirements
  Int default_ram_mb = 1024 * 3
  # WARNING: In the workflow, you should calculate the disk space as an input to this task (disk_space_gb).  Please see [TODO: Link from Jose] for examples.
  Int default_disk_space_gb = 100

  Int default_boot_disk_size_gb = 15

  # Mem is in units of GB but our command and memory runtime values are in MB
  Int machine_mem = if defined(mem) then mem * 1024 else default_ram_mb
  Int command_mem = machine_mem - 1024

  # ------------------------------------------------
  # Run our command:
  command <<<

    set -e
    samtools faidx ${input_fasta_file}

  >>>

  # ------------------------------------------------
  # Runtime settings:
  runtime {
      docker: gatk_docker
      memory: machine_mem + " MB"
      disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
      bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
      preemptible: 0
      cpu: select_first([cpu, 1])
  }

  # ------------------------------------------------
  # Outputs:
  output {
      File vcf_index         = "${input_fasta_file}.idx"
  }

}