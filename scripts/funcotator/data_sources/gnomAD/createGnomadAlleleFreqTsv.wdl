# Create a TSV containing genomic position, dbSNP ID, alleles, and the allele frequency from v2.1 of gnomAD (hg19/b37).
#
# NOTE: This will by default download all of gnomAD to disk.  This is a big file, so be careful!
#
# Description of inputs:
#
#   Required:
#     String gatk_docker             -  GATK Docker image in which to run
#     File gnomAD_file               -  gnomAD VCF file to process
#     String out_file_name           -  Output file name.
#
#   Optional:
#     File gatk4_jar_override        -  Override Jar file containing GATK 4.  Use this when overriding the docker JAR or when using a backend without docker.
#     Int  mem                       -  Amount of memory to give to the machine running each task in this workflow.
#     Int  preemptible_attempts      -  Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb             -  Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                       -  Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb         -  Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
# This WDL needs to decide whether to use the ``gatk_jar`` or ``gatk_jar_override`` for the jar location.  As of cromwell-0.24,
# this logic *must* go into each task.  Therefore, there is a lot of duplicated code.  This allows users to specify a jar file
# independent of what is in the docker file.  See the README.md for more info.
#
workflow CreateGnomadAlleleFreqTsv {

    File gnomAD_file = "gs://gnomad-public/release/2.1/vcf/genomes/gnomad.genomes.r2.1.sites.vcf.bgz"
    String out_file_name

    String gatk_docker

    Int?  mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    call CreateGnomadAlleleFreqTsvTask {
        input:
            gnomAD_file               = gnomAD_file,
            out_file_name             = out_file_name,
            gatk_docker               = gatk_docker,
            mem_gb                    = mem_gb,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb         = boot_disk_size_gb
    }

    output {
        File gnomadAlleleFreqTsv = CreateGnomadAlleleFreqTsvTask.gnomadAlleleFreqTsv
    }
}


task CreateGnomadAlleleFreqTsvTask {

    File gnomAD_file
    String out_file_name

     # ------------------------------------------------
     # runtime
     String gatk_docker
     Int? mem_gb
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
     Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

     # ------------------------------------------------
     # Run our command:
     command <<<
         set -e

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > timingInformation.txt

         cat ${gnomAD_file} | sed 's#^\([0-9X]*\)\t\([0-9]*\)\t\(.*\)\t\([ATGCN]*\)\t\([ATGCN,]*\)\t.*;AF=\([e0-9\.+\-]*\);.*#\1 \2 \3 \4 \5 \6#g' > ${out_file_name}

         endTime=`date +%s.%N`
         echo "EndTime: $endTime" >> timingInformation.txt
         elapsedTime=`echo "scale=5;$endTime - $startTime" | bc`
         echo "Elapsed Time: $elapsedTime" >> timingInformation.txt
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
         File gnomadAlleleFreqTsv = "${out_file_name}"
     }
 }
