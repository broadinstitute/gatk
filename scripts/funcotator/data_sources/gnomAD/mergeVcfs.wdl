# Run MergeVCFs on a list of VCF files. 
#
# Description of inputs:
#
#   Required:
#     gatk_docker                    -  GATK Docker image in which to run
#     variant_vcfs                   -  Array of Variant Context Files (VCF) containing the variants.
#     output_vcf_file_name           -  Desired name of the resulting VCF output file.
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
workflow MergeVcfsWorkflow {
    String gatk_docker
    Array[File] variant_vcfs
    String output_vcf_file_name

    File? gatk4_jar_override
    Int?  mem
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    call MergeVcfs {
        input:
            input_vcfs                = variant_vcfs,
            output_vcf_file           = output_vcf_file_name,
            gatk_docker               = gatk_docker,
            gatk_override             = gatk4_jar_override,
            mem                       = mem,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb         = boot_disk_size_gb
    }

    output {
        File vcf_file = MergeVcfs.vcf_file
    }
}


task MergeVcfs {

    # ------------------------------------------------
    # Input args:
     Array[File] input_vcfs

     # Output Names:
     String output_vcf_file

     # Runtime Options:
     String gatk_docker
     File? gatk_override
     Int? mem
     Int? preemptible_attempts
     Int? disk_space_gb
     Int? cpu
     Int? boot_disk_size_gb

    String dollar = "$"

    # ------------------------------------------------
    # Process input args:
    String timing_output_file = basename(output_vcf_file) + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3 * 1024
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
         export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        fileListArgs=""

        # Ensure that the names of the files end in the correct suffixes:
        # (MergeVCFs requires compressed vcfs to end in '.vcf.gz')
        for f in ${sep=' ' input_vcfs} ; do
            base=$( basename $f )
            d=$( dirname $f )
            echo "$base" | grep -q ".vcf.bgz$"
            r=$?
            if [ $r -eq 0 ] ; then
                newName=$( echo $base | sed 's#.vcf.bgz$#.vcf.gz#g' )
                mv $f ${dollar}{d}/${dollar}{newName}
                fileListArgs="${dollar}{fileListArgs} -I ${dollar}{d}/${dollar}{newName}"
            else
                fileListArgs="${dollar}{fileListArgs} -I $f"
            fi
        done

        echo "Using file list: ${dollar}{fileListArgs}"

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

         gatk --java-options "-Xmx${command_mem}m -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            MergeVcfs \
             ${dollar}{fileListArgs} \
             -O ${output_vcf_file}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}
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
         File vcf_file = "${output_vcf_file}"
         File timing_info = timing_output_file
     }
 }
