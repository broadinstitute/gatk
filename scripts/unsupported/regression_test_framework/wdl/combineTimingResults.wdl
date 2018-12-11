# Description of inputs:
#
#   Required:
#     gatk_docker                    - GATK Docker image in which to run
#
#     inputTimingResutls             - Array of outputs for the test .
#     inputTimingResultTaskNames     - Base file name for generated output files.
#
# This WDL is used to combine the various outputs from a run of

workflow CombineTiming {

    # ------------------------------------------------
    # Input args:
    String gatk_docker = "broadinstitute/gatk:gatkbase-2.0.2"

    Array[File]   inputTimingResults
    Array[String] inputTimingResultTaskNames

    # ------------------------------------------------
    # Call our tasks:
    call CombineTimingTask {
        input:
            timing_files         = inputTimingResults,
            timing_names          = inputTimingResultTaskNames,
    }

    # ------------------------------------------------
    # Outputs:
    output {
        File summary_metrics = CombineTimingTask.combined_file
    }

}


task CombineTimingTask {

    ####################################################################################
    # Inputs:
    Array[File] timing_files
    Array[String] timing_names

    String? gatk_docker = "broadinstitute/gatk:gatkbase-2.0.2"

    String output_name = "combined.txt"

    Int num_files = length(timing_files)

    String dollar = "$"

    ####################################################################################
    # Do the work:
    command <<<
        echo "File summary:\n" >> ${output_name}

        FILES=(${sep=" " timing_files} )
        NAMES=(${sep=" " timing_names} )

        for ((i=0; i<${num_files}; i++)); do
            printf "%s\n$s\n" ${dollar}{NAMES[i]} `cat ${dollar}{FILES[i]}` >> ${output_name}
        done
    >>>

    ####################################################################################
    # Runtime Environment:
    runtime {
        cpu: 1
        disks: "local-disk 30 HDD"
        docker: "${gatk_docker}"
    }

    ####################################################################################
    # Outputs:
    output {
        File combined_file        = output_name
    }

}