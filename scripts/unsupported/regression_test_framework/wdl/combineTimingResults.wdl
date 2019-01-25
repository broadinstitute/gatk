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
    Array[File]   inputTimingCSV
    Array[String] inputTimingResultTaskNames

    # ------------------------------------------------
    # Call our tasks:
    call CombineTimingTask {
        input:
            timing_files          = inputTimingResults,
            timing_csvs           = inputTimingCSV,
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
    Array[File] timing_csvs
    Array[String] timing_names

    String? gatk_docker = "broadinstitute/gatk:gatkbase-2.0.2"

    String output_name = "combined.txt"
    String output_csv = "combined.csv"

    Int num_files = length(timing_files)

    String dollar = "$"

    ####################################################################################
    # Do the work:
    command <<<
        echo "File summary:\n" >> ${output_name}

        FILES=(${sep=" " timing_files} )
        CSVS=(${sep=" " timing_csvs} )
        NAMES=(${sep=" " timing_names} )

        echo "Type,Title,Runtime" > ${output_csv}

        for ((i=0; i<${num_files}; i++)); do
            cat ${dollar}{CSVS[i]} >> ${output_csv}
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
        File combined_csv         = output_csv
        File combined_file        = output_name
    }

}

task PlotRuntimeTimingResults {

    ####################################################################################
    # Inputs:
    File timing_csv
    String bam_name

    String? gatk_docker = "broadinstitute/gatk:gatkbase-2.0.2"

    String boxplot = "plotboxplot.pdf"
    String histogram = "plothisogram.pdf"

    String dollar = "$"

    ####################################################################################
    # Do the work:
    command <<<
        gsutil cp gs://emeryj-testing/generatePerformanceGraphs.R .
        Rscript generatePerformanceGraphs.R ${timing_csv} bamname
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
        File combined_csv         = boxplot
        File combined_file        = histogram
    }

}