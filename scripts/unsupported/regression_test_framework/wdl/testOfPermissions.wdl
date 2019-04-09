workflow test_of_permissions {
    File script

    call testScript as test {
      input:
        script_location = script
    }
}


# Comparison of mark duplicates marked reads between files
task testScript {
    File script_location

    command {
        bash ${script_location} gs://haplotypecallerspark-evaluation/inputData/NexPond-359781.bam gs://haplotypecallerspark-evaluation/testinput/NIST.interval_list gs://haplotypecallerspark-evaluation/testoutput/test1.bam
    }
    runtime {
        docker:  "docker.io/jamesemery/gatk-nightly:markDuplicatesOpticalFix8"
    }
}