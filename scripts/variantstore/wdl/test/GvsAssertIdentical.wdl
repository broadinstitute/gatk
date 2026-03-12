version 1.0

workflow GvsAssertIdentical {
    input {
        String expected_output_prefix
        String expected_output_suffix
        Array[File] actual_vcfs
        String gatk_docker
    }

    call AssertIdenticalOutputs {
        input:
            expected_output_prefix = expected_output_prefix,
            expected_output_suffix = expected_output_suffix,
            actual_vcfs = actual_vcfs,
            gatk_docker = gatk_docker,
    }
}


task AssertIdenticalOutputs {
    input {
        String expected_output_prefix
        String expected_output_suffix
        Array[File] actual_vcfs
        String gatk_docker
    }
    parameter_meta {
        actual_vcfs: {
            localization_optional: true
        }
    }
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        failures=()

        # Where the current set of expected results lives in the cloud
        expected_prefix="~{expected_output_prefix}"
        # Remove a trailing slash if there is one
        expected_prefix=${expected_prefix%/}

        # Download all the expected data
        mkdir expected
        cd expected
        gcloud storage cp -r "${expected_prefix}"'/*.vcf~{expected_output_suffix}' .
        gzip -S ~{expected_output_suffix} -d *~{expected_output_suffix}
        cd ..

        mkdir actual
        cd actual
        touch actual_manifest.txt
        # Making the manifest is pretty uninteresting and very noisy so turn off xtrace temporarily.
        set +o xtrace
        for actual in ~{sep=' ' actual_vcfs}
        do
            echo $actual >> actual_manifest.txt
        done
        set -o xtrace

        cat actual_manifest.txt | gcloud storage cp -I .
        # Unzip actual result data.
        ls -1 | grep -E '\.vcf\~{expected_output_suffix}$' | xargs gzip -S ~{expected_output_suffix} -d
        cd ..

        echo "Header Check"
        # Headers first, these can yield useful diagnostics when there are mismatches.
        for vcf in $(ls -1 actual | grep -E '\.vcf$')
        do
          actual="actual/$vcf"
          expected="expected/$vcf"
          set +o errexit
          cmp <(grep '^#' $actual | grep -E -v '^##GATKCommandLine=') <(grep '^#' $expected | grep -E -v '^##GATKCommandLine=')
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
            # If there is a mismatch add it to a list of failures but keep on looking for mismatches.
            failures+=( $vcf )
          fi
        done

        echo "Header Failure Check"
        if [[ ${#failures[@]} -ne 0 ]]; then
          echo "Error: headers for the following files do not match:"
          for failure in ${failures[@]}; do
            echo $failure
            expected="expected/$failure"
            actual="actual/$failure"
            diff <(grep '^#' $actual) <(grep '^#' $expected)
          done
          exit 1
        fi

        echo "Overall Check"
        # If the headers all matched look for any mismatches in overall file content.
        fail=0
        for vcf in $(ls -1 actual | grep -E '\.vcf$')
        do
          expected="expected/$vcf"
          actual="actual/$vcf"
          set +o errexit
          cmp <(grep -E -v '^##GATKCommandLine=' $actual) <(grep -E -v '^##GATKCommandLine=' $expected)
          rc=$?
          set -o errexit
          if [[ $rc -ne 0 ]]; then
            echo "Error: file contents of expected and actual do not match: $vcf"
            fail=1
          fi
        done

        if [[ $fail -ne 0 ]]; then
          exit 1
        fi

        echo "All vcfs compared and matched!"
    >>>

    runtime {
        docker: gatk_docker
        disks: "local-disk 500 HDD"
    }

    output {
        Boolean done = true
    }
}
