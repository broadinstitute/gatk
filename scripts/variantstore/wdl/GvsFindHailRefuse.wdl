version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsFindHailRefuse {
    input {
        String? submission_id_to_search
        Boolean perform_deletion = false
    }
    parameter_meta {
        submission_id_to_search: "If defined, the submission id under which to search for GvsExtractAvroFilesForHail and GvsCreateVDS workflow runs, useful for testing against integration test runs. If undefined, search all submissions in the current workspace, which will only find top-level invocations of GvsExtractAvroFilesForHail and GvsCreateVDS and not invocations as subworkflows of GvsQuickstartIntegration."
    }

    call Utils.GetToolVersions

    call FindAvroExtractDirectories {
        input:
            workspace_bucket = GetToolVersions.workspace_bucket,
            variants_docker = GetToolVersions.variants_docker,
            submission_id = submission_id_to_search,
            perform_deletion = perform_deletion,
    }

    call FindHailTempDirectories {
        input:
            workspace_bucket = GetToolVersions.workspace_bucket,
            variants_docker = GetToolVersions.variants_docker,
            perform_deletion = perform_deletion,
    }
    output {
        Array[String] nonempty_avro_directories = FindAvroExtractDirectories.nonempty_avro_directories
        Array[String] nonempty_temp_dirs = FindHailTempDirectories.nonempty_temp_dirs
    }
}


task FindAvroExtractDirectories {
    input {
        String variants_docker
        String workspace_bucket
        String? submission_id
        Boolean perform_deletion
    }
    meta {
        # Do not cache in case new runs are made in a workspace that was previously cleaned up.
        volatile: true
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Searching for Avros might fail and that's okay
        set +o errexit
        if [[ -n "~{submission_id}" ]]
        then
            gsutil ls '~{workspace_bucket}/submissions/~{submission_id}/**/GvsExtractAvroFilesForHail/*/call-OutputPath/avro/**/*.avro' > avros.txt 2> err.txt
        else
            gsutil ls '~{workspace_bucket}/submissions/**/GvsExtractAvroFilesForHail/*/call-OutputPath/avro/**/*.avro' > avros.txt 2> err.txt
        fi

        if [[ $? -eq 1 ]]
        then
            if ! grep "CommandException: One or more URLs matched no objects" err.txt
            then
                echo "gsutil ls failed for unknown reason, failing."
                exit 1
            fi
        fi

        set -o errexit

        # post process the above to get just the Avro directory paths
        sed -n -E 's!(.*/call-OutputPath/avro).*!\1!p' avros.txt | sort -u > avro_directories.txt

        if [[ "~{perform_deletion}" == "true" ]]
        then
            for avro_directory in $(cat avro_directories.txt)
            do
                gsutil rm -rf "${avro_directory}"
            done
        fi
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Array[String] nonempty_avro_directories = read_lines("avro_directories.txt")
        File err = "err.txt"
    }
}

task FindHailTempDirectories {
    input {
        String variants_docker
        String workspace_bucket
        Boolean perform_deletion
    }
    meta {
        # Don't cache this, we might clean up the directory and then run a workflow that repopulates it.
        volatile: true
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        set +o errexit
        gsutil ls "~{workspace_bucket}/hail-temp/hail-temp-*" > temp_dirs.txt 2> err.txt
        if [[ $? -eq 1 ]]
        then
            if ! grep "CommandException: One or more URLs matched no objects" err.txt
            then
                echo "gsutil ls failed for unknown reason, failing."
                exit 1
            fi
        fi
        set -o errexit

        if [[ "~{perform_deletion}" == "true" ]]
        then
            for avro_directory in $(cat temp_dirs.txt)
            do
                gsutil rm -rf "${avro_directory}"
            done
        fi
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Array[String] nonempty_temp_dirs = read_lines("temp_dirs.txt")
    }
}
