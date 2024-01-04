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

    call Utils.GetWorkspaceInfo {
        input:
            variants_docker = GetToolVersions.variants_docker,
            workspace_id = GetToolVersions.workspace_id,
    }

    if (!defined(submission_id_to_search)) {
        call FindAvroExtractDirectories {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                variants_docker = GetToolVersions.variants_docker,
                perform_deletion = perform_deletion,
        }

        call FindHailTempDirectories {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_namespace = GetWorkspaceInfo.workspace_namespace,
                workspace_name = GetWorkspaceInfo.workspace_name,
                variants_docker = GetToolVersions.variants_docker,
                perform_deletion = perform_deletion,
        }
    }

    if (defined(submission_id_to_search)) {
        call FindAvroExtractDirectories as FindAvroExtractDirectoriesInSpecifiedSubmission {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                submission_id = select_first([submission_id_to_search]),
                variants_docker = GetToolVersions.variants_docker,
                perform_deletion = perform_deletion,
        }

        call FindHailTempDirectories as FindHailTempDirectoriesInSpecifiedSubmission {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_namespace = GetWorkspaceInfo.workspace_namespace,
                workspace_name = GetWorkspaceInfo.workspace_name,
                variants_docker = GetToolVersions.variants_docker,
                submission_id = select_first([submission_id_to_search]),
                perform_deletion = perform_deletion,
        }
    }

    output {
        Array[String] nonempty_avro_directories = select_first([FindAvroExtractDirectories.nonempty_avro_directories, FindAvroExtractDirectoriesInSpecifiedSubmission.nonempty_avro_directories])
        Array[String] nonempty_temp_dirs = select_first([FindHailTempDirectories.nonempty_temp_dirs, FindHailTempDirectoriesInSpecifiedSubmission.nonempty_temp_dirs])
    }
}


task FindHailTempDirectories {
    input {
        String variants_docker
        String workspace_bucket
        String workspace_namespace
        String workspace_name
        String? submission_id
        Boolean perform_deletion
    }
    meta {
        # Don't cache this, we might clean up the directory and then re-run a configuration that repopulates it.
        volatile: true
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # This gets paths like <workspace bucket>/submissions/<submission id>/GvsCreateVDS/<workflow id>'.
        # We will need both the submission id and workflow id when we call into FISS.
        gsutil ls '~{workspace_bucket}/submissions/*/GvsCreateVDS/' > create_vds_raw.txt

        # Emit submission / workflow id pairs delimited by a space
        if [[ -n "~{submission_id}" ]]
        then
            sed -n -E 's!.*/submissions/(~{submission_id})/([-a-f0-9]+)/$!\1 \2!p' create_vds_raw.txt > pairs.txt
        else
            sed -n -E 's!.*/submissions/([-a-f0-9]+).*/([-a-f0-9]+)/$!\1 \2!p' create_vds_raw.txt > pairs.txt
        fi

        # Iterate over submission / workflow id pairs, fetch workflow metadata and extracting hail temp dir if present.
        while read -r submission workflow
        do
            python -c "from firecloud import fiss; resp = fiss.fapi.get_workflow_metadata('~{workspace_namespace}', '~{workspace_name}', '${submission}', '${workflow}'); print(resp.text)" > temp.json
            jq -M -r '.. | .inputs?.hail_temp_path? //empty' temp.json >> temp_dirs.txt
        done < pairs.txt

        # Uniquify and delete blank lines
        sort -u temp_dirs.txt | sed '/^[[:space:]]*$/d' > unique_temp_dirs.txt

        if [[ "~{perform_deletion}" == "true" ]]
        then
            for dir in $(cat unique_temp_dirs.txt)
            do
                gsutil rm -rf "${dir}"
            done
        fi
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        File pairs = "pairs.txt"
        File temp_dirs = "temp_dirs.txt"
        File unique_temp_dirs = "unique_temp_dirs.txt"
        Array[String] nonempty_temp_dirs = read_lines("unique_temp_dirs.txt")
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
            gsutil ls '~{workspace_bucket}/submissions/*/GvsExtractAvroFilesForHail/*/call-OutputPath/avro/**/*.avro' > avros.txt 2> err.txt
        fi

        if [[ $? -eq 1 ]]
        then
            if !grep "CommandException: One or more URLs matched no objects" err.txt
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
