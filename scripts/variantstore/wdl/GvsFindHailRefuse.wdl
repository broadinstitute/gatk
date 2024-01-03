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

    if (!defined(submission_id_to_search)) {
        call FindAvroExtractDirectories {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                variants_docker = GetToolVersions.variants_docker,
                perform_deletion = perform_deletion,
        }

        call FindHailTempDirectoriese {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                workspace_namespace = GetToolVersions.workspace_namespace,
                workspace_name = GetToolVersions.workspace_name,
                variants_docker = GetToolVersions.variants_docker,
                perform_deletion = perform_deletion,
        }
    }

    if (defined(submission_id_to_search)) {
        call FindAvroExtractDirectoriesInSpecifiedSubmission {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                submission_id = select_first([submission_id_to_search]),
                variants_docker = GetToolVersions.variants_docker,
                perform_deletion = perform_deletion,
        }
    }


}


task FindHailTempDirectoriese {
    input {
        String variants_docker
        String workspace_bucket
        String workspace_namespace
        String workspace_name
        Boolean perform_deletion
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # This gets paths like <workspace bucket>/submissions/<submission id>/GvsCreateVDS/<workflow id>'.
        # We will need both the submission id and workflow id when we call into FISS.
        gsutil ls '~{workspace_bucket}/submissions/*/GvsCreateVDS/' > create_vds_raw.txt

        # Emit submission / workflow id pairs delimited by a space
        sed -n -E 's!.*/submissions/([-a-f0-9]+).*/([-a-f0-9]+)/$!\1 \2!p' create_vds_raw.txt > pairs.txt
        while read -r submission workflow
        do
            python -c "from firecloud import fiss; fiss.fapi.get_workflow_metadata(...)"
        done < pairs.rxt


    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Array[String] avro_paths = read_lines("avro_directories.txt")
    }
}


task FindAvroExtractDirectories {
    input {
        String variants_docker
        String workspace_bucket
        Boolean perform_deletion
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        gsutil ls '~{workspace_bucket}/submissions/*/GvsExtractAvroFilesForHail/*/call-OutputPath/avro/**/*.avro' > avros.txt

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
    }
}


task FindAvroExtractDirectoriesInSpecifiedSubmission {
    input {
        String variants_docker
        String workspace_bucket
        String submission_id
        Boolean perform_deletion
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        gsutil ls '~{workspace_bucket}/submissions/~{submission_id}/**/GvsExtractAvroFilesForHail/*/call-OutputPath/avro/**/*.avro' > avros.txt

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
    }
}
