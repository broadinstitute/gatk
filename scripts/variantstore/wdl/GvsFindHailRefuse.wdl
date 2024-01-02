version 1.0

import "GvsUtils.wdl" as Utils

workflow GvsFindHailRefuse {
    input {
        String? submission_id_to_search
    }
    parameter_meta {
        submission_id_to_search: "If defined, the submission id under which to search for GvsExtractAvroFilesForHail and GvsCreateVDS workflow runs, useful for testing against integration test runs. If undefined, search all submissions in the current workspace, which will only find top-level invocations of GvsExtractAvroFilesForHail and GvsCreateVDS and not invocations as subworkflows of GvsQuickstartIntegration."
    }

    call Utils.GetToolVersions

    if (!defined(submission_id_to_search)) {
        call FindAvroExtractDirectoriesInThisWorkspace {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                variants_docker = GetToolVersions.variants_docker,
        }
    }
    if (defined(submission_id_to_search)) {
        call FindAvroExtractDirectoriesInSpecifiedSubmission {
            input:
                workspace_bucket = GetToolVersions.workspace_bucket,
                submission_id = select_first([submission_id_to_search]),
                variants_docker = GetToolVersions.variants_docker,
        }
    }
}

task FindAvroExtractDirectoriesInThisWorkspace {
    input {
        String variants_docker
        String workspace_bucket
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        gsutil ls '~{workspace_bucket}/submissions/*/GvsExtractAvroFilesForHail/*/call-OutputPath/avro/**/*.avro' > avros.txt

        # post process the above to get just the Avro directory paths
        sed -n -E 's!(.*/call-OutputPath/avro).*!\1!p' avros.txt | sort -u > avro_directories.txt
    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Array[String] avro_paths = read_lines("avro_directories.txt")
    }
}


task FindAvroExtractDirectoriesInSpecifiedSubmission {
    input {
        String variants_docker
        String workspace_bucket
        String submission_id
    }
    command <<<
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        gsutil ls '~{workspace_bucket}/submissions/~{submission_id}/**/GvsExtractAvroFilesForHail/*/call-OutputPath/avro/**/*.avro' > avros.txt

        # post process the above to get just the Avro directory paths
        sed -n -E 's!(.*/call-OutputPath/avro).*!\1!p' avros.txt | sort -u > avro_directories.txt

    >>>
    runtime {
        docker: variants_docker
    }
    output {
        Array[String] avro_paths = read_lines("avro_directories.txt")
    }
}
