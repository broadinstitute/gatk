version 1.0

workflow WorkspaceId {
    call WorkspaceIdTask

    output {
        String workspace_id = WorkspaceIdTask.workspace_id
    }
}

task WorkspaceIdTask {
    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        # Sniff the workspace bucket out of the delocalization script and extract the workspace id from that.
        sed -n -E 's!.*gs://fc-(secure-)?([^\/]+).*!\2!p' /cromwell_root/gcs_delocalization.sh | sort -u > workspace_id.txt
    >>>

    runtime {
        docker: "ubuntu:latest"
    }

    output {
        String workspace_id = read_string("workspace_id.txt")
    }
}
