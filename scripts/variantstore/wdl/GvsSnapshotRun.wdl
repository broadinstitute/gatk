version 1.0

import "GvsUtils.wdl" as Utils
 #
workflow GvsSnapshotTest {
    input {
        String? git_branch_or_tag
        String? cloud_sdk_docker
    }

#    String project_id
#    String dataset_name
#    String snapshot_dataset
#    String storage_key
#    Array[String] table_patterns
#    Boolean snapshot = true
#    String cloud_sdk_docker

    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])


    call Utils.SnapshotTables {
        input:
            project_id = "gvs-dev",
            dataset_name = "hatcher_microsoft_walkthrough",
            cloud_sdk_docker = effective_cloud_sdk_docker,

    }
}