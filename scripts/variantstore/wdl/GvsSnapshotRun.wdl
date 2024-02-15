version 1.0

import "GvsUtils.wdl" as Utils
 #
workflow GvsSnapshotTest {
    input {
        String? git_branch_or_tag
        String? cloud_sdk_docker
    }


#    Array[String] table_patterns

#    call Utils.GetToolVersions {
#        input:
#            git_branch_or_tag = git_branch_or_tag,
#    }
#
#    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])


    call Utils.SnapshotTables {
        input:
            project_id = "gvs-internal",
            dataset_name = "hatcher_microsoft_walkthrough_v2",
            snapshot_dataset = "hatcher_vs_299_test_storage",
            retrieval_key = "test_run",
            snapshot = true,
            cloud_sdk_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:435.0.0-alpine",
    }
}