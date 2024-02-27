version 1.0

import "GvsUtils.wdl" as Utils
 #
workflow GvsSnapshotTest {
    input {
        String run_name = "test_run"
        String retrieval_key = "post-extract"
        String target_restore_dataset
        String comparison_key
    }


#    Array[String] table_patterns

#    call Utils.GetToolVersions {
#        input:
#            git_branch_or_tag = git_branch_or_tag,
#    }
#
#    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])

#    call Utils.GetTableListFromDataset as TableNames {
#        input:
#        project_id = "gvs-internal",
#        dataset_name = "hatcher_microsoft_walkthrough_v2",
#        cloud_sdk_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:435.0.0-alpine",
#    }
#    # run it again, filtering out vet data for testing purposes
#    call Utils.GetTableListFromDataset as TableNamesFiltered {
#        input:
#        project_id = "gvs-internal",
#        dataset_name = "hatcher_microsoft_walkthrough_v2",
#        exclude_regex = "vet",
#        cloud_sdk_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:435.0.0-alpine",
#    }

#    call Utils.SnapshotTables as SnapshotTables {
#        input:
#            project_id = "gvs-internal",
#            dataset_name = "hatcher_microsoft_walkthrough_v2",
#            snapshot_dataset = "hatcher_vs_299_test_storage",
#            run_name = run_name,
#            retrieval_key = retrieval_key,
#            cloud_sdk_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:435.0.0-alpine",
#    }


#    call Utils.RestoreSnapshotForRun as RestoreSnapshot {
#        input:
#            project_id = "gvs-internal",
#            dest_dataset = target_restore_dataset,
#            snapshot_dataset = "hatcher_vs_299_test_storage",
#            run_name = run_name,
#            start_retrieval_key = start_retrieval_key,
#            end_retrieval_key = end_retrieval_key,
#            cloud_sdk_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:435.0.0-alpine",
#    }


    call Utils.CheckResultsAgainstStoredState as CheckResultsOfRun {
        input:
            project_id = "gvs-internal",
            run_dataset = target_restore_dataset,
            snapshot_dataset = "hatcher_vs_299_test_storage",
            run_name = run_name,
            comparison_key = comparison_key,
            cloud_sdk_docker = "gcr.io/google.com/cloudsdktool/cloud-sdk:435.0.0-alpine",
    }


}