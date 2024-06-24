version 1.0

import "GvsAssignIds.wdl" as GvsAssignIds
import "GvsImportGenomes.wdl" as GvsImportGenomes
import "GvsUtils.wdl" as Utils

workflow GvsIngestTieout {
    input {
        String project
        String reference_dataset_name
        String? git_branch_or_tag
        File sample_names
        File input_vcfs
        File input_vcf_indexes
        String? cloud_sdk_docker
    }

    # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
    # no calling WDLs that might supply `git_hash`).
    call Utils.GetToolVersions {
        input:
            git_branch_or_tag = git_branch_or_tag,
    }

    String effective_cloud_sdk_docker = select_first([cloud_sdk_docker, GetToolVersions.cloud_sdk_docker])

    call Utils.BuildGATKJarAndCreateDataset {
        input:
            git_branch_or_tag = git_branch_or_tag,
            dataset_prefix = "ingest_tieout"
    }

    call GvsAssignIds.GvsAssignIds  {
        input:
            project_id = project,
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            external_sample_names = sample_names,
    }

    call GvsImportGenomes.GvsImportGenomes {
        input:
            go = GvsAssignIds.done,
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            project_id = project,
            external_sample_names = sample_names,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
    }

    call IngestTieout {
        input:
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            reference_dataset_name = reference_dataset_name,
            project = project,
            stderrs = GvsImportGenomes.load_data_stderrs,
            cloud_sdk_docker = effective_cloud_sdk_docker,
    }
    output {
        String recorded_git_hash = GetToolVersions.git_hash
    }
}


task IngestTieout {
    input {
        Boolean? go
        String dataset_name
        String reference_dataset_name
        String project
        Array[File] stderrs
        String cloud_sdk_docker
    }
    meta {
        # Do not call cache, dataset may have been updated.
        volatile: true
    }

    parameter_meta {
        stderrs: {
            localization_optional: true
        }
    }

    command <<<
        set -o xtrace
        fail=0

        check_table() {
            local table_name=$1

            # bq query --max_rows check: ok, anything > 0 is an error, error message explicit about 100 row limit.
            bq --apilog=false query --project_id=~{project} --format=csv --use_legacy_sql=false \
                '(SELECT
                        sample_id, count(*) AS count
                    FROM
                        `gvs-internal.~{dataset_name}.'"${table_name}"'`
                    GROUP BY sample_id) actual
                FULL OUTER JOIN
                (SELECT
                        sample_id, count(*) AS count
                    FROM
                        `gvs-internal.~{reference_dataset_name}.'"${table_name}"'`
                    GROUP BY sample_id) expected
                ON actual.sample_id = expected.sample_id' > differences.txt

            if [[ -s differences.txt ]]; then
                fail=1
                echo "${table_name} row counts are mismatched for the following samples (100 samples max in output):"
                cat differences.txt
            fi
        }

        # This task is currently being called with a sample set of 2000 so it only needs to check a single vet and
        # ref_ranges table each. If this code were to be called with a sample set of more than 4000 it should test all
        # the additional vet and ref_ranges tables that would be introduced.
        if [[ ~{length(stderrs)} -gt 4000 ]]; then
            echo "IngestTieout invoked with a sample set of size ~{length(stderrs)} but is currently limited to sample sets no larger than 4000."
            exit 1
        fi

        check_table "ref_ranges_001"
        check_table "vet_001"

        mkdir logs
        cd logs

        # Get everything up to and including the /call-LoadData/ part of the GCS path.
        gcs_prefix=$(echo ~{stderrs[0]} | sed -E 's!(.*/call-LoadData/).*!\1!')

        # Recursively copy everything under this path.
        gsutil -m cp -r ${gcs_prefix} .

        # Find the stderr files in the call execution directories but not the ones in the pipelines-logs directories.
        find call-LoadData -name pipelines-logs -prune -o -name stderr -print > stderr_fofn.txt
        echo "Found $(wc -l stderr_fofn.txt | awk '{print $1}') stderr files in call directories"

        cat stderr_fofn.txt | xargs grep StatusRuntimeException | tee exceptions.txt
        if [[ $? -ne 0 ]]; then
            fail=1
            echo "Did not find any StatusRuntimeExceptions among the stderr files!"
        else
            grep UNAVAILABLE exceptions.txt
            if [[ $? -ne 0 ]]; then
                fail=1
                echo "No UNAVAILABLE StatusRuntimeExceptions found for run!"
            fi
            grep -E 'ABORTED|INTERNAL|CANCELLED' exceptions.txt
            if [[ $? -ne 0 ]]; then
                fail=1
                echo "No retryable StatusRuntimeExceptions (ABORTED, INTERNAL, CANCELLED) found for run!"
            fi
        fi

        cd ..

        if [[ $fail -ne 0 ]]; then
            exit 1;
        fi
    >>>

    output {
        File stderr_fofn = "logs/stderr_fofn.txt"
        File exceptions = "logs/exceptions.txt"
        Boolean done = true
    }

    runtime {
        docker: cloud_sdk_docker
        memory: "14 GB"
        disks: "local-disk 2000 HDD"
        preemptible: 3
        cpu: 4
    }
}
