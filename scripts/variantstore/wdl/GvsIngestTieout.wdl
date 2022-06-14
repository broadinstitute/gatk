version 1.0

import "GvsAssignIds.wdl" as GvsAssignIds
import "GvsImportGenomes.wdl" as GvsImportGenomes
import "GvsUtils.wdl" as Utils

workflow GvsIngestTieout {
    input {
        String project
        String reference_dataset_name
        String branch_name
        Array[String] sample_names
        Array[File] input_vcfs
        Array[File] input_vcf_indexes
        String? service_account_json_path
    }

    call Utils.BuildGATKJarAndCreateDataset {
        input:
            branch_name = branch_name,
            dataset_prefix = "ingest_tieout"
    }

    call GvsAssignIds.GvsAssignIds  {
        input:
            project_id = project,
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            external_sample_names = sample_names,
            assign_ids_gatk_override = BuildGATKJarAndCreateDataset.jar,
            service_account_json_path = service_account_json_path
    }

    call GvsImportGenomes.GvsImportGenomes {
        input:
            go = GvsAssignIds.done,
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            project_id = project,
            external_sample_names = sample_names,
            input_vcfs = input_vcfs,
            input_vcf_indexes = input_vcf_indexes,
            load_data_gatk_override = BuildGATKJarAndCreateDataset.jar,
            service_account_json_path = service_account_json_path
    }

    call IngestTieout {
        input:
            dataset_name = BuildGATKJarAndCreateDataset.dataset_name,
            reference_dataset_name = reference_dataset_name,
            project = project,
            stderrs = GvsImportGenomes.load_data_stderrs
    }
}


task IngestTieout {
    input {
        Boolean? go
        String dataset_name
        String reference_dataset_name
        String project
        Array[File] stderrs
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

            bq query --location=US --project_id=~{project} --format=csv --use_legacy_sql=false \
                "select actual.sample_id, expected.sample_id from
                (select sample_id, count(*) as count from \`gvs-internal.~{dataset_name}.${table_name}\` group by sample_id) actual full outer join
                (select sample_id, count(*) as count from \`gvs-internal.~{reference_dataset_name}.${table_name}\` group by sample_id) expected on actual.sample_id = expected.sample_id
                where actual.count != expected.count OR actual.sample_id is null OR expected.sample_id is null" > differences.txt

            if [[ -s differences.txt ]]; then
                fail=1
                echo "${table_name} row counts are mismatched for the following samples:"
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
        docker: "gcr.io/google.com/cloudsdktool/cloud-sdk:latest"
        memory: "14 GB"
        disks: "local-disk 2000 HDD"
        preemptible: 3
        cpu: 4
    }
}
