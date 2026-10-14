version 1.0

# Runs an arbitrary Hail Python script on an autoscaling Dataproc cluster via run_in_hail_cluster.py. For one-off
# analyses that do not merit a dedicated workflow; anything run routinely should get its own WDL.
import "GvsUtils.wdl" as Utils

workflow GvsRunHailScript {
    input {
        File script
        # Shipped with the job via `--py-files`, so the main script can import them as modules.
        Array[File] secondary_scripts = []
        # Each entry is rendered as `--key value`.
        Map[String, String] script_arguments

        String cluster_prefix = "gvs-hail-script"
        String region = "us-central1"
        # Compute Engine stockouts are per-zone, and Dataproc's default placement picks one zone and fails rather than
        # moving on. "auto" tries every zone in the region.
        String cluster_zones = "auto"
        Int num_primary_workers = 2
        Int max_secondary_workers = 200
        Int num_local_ssds = 1
        String worker_machine_type = "n1-highmem-8"
        String master_machine_type = "n1-highmem-32"
        Boolean use_tiny_dataproc_cluster = false

        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Boolean leave_cluster_running_at_end = false
        Float? master_memory_fraction

        String? git_branch_or_tag
        String? hail_version
        File? hail_wheel
        String? basic_docker
        String? variants_docker
        String? workspace_project
        String? cloud_sdk_slim_docker
    }

    parameter_meta {
        script: {
            help: "GCS path to the Hail Python script to run."
        }
        secondary_scripts: {
            help: "GCS paths to Python modules the script imports. They are shipped as given, not drawn from the Variants image."
        }
        script_arguments: {
            help: "Arguments for the script, e.g. {\"mode\": \"counts\", \"output\": \"gs://...\"}. Each entry becomes `--key value`, so keys are written without the leading dashes, every key needs a value, a key cannot repeat, and values must not contain spaces."
        }
        cluster_zones: {
            help: "Comma-separated zones to try in order, or auto for every zone in the region."
        }
        num_primary_workers: {
            help: "Non-preemptible workers, present for the life of the cluster. Raise this for a job that cannot tolerate preemption; otherwise leave it low and let the autoscaler add secondary workers."
        }
        max_secondary_workers: {
            help: "Ceiling on the preemptible secondary workers the autoscaling policy may add."
        }
        hail_version: {
            help: "Optional Hail version. Cannot define both this parameter and `hail_wheel`."
        }
        hail_wheel: {
            help: "Optional Hail wheel file. Cannot define both this parameter and `hail_version`."
        }
    }

    if (!defined(variants_docker) || !defined(basic_docker) || !defined(cloud_sdk_slim_docker) || !defined(workspace_project) || !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_google_project = select_first([workspace_project, GetToolVersions.google_project])
    String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])

    if (defined(hail_version) && defined(hail_wheel)) {
        call Utils.TerminateWorkflow as BothHailVersionAndHailWheelDefined {
            input:
                message = "Cannot define both `hail_version` and `hail_wheel`, exiting.",
                basic_docker = effective_basic_docker,
        }
    }

    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    call Utils.GetHailScripts {
        input:
            variants_docker = effective_variants_docker,
    }

    call RunHailScript {
        input:
            script = script,
            secondary_scripts = secondary_scripts,
            script_arguments = script_arguments,
            prefix = cluster_prefix,
            region = region,
            cluster_zones = cluster_zones,
            num_primary_workers = num_primary_workers,
            max_secondary_workers = max_secondary_workers,
            num_local_ssds = num_local_ssds,
            worker_machine_type = worker_machine_type,
            master_machine_type = master_machine_type,
            use_tiny_dataproc_cluster = use_tiny_dataproc_cluster,
            cluster_max_idle_minutes = cluster_max_idle_minutes,
            cluster_max_age_minutes = cluster_max_age_minutes,
            leave_cluster_running_at_end = leave_cluster_running_at_end,
            master_memory_fraction = master_memory_fraction,
            hail_version = effective_hail_version,
            hail_wheel = hail_wheel,
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
            workspace_project = effective_google_project,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
    }

    output {
        String cluster_name = RunHailScript.cluster_name
        File log = RunHailScript.log
    }
}

task RunHailScript {
    input {
        File script
        Array[File] secondary_scripts
        Map[String, String] script_arguments

        String prefix
        String region
        String cluster_zones
        Int num_primary_workers
        Int max_secondary_workers
        Int num_local_ssds
        String worker_machine_type
        String master_machine_type
        Boolean use_tiny_dataproc_cluster
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Boolean leave_cluster_running_at_end
        Float? master_memory_fraction

        String? hail_version
        File? hail_wheel
        File run_in_hail_cluster_script
        String workspace_project
        String cloud_sdk_slim_docker
    }

    meta {
        # Never cache: the script's outputs are written by the script itself, outside Cromwell's view.
        volatile: true
    }

    File script_arguments_json = write_json(script_arguments)

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        apt-get update
        apt-get install --assume-yes python3.11-venv
        python3 -m venv ./localvenv
        . ./localvenv/bin/activate

        pip3 install --upgrade pip

        if [[ ! -z "~{hail_wheel}" ]]
        then
            pip3 install ~{hail_wheel}
        else
            pip3 install hail~{'==' + hail_version}
        fi

        pip3 install --upgrade google-cloud-dataproc ijson

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | od -h -An | tr -d '[:space:]')-$(head -c2 < /dev/urandom | od -h -An | tr -d '[:space:]')"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        secondary_script_args=()
        for secondary_script in ~{sep=' ' secondary_scripts}
        do
            secondary_script_args+=(--secondary-script-path-list "${secondary_script}")
        done

        cat ~{script_arguments_json}

        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{script} \
            ${secondary_script_args[@]+"${secondary_script_args[@]}"} \
            --script-arguments-json-path ~{script_arguments_json} \
            --account ${account_name} \
            ~{true='--use-tiny-dataproc-cluster' false='' use_tiny_dataproc_cluster} \
            --num-primary-workers ~{num_primary_workers} \
            --max-secondary-workers ~{max_secondary_workers} \
            --region ~{region} \
            --zones ~{cluster_zones} \
            --num-local-ssds ~{num_local_ssds} \
            --worker-machine-type ~{worker_machine_type} \
            --master-machine-type ~{master_machine_type} \
            --workspace-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end} \
            2>&1 | tee script.log
    >>>

    runtime {
        memory: "6.5 GB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: 0
        docker: cloud_sdk_slim_docker
        bootDiskSizeGb: 10
    }

    output {
        String cluster_name = read_string("cluster_name.txt")
        File log = "script.log"
    }
}
