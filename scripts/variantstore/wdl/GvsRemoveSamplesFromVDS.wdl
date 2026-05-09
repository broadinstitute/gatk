version 1.0

import "GvsUtils.wdl" as Utils


workflow GvsRemoveSamplesFromVDS {
    meta {
        description: "Removes a specified set of samples from a Hail VDS and writes the result to a new path. Implemented by remove_samples_from_vds.py, which: (1) validates the removal file for duplicates and checks that at least one listed research ID exists in the VDS; (2) filters out the specified samples and removes dead alleles; (3) drops monomorphic reference rows in the variant data left behind after sample removal (important when a withdrawn sample was the sole carrier of an alt allele at a site). GVS-produced VDSes do not contain a GT field, so no GT recalculation is needed or performed. Phasing information is passed through unchanged."
    }
    input {
        Boolean use_tiny_dataproc_cluster = false
        String input_vds_path
        String samples_to_remove_path
        String output_vds_path

        String cluster_prefix = "vds-cluster"
        String? hail_temp_path
        String region = "us-central1"

        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Boolean leave_cluster_running_at_end = false
        Float? master_memory_fraction

        String? git_branch_or_tag
        String? hail_version
        File? hail_wheel
        String? basic_docker
        String? variants_docker
        String? workspace_bucket
        String? workspace_project
        String? cloud_sdk_slim_docker
    }

    parameter_meta {
        use_tiny_dataproc_cluster: {
            help: "If true, use a small Dataproc autoscaling configuration suited for integration tests. Defaults to false (large configuration for production callsets)."
        }
        input_vds_path: {
            help: "Path to the input VDS from which samples will be removed."
        }
        samples_to_remove_path: {
            help: "Path to a single-column file listing research IDs of samples to remove. Must have a header of 'research_id', with one research ID per line."
        }
        output_vds_path: {
            help: "Path to write the output VDS with the specified samples removed."
        }
        cluster_prefix: {
            help: "Prefix of the Dataproc cluster name."
        }
        hail_temp_path: {
            help: "Optional Hail temp path for intermediate files."
        }
        region: {
            help: "GCP region, e.g. us-central1."
        }
        hail_version: {
            help: "Optional Hail version. Cannot define both this parameter and `hail_wheel`."
        }
        hail_wheel: {
            help: "Optional Hail wheel. Cannot define both this parameter and `hail_version`."
        }
    }

    if (!defined(variants_docker) || !defined(basic_docker) || !defined(cloud_sdk_slim_docker) || !defined(workspace_bucket) || !defined(workspace_project) || !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
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

    call RemoveSamplesFromVDS {
        input:
            prefix = cluster_prefix,
            use_tiny_dataproc_cluster = use_tiny_dataproc_cluster,
            input_vds_path = input_vds_path,
            samples_to_remove_path = samples_to_remove_path,
            output_vds_path = output_vds_path,
            hail_version = effective_hail_version,
            hail_wheel = hail_wheel,
            hail_temp_path = hail_temp_path,
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
            remove_samples_from_vds_script = GetHailScripts.remove_samples_from_vds_script,
            workspace_project = effective_google_project,
            region = region,
            workspace_bucket = effective_workspace_bucket,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            leave_cluster_running_at_end = leave_cluster_running_at_end,
            cluster_max_idle_minutes = cluster_max_idle_minutes,
            cluster_max_age_minutes = cluster_max_age_minutes,
            master_memory_fraction = master_memory_fraction,
    }

    output {
        String cluster_name = RemoveSamplesFromVDS.cluster_name
        Boolean done = true
    }
}

task RemoveSamplesFromVDS {
    input {
        String prefix
        Boolean use_tiny_dataproc_cluster
        String input_vds_path
        String samples_to_remove_path
        String output_vds_path
        Boolean leave_cluster_running_at_end
        File remove_samples_from_vds_script
        File run_in_hail_cluster_script
        String? hail_version
        File? hail_wheel
        String? hail_temp_path
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction

        String workspace_project
        String workspace_bucket
        String region

        String cloud_sdk_slim_docker
    }

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

        if [[ -z "~{hail_temp_path}" ]]
        then
            hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"
        else
            hail_temp_path="~{hail_temp_path}"
        fi


        # Construct a JSON of arguments for the python script to be run in the Hail cluster.
        cat > script-arguments.json <<FIN
        {
            "input-vds-path": "~{input_vds_path}",
            "samples-to-remove-path": "~{samples_to_remove_path}",
            "output-vds-path": "~{output_vds_path}",
            "temp-path": "${hail_temp_path}"
        }
        FIN

        # Run the Hail python script to remove samples from the VDS.
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{remove_samples_from_vds_script} \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            ~{true='--use-tiny-dataproc-cluster' false='' use_tiny_dataproc_cluster} \
            --region ~{region} \
            --workspace-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end}
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
        Boolean done = true
    }
}
