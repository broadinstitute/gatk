version 1.0

# This WDL will create a VDS in Hail running in a Dataproc cluster.
import "GvsUtils.wdl" as Utils

workflow GvsCreateVDS {
    input {
        Boolean use_tiny_dataproc_cluster = false
        String avro_path
        String vds_destination_path

        String cluster_prefix = "vds-cluster"
        String? hail_temp_path
        Int? intermediate_resume_point
        String region = "us-central1"

        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction
        String worker_machine_type = "n2d-highmem-8"
        String master_machine_type = "n2d-highmem-32"
        Int? num_master_local_ssds
        # Compute Engine stockouts are per-zone, and Dataproc's default placement picks one zone and fails rather than
        # moving on. "auto" tries every zone in the region.
        String cluster_zones = "auto"

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
        avro_path : {
            help: "Input location for the avro files"
        }
        vds_destination_path: {
            help: "Location for the final created VDS"
        }
        cluster_prefix: {
            help: "Prefix of the Dataproc cluster name"
        }
        hail_temp_path: {
            help: "Hail temp path to use, specify if resuming from a run that failed midway through creating intermediate VDSes."
        }
        intermediate_resume_point: {
            help: "Index at which to resume creating intermediate VDSes."
        }
        region: {
            help: "us-central1"
        }
        hail_version: {
            help: "Optional Hail version, defaults to 0.2.130.post1. Cannot define both this parameter and `hail_wheel`."
        }
        hail_wheel: {
            help: "Optional Hail wheel. Cannot define both this parameter and `hail_version`."
        }
        worker_machine_type: {
            help: "Machine type for the Dataproc cluster's primary and secondary workers. Defaults to n2d-highmem-8. Must be a type listed in hailctl's MACHINE_MEM table for the Hail version in use, or cluster creation fails with a KeyError; as of Hail 0.2.134 that covers the n1, n2, n2d, e2, c2, m1 and m2 families but not C3, C3D or any 4th-generation family."
        }
        master_machine_type: {
            help: "Machine type for the Dataproc cluster's master. Defaults to n2d-highmem-32. Subject to the same hailctl MACHINE_MEM restriction as `worker_machine_type`. Not affected by `use_tiny_dataproc_cluster`, which only changes worker counts."
        }
        num_master_local_ssds: {
            help: "Local SSDs on the master, overriding the workers' 1. Valid counts depend on machine family and vCPU count, and larger machines need more: N2D machines with 32 or 48 vCPUs reject 1 and accept 2, 4, 8, 16 or 24, but an N2 machine with 32 vCPUs needs at least 4. Left unset, the master gets the workers' count, except that 1 becomes 2, which suits the default n2d-highmem-32 master and n1-highmem-32. For any other master type, check the allowed counts in the Compute Engine local SSD documentation and set this explicitly."
        }
        cluster_zones: {
            help: "Comma-separated zones to try in order when creating the Dataproc cluster, or auto for every zone in the region. A stockout is a per-zone condition, so a create that fails for capacity is retried in the next zone rather than failing the workflow. Set to an empty string to keep Dataproc's own placement, which picks one zone and does not retry."
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

    if (defined(intermediate_resume_point) && !defined(hail_temp_path)) {
        call Utils.TerminateWorkflow as NeedHailTempPath {
            input:
                message = "GvsCreateVDS called with an intermediate resume point but no specified hail temp path from which to resume",
                basic_docker = effective_basic_docker,
        }
    }

    if (!defined(intermediate_resume_point) && defined(hail_temp_path)) {
        call Utils.TerminateWorkflow as NeedIntermediateResumePoint {
            input:
                message = "GvsCreateVDS called with no intermediate resume point but a specified hail temp path, which isn't a known use case",
                basic_docker = effective_basic_docker,
        }
    }

    call Utils.GetHailScripts {
        input:
            variants_docker = effective_variants_docker,
    }

    call CreateVds {
        input:
            prefix = cluster_prefix,
            use_tiny_dataproc_cluster = use_tiny_dataproc_cluster,
            vds_path = vds_destination_path,
            avro_path = avro_path,
            hail_version = effective_hail_version,
            hail_wheel = hail_wheel,
            hail_temp_path = hail_temp_path,
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
            gvs_import_script = GetHailScripts.gvs_import_script,
            gvs_import_ploidy_script = GetHailScripts.gvs_import_ploidy_script,
            hail_gvs_import_script = GetHailScripts.hail_gvs_import_script,
            hail_gvs_util_script = GetHailScripts.hail_gvs_util_script,
            vds_validation_script = GetHailScripts.vds_validation_script,
            intermediate_resume_point = intermediate_resume_point,
            workspace_project = effective_google_project,
            region = region,
            workspace_bucket = effective_workspace_bucket,
            cloud_sdk_slim_docker = effective_cloud_sdk_slim_docker,
            cluster_max_idle_minutes = cluster_max_idle_minutes,
            cluster_max_age_minutes = cluster_max_age_minutes,
            master_memory_fraction = master_memory_fraction,
            worker_machine_type = worker_machine_type,
            master_machine_type = master_machine_type,
            num_master_local_ssds = num_master_local_ssds,
            cluster_zones = cluster_zones,
    }

    output {
        String cluster_name = CreateVds.cluster_name
        Boolean done = true
    }
}

task CreateVds {
    input {
        String prefix
        Boolean use_tiny_dataproc_cluster
        String vds_path
        String avro_path
        Boolean run_validation = true
        File run_in_hail_cluster_script
        File hail_gvs_import_script
        File hail_gvs_util_script
        File gvs_import_script
        File gvs_import_ploidy_script
        File vds_validation_script
        String? hail_version
        File? hail_wheel
        String? hail_temp_path
        Int? intermediate_resume_point
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction
        String worker_machine_type
        String master_machine_type
        Int? num_master_local_ssds
        String cluster_zones

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
        apt install --assume-yes python3.11-venv
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


        # construct a JSON of arguments for python script to be run in the hail cluster
        cat > script-arguments.json <<FIN
        {
            "vds-path": "~{vds_path}",
            "temp-path": "${hail_temp_path}",
            "avro-path": "~{avro_path}"
            ~{if (run_validation) then ', "run-validation": ""' else ''}
            ~{', "intermediate-resume-point": ' + intermediate_resume_point}
        }
        FIN

        # Run the hail python script to make a VDS
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{hail_gvs_import_script} \
            --secondary-script-path-list ~{gvs_import_script} \
            --secondary-script-path-list ~{hail_gvs_util_script} \
            --secondary-script-path-list ~{gvs_import_ploidy_script} \
            --secondary-script-path-list ~{vds_validation_script} \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            ~{true='--use-tiny-dataproc-cluster' false='' use_tiny_dataproc_cluster} \
            --region ~{region} \
            --worker-machine-type ~{worker_machine_type} \
            --master-machine-type ~{master_machine_type} \
            ~{'--num-master-local-ssds ' + num_master_local_ssds} \
            --zones '~{cluster_zones}' \
            --workspace-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction}
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
