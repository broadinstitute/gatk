version 1.0

# This WDL will create a VDS in Hail running in a Dataproc cluster.
import "GvsUtils.wdl" as Utils


workflow GvsCreateVDS {
    input {
        String? git_branch_or_tag
        String vds_destination_path
        String avro_path
        Boolean use_classic_VQSR = false
        Boolean leave_cluster_running_at_end = false
        Int? intermediate_resume_point
        String? hail_version
        String? hail_temp_path
        String cluster_prefix = "vds-cluster"
        String gcs_subnetwork_name = "subnetwork"
        String region = "us-central1"
        String? gcs_project
        String? workspace_bucket
        String? variants_docker
    }
    parameter_meta {
        # Analysis parameters, i.e., parameters that go to the Hail python code (submission_script below)
        avro_path : {
            help: "Input location for the avro files"
        }
        vds_output_path: {
            help: "Location for the final created VDS"
        }
        # Cluster parameters
        cluster_prefix: {
            help: "Prefix of the Dataproc cluster name"
        }
        gcs_subnetwork_name: {
            help: "Set to 'subnetwork' if running in Terra Cromwell"
        }
        region: {
            help: "us-central1"
        }
        hail_temp_path: {
            help: "Hail temp path to use, specify if resuming from a run that failed midway through creating intermediate VDSes."
        }
        intermediate_resume_point: {
            help: "Index at which to resume creating intermediate VDSes."
        }
    }


    if (!defined(variants_docker) || !defined(workspace_bucket) || !defined(gcs_project) || !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
    String effective_workspace_bucket = select_first([workspace_bucket, GetToolVersions.workspace_bucket])
    String effective_google_project = select_first([gcs_project, GetToolVersions.google_project])
    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    if (defined(intermediate_resume_point) && !defined(hail_temp_path)) {
        call Utils.TerminateWorkflow as NeedHailTempPath {
            input:
                message = "GvsCreateVDS called with an intermediate resume point but no specified hail temp path from which to resume",
                basic_docker = effective_variants_docker, # intentionally mismatched as basic was not already defined and this should be a superset
        }
    }

    if (!defined(intermediate_resume_point) && defined(hail_temp_path)) {
        call Utils.TerminateWorkflow as NeedIntermediateResumePoint {
            input:
                message = "GvsCreateVDS called with no intermediate resume point but a specified hail temp path, which isn't a known use case",
                basic_docker = effective_variants_docker, # intentionally mismatched as basic was not already defined and this should be a superset
        }
    }

    call create_vds {
        input:
            prefix = cluster_prefix,
            vds_path = vds_destination_path,
            avro_path = avro_path,
            use_classic_VQSR = use_classic_VQSR,
            hail_version = effective_hail_version,
            hail_temp_path = hail_temp_path,
            intermediate_resume_point = intermediate_resume_point,
            gcs_project = effective_google_project,
            region = region,
            workspace_bucket = effective_workspace_bucket,
            gcs_subnetwork_name = gcs_subnetwork_name,
            variants_docker = effective_variants_docker,
            leave_cluster_running_at_end = leave_cluster_running_at_end,
    }

    call validate_vds {
        input:
            go = select_first([create_vds.done, true]),
            prefix = cluster_prefix,
            vds_path = vds_destination_path,
            hail_version = effective_hail_version,
            gcs_project = effective_google_project,
            region = region,
            workspace_bucket = effective_workspace_bucket,
            gcs_subnetwork_name = gcs_subnetwork_name,
            variants_docker = effective_variants_docker,
    }

    output {
        String cluster_name = create_vds.cluster_name
        Boolean done = true
    }

}

task create_vds {
    input {
        String prefix
        String vds_path
        String avro_path
        Boolean use_classic_VQSR
        Boolean leave_cluster_running_at_end
        String? hail_version
        String? hail_temp_path
        Int? intermediate_resume_point

        String gcs_project
        String workspace_bucket
        String region
        String gcs_subnetwork_name

        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        pip3 install --upgrade pip
        pip3 install hail~{'==' + hail_version}
        pip3 install --upgrade google-cloud-dataproc

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | xxd -p)-$(head -c2 < /dev/urandom | xxd -p)"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        if [[ -z "~{hail_temp_path}" ]]
        then
            hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"
        else
            hail_temp_path="~{hail_temp_path}"
        fi

        # Set up the autoscaling policy
        cat > auto-scale-policy.yaml <<FIN
        workerConfig:
            minInstances: 2
            maxInstances: 2
        secondaryWorkerConfig:
            maxInstances: 500
        basicAlgorithm:
            cooldownPeriod: 120s
            yarnConfig:
                scaleUpFactor: 1.0
                scaleDownFactor: 1.0
                gracefulDecommissionTimeout: 120s
        FIN
        gcloud dataproc autoscaling-policies import gvs-autoscaling-policy --project=~{gcs_project} --source=auto-scale-policy.yaml --region=~{region} --quiet

        # Run the hail python script to make a VDS
        python3 /app/run_in_hail_cluster.py \
            --script-path /app/hail_gvs_import.py \
            --secondary-script-path /app/import_gvs.py \
            --account ${account_name} \
            --autoscaling-policy gvs-autoscaling-policy \
            --region ~{region} \
            --gcs-project ~{gcs_project} \
            --cluster-name ${cluster_name} \
            --avro-path ~{avro_path} \
            --vds-path ~{vds_path} \
            --temp-path ${hail_temp_path} \
            ~{'--intermediate-resume-point ' + intermediate_resume_point} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end} \
            ~{true='--use-classic-vqsr' false='' use_classic_VQSR}
    >>>

    runtime {
        memory: "6.5 GB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: 0
        docker: variants_docker
        bootDiskSizeGb: 10
    }

    output {
        String cluster_name = read_string("cluster_name.txt")
        Boolean done = true
    }
}

task validate_vds {
    input {
        Boolean go
        String prefix
        String vds_path
        String? hail_version
        String gcs_project
        String workspace_bucket
        String region
        String gcs_subnetwork_name
        String variants_docker
    }

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        pip3 install --upgrade pip
        pip3 install hail~{'==' + hail_version}
        pip3 install --upgrade google-cloud-dataproc

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | xxd -p)-$(head -c2 < /dev/urandom | xxd -p)"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt
        hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"

        # The autoscaling policy gvs-autoscaling-policy will exist already from the VDS creation

        # Run the hail python script to validate a VDS
        python3 /app/run_in_hail_cluster.py \
        --script-path /app/vds_validation.py \
        --account ${account_name} \
        --autoscaling-policy gvs-autoscaling-policy \
        --region ~{region} \
        --gcs-project ~{gcs_project} \
        --cluster-name ${cluster_name} \
        --vds-path ~{vds_path} \
        --temp-path ${hail_temp_path}
    >>>

    runtime {
        memory: "6.5 GB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: 0
        docker: variants_docker
        bootDiskSizeGb: 10
    }
    output {
        String cluster_name = read_string("cluster_name.txt")
    }
}
