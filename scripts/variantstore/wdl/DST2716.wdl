version 1.0

import "GvsUtils.wdl" as Utils

workflow DST2716 {
    input {
        Boolean leave_cluster_running_at_end = false
        String region = "us-central1"
        String gcs_subnetwork_name = "subnetwork"
        String cluster_prefix = "vds-cluster"
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction
        String vds_path
        String chrom
        Int position
    }

    call Utils.GetToolVersions

    call Utils.GetHailScripts {
        input:
            variants_docker = GetToolVersions.variants_docker,
    }

    call InspectVds {
        input:
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
            site_specific_variant_qc_script = GetHailScripts.site_specific_variant_qc_script,
            prefix = cluster_prefix,
            vds_path = vds_path,
            chrom = chrom,
            position = position,
            hail_version = GetToolVersions.hail_version,
            workspace_project = GetToolVersions.google_project,
            region = region,
            workspace_bucket = GetToolVersions.workspace_bucket,
            gcs_subnetwork_name = gcs_subnetwork_name,
            leave_cluster_running_at_end = leave_cluster_running_at_end,
            cluster_max_idle_minutes = cluster_max_idle_minutes,
            cluster_max_age_minutes = cluster_max_age_minutes,
            master_memory_fraction = master_memory_fraction,
            cloud_sdk_slim_docker = GetToolVersions.cloud_sdk_slim_docker,
    }
}



task InspectVds {
    input {
        Boolean go = true
        File run_in_hail_cluster_script
        File site_specific_variant_qc_script
        String prefix
        String vds_path
        String chrom
        Int position
        String? hail_version
        File? hail_wheel
        String? hail_temp_path
        String workspace_project
        String workspace_bucket
        String region
        String gcs_subnetwork_name
        String cloud_sdk_slim_docker
        Boolean leave_cluster_running_at_end
        Int? cluster_max_idle_minutes
        Int? cluster_max_age_minutes
        Float? master_memory_fraction
    }

    meta {
        # should always be run
        volatile: true
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
            "chrom": "~{chrom}",
            "position": ~{position},
            "temp-path": "${hail_temp_path}"
        }
        FIN

        # Run the hail python script to validate a VDS
        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{site_specific_variant_qc_script} \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            --autoscaling-policy gvs-autoscaling-policy \
            --region ~{region} \
            --gcs-project ~{workspace_project} \
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
    }
}