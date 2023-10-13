version 1.0

# This WDL will create a VDS in Hail running in a Dataproc cluster.
import "GvsUtils.wdl" as Utils


workflow GvsCreateVDS {
    input {
        String vds_output_url
        String avro_path
        String? hail_version='0.2.124'

        String prefix = "vds-cluster"
        String gcs_subnetwork_name='subnetwork'
        String region = "us-central1"
    }
    parameter_meta {
        ## ANALYSIS PARAMETERS
        # i.e., parameters that go to the Hail python code (submission_script below)
        avro_path : {
            help: "Input location for the avro files"
        }
        vds_output_url: {
            help: "location for the final created VDS"
        }
        hail_version: {
            help: "0.2.124"
        }

        ## CLUSTER PARAMETERS
        prefix: {
            help: "used in construction of cluster name"
        }
        gcs_subnetwork_name: {
            help: "set to subnetwork if running in Terra Cromwell"
        }
        region: {
            help: "us-central1"
        }
    }

    call Utils.GetToolVersions

    call create_vds {
        input:
            prefix = prefix,
            vds_url = vds_output_url,
            avro_path = avro_path,
            hail_version = hail_version,
            gcs_project = GetToolVersions.google_project,
            workspace_bucket = GetToolVersions.workspace_bucket,
            gcs_subnetwork_name = gcs_subnetwork_name,
            variants_docker = GetToolVersions.variants_docker,
    }
}

task create_vds {
    input {
        String prefix
        String vds_url
        String avro_path
        String? hail_version

        # Cluster params
        String gcs_project
        String workspace_bucket
        String region = "us-central1"
        String gcs_subnetwork_name

        String variants_docker
    }

    String temp_path = "~{workspace_bucket}/quickstart-vds-for-wdl-tieout/temp-dir/"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)") ## do we need a service account anymore?

        pip3 install --upgrade pip
        pip3 install hail~{'==' + hail_version}
        pip3 install --upgrade google-cloud-dataproc

        # Generate a UUIDish random hex string of <8 hex chars (4 bytes)>-<4 hex chars (2 bytes)>
        hex="$(head -c4 < /dev/urandom | xxd -p)-$(head -c2 < /dev/urandom | xxd -p)"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        # Set up the autoscaling policy
        cat > auto-scale-policy.yaml <<FIN
        workerConfig:
            minInstances: 2
            maxInstances: 20
        secondaryWorkerConfig:
            maxInstances: 50
        basicAlgorithm:
            cooldownPeriod: 4m
            yarnConfig:
                scaleUpFactor: 0.05
                scaleDownFactor: 1.0
                gracefulDecommissionTimeout: 1h
        FIN
        gcloud dataproc autoscaling-policies import rc-example-autoscaling-policy --project=~{gcs_project} --source=auto-scale-policy.yaml --region=~{region}

        # Run the hail python script to make a VDS
        python3 /app/run_in_hail_cluster.py \
            --script-path /app/hail_gvs_import.py \
            --secondary-script-path /app/import_gvs.py \
            --account ${account_name} \
            --autoscaling-policy rc-example-autoscaling-policy \
            --region ~{region} \
            --gcs-project ~{gcs_project} \
            --cluster-name ${cluster_name} \
            --avro-path ~{avro_path} \
            --vds-path ~{vds_url} \
            --temp-path ~{temp_path}

        echo "Goodbye cluster"
    >>>

    output {
        String cluster_name = read_string("cluster_name.txt")
    }

    runtime {
        memory: "6.5 GB"
        disks: "local-disk 100 SSD"
        cpu: 1
        preemptible: false
        docker: variants_docker
        bootDiskSizeGb: 10
    }
}
