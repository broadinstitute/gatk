version 1.0
#
# This script will run the VDS creation in Hail
# the input will need to be a path to avro files and many parameters
#
# This has not been tested on any reference other than hg38.
#
# Inputs:
#
#         ## ANALYSIS PARAMETERS
#        #  ie, parameters that go to the Hail python code (submission_script below)
#        # location for the avro files
#        String avro_path
#        # location for the final, created VDS
#        String vds_output_url
#
#        ## CLUSTER PARAMETERS
#        # Number of workers (per shard) to use in the Hail cluster.
#        Int num_workers
#
#        # Set to 'subnetwork' if running in Terra Cromwell
#        String gcs_subnetwork_name='subnetwork'
#
#        # Set to "us-central1" if running in Terra Cromwell
#        String region = "us-central1"
#
#        ## VM PARAMETERS
#        # Please note that there is a RuntimeAttr struct and a task parameter that can be used to override the defaults
#        #  of the VM.  These are task parameters.
#        #  However, since this can be a lightweight VM, overriding is unlikely to be necessary.
#
#
# Important notes:
#   - This WDL script is still dependent on the python/Hail script that it calls.  You will see this when the parameters
#    are passed into the script.
#   - This WDL is boilerplate, except for input parameters, output parameters, and where marked in the main task.
#   - We HIGHLY recommend that the WDL do NOT run on a preemptible VM
#    (reminder, this is a single VM that spins up the dataproc cluster and submits jobs -- it is not doing any of the
#     actual computation.  In other words, it does not need to be a heavy machine.)
#     In other words, always set `preemptible_tries` to zero (default).
#

import "GvsUtils.wdl" as Utils

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow GvsCreateVDS {
    ### Change here:  You will need to specify all parameters (both analysis and runtime) that need to go to the
    # cluster, VM spinning up the cluster, and the script being run on the cluster.
    input {

        ## ANALYSIS PARAMETERS
        #  ie, parameters that go to the Hail python code (submission_script below)
        String vds_output_url
        String avro_path
        String? hail_version='0.2.122'

        # String used in construction of output filename
        #  Cannot contain any special characters, ie, characters must be alphanumeric or "-"
        String prefix = "rc-test-cluster" ## TODO: what are we going to want to hardcode this as?

        ## CLUSTER PARAMETERS
        # Number of workers (per shard) to use in the Hail cluster.
        Int num_workers

        # Set to 'subnetwork' if running in Terra Cromwell
        String gcs_subnetwork_name='subnetwork'

        # Set to "us-central1" if running in Terra Cromwell
        String region = "us-central1"

        ## VM PARAMETERS
        # Please note that there is a RuntimeAttr struct and a task parameter that can be used to override the defaults
        #  of the VM.  These are task parameters.
        #  However, since this can be a lightweight VM, overriding is unlikely to be necessary.
    }

    call Utils.GetToolVersions

    call create_vds {
        input:
            vds_url = vds_output_url,
            avro_path = avro_path,
            hail_version = hail_version,
            prefix = prefix,
            gcs_project = GetToolVersions.google_project,
            gcs_bucket = GetToolVersions.workspace_bucket,
            num_workers = num_workers,
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
        String gcs_project  # The Google project ID information is necessary when spinning up dataproc.
        String gcs_bucket
        String region = "us-central1"
        Int num_workers
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String variants_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 6.5,
                                      disk_gb: 100,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    String temp_path = "gs://fc-eada2674-7c2b-42a6-8db3-0246872596dc/quickstart-vds-for-wdl-tieout/temp-dir/" ## use the gcs_bucket here

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

        cluster_name="~{prefix}-hail-${hex}"
        echo ${cluster_name} > cluster_name.txt

        gcloud storage cp gs://fc-d5e319d4-b044-4376-afde-22ef0afc4088/auto-scale-policy.yaml auto-scale-policy.yaml
        gcloud dataproc autoscaling-policies import rc-example-autoscaling-policy --project=~{gcs_project} --source=auto-scale-policy.yaml --region=~{region}


        # get the scripts that we will use for the VDS creation
        curl --silent --location --remote-name https://raw.githubusercontent.com/broadinstitute/gatk/rc-vs-1025-hail-version-120/scripts/variantstore/wdl/extract/hail_gvs_import.py --output hail_gvs_import.py
        curl --silent --location --remote-name https://raw.githubusercontent.com/broadinstitute/gatk/rc-vs-1025-hail-version-120/scripts/variantstore/wdl/extract/import_gvs.py --output import_gvs.py
        curl --silent --location --remote-name https://raw.githubusercontent.com/broadinstitute/gatk/13b93fb14df32b3d1fcd88a290d9873f13bebc18/scripts/variantstore/wdl/extract/run_in_hail_cluster.py --output run_in_hail_cluster.py

        # Must be local filepath once a script is finally set

        python3 run_in_hail_cluster.py \
            --script-path hail_gvs_import.py \
            --account ${account_name} \
            --num-workers ~{num_workers} \
            --region ~{region} \
            --gcs-project ~{gcs_project} \
            --cluster-name ${cluster_name} \
            --prefix ~{prefix} \
            --avro-path ~{avro_path} \
            --vds-path ~{vds_url} \
            --temp-path ~{temp_path}

        echo "Goodbye cluster"
    >>>

    output {
        String cluster_name = read_string("cluster_name.txt")
    }

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: variants_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}