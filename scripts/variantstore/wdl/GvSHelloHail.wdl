version 1.0
#
# For the actual script that will run the VDS creation in Hail--the input will need to be a path to avro files and many parameters
# For now I just want a cluster and a print out in the script
# This has not been tested on any reference other than hg38.
#
# Inputs:
#
#         ## ANALYSIS PARAMETERS
#        #  ie, parameters that go to the Hail python code (submission_script below)
#        String vds_url
#
#        # String used in construction of output filename
#        #  Cannot contain any special characters, ie, characters must be alphanumeric or "_"
#        String prefix
#
#        ## CLUSTER PARAMETERS
#        # Number of workers (per shard) to use in the Hail cluster.
#        Int num_workers
#
#        # The Google project ID information is necessary when spinning up dataproc.
#        #  ie, terra-<hex>
#        # This must match the workspace that this workflow is being run from.
#        #     eg, "terra-491d5f31"
#        # once Miguel's code has been merged, this will be available as a lookup task in GvsUtils
#        String gcs_project
#
#        # Set to 'subnetwork' if running in Terra Cromwell
#        String gcs_subnetwork_name='subnetwork'
#
#        # The script that is run on the cluster
#        #  See filter_VDS_and_shard_by_contig.py for an example.
#        # File submission_script
#
#        # Set to "us-central1" if running in Terra Cromwell
#        String region = "us-central1"
#
#        ## VM PARAMETERS
#        # Please note that there is a RuntimeAttr struct and a task parameter that can be used to override the defaults
#        #  of the VM.  These are task parameters.
#        #  However, since this can be a lightweight VM, overriding is unlikely to be necessary.
#
#        # The docker to be used on the VM.  This will need both Hail and Google Cloud SDK installed.
#        String hail_docker="us.gcr.io/broad-dsde-methods/lichtens/hail_dataproc_wdl:1.0"
#
# Important notes:
#   - Hail will save the VCFs in the cloud.  You will need to provide this storage space.  In other words, the runtime
#     parameters must have enough storage space to support a single contig
#   - This WDL script is still dependent on the python/Hail script that it calls.  You will see this when the parameters
#    are passed into the script.
#   - This WDL is boilerplate, except for input parameters, output parameters, and where marked in the main task.
#   - We HIGHLY recommend that the WDL do NOT run on a preemptible VM
#    (reminder, this is a single VM that spins up the dataproc cluster and submits jobs -- it is not doing any of the
#     actual computation.  In other words, it does not need to be a heavy machine.)
#     In other words, always set `preemptible_tries` to zero (default).
#
struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow GvSHelloHail {
    ### Change here:  You will need to specify all parameters (both analysis and runtime) that need to go to the
    # cluster, VM spinning up the cluster, and the script being run on the cluster.
    input {

        ## ANALYSIS PARAMETERS
        #  ie, parameters that go to the Hail python code (submission_script below)
        # String vds_url

        # String used in construction of output filename
        #  Cannot contain any special characters, ie, characters must be alphanumeric or "-"
        String prefix = "rc-test-cluster"

        ## CLUSTER PARAMETERS
        # Number of workers (per shard) to use in the Hail cluster.
        Int num_workers

        # The Google project ID information is necessary when spinning up dataproc.
        String gcs_project

        # Set to 'subnetwork' if running in Terra Cromwell
        String gcs_subnetwork_name='subnetwork'

        # The script that is run on the cluster
        #  See filter_VDS_and_shard_by_contig.py for an example.
        # File submission_script0

        # Set to "us-central1" if running in Terra Cromwell
        String region = "us-central1"

        ## VM PARAMETERS
        # Please note that there is a RuntimeAttr struct and a task parameter that can be used to override the defaults
        #  of the VM.  These are task parameters.
        #  However, since this can be a lightweight VM, overriding is unlikely to be necessary.

        # The docker to be used on the VM.  This will need both Hail and Google Cloud SDK installed.
        String hail_docker="us.gcr.io/broad-dsde-methods/lichtens/hail_dataproc_wdl:1.1"
    }

    call say_hello_hail {
        input:
            # vds_url = vds_url,
            prefix = prefix,
            gcs_project = gcs_project,
            num_workers = num_workers,
            gcs_subnetwork_name = gcs_subnetwork_name,
            # submission_script1 = submission_script2,
            hail_docker = hail_docker,
    }

}

task say_hello_hail {
    input {
        # Goodnight Gracie was already taken

        # We will eventually want a script here
        # File submission_script3
        String prefix

        # Cluster params
        String gcs_project
        String region = "us-central1"
        Int num_workers
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
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

    command <<<
        set -euxo pipefail

        gcloud config list account --format "value(core.account)" 1> account.txt

        #### TEST:  Make sure that this docker image is configured for python3
        if which python3 > /dev/null 2>&1; then
            pt3="$(which python3)"
            echo "** python3 located at $pt3"
            echo "** magic: $(file $pt3)"
            echo "** Version info:"
            echo "$(python3 -V)"
            echo "** -c test"
            python3 -c "print('hello world-- no cluster yet tho')"
        else
            echo "!! No 'python3' in path."
            exit 1
        fi
        #### END TEST

        python3 <<EOF
        print("Running python code...")
        import hail as hl
        import os
        import uuid
        import time
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-hail-{str(uuid.uuid4())[0:13]}'

        # Must be local filepath once a script is finally set

        with open("account.txt", "r") as account_file:
            account = account_file.readline().strip()
        print("account: " + account)

        try:
            cluster_start_cmd = "hailctl dataproc start --num-workers ~{num_workers} --region {} --project {} --service-account {} --num-master-local-ssds 1 --num-worker-local-ssds 1 --max-idle=60m --max-age=1440m --subnet={} {}".format("~{region}", "~{gcs_project}", account, "projects/~{gcs_project}/regions/~{region}/subnetworks/~{gcs_subnetwork_name}", cluster_name)
            print("Starting cluster...")
            print(cluster_start_cmd)
            print(os.popen(cluster_start_cmd).read())

            cluster_client = dataproc.ClusterControllerClient(
                client_options={"api_endpoint": f"~{region}-dataproc.googleapis.com:443"}
            )

            print("Hello cluster!")

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_staging_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    ## Maybe next step is to make a simple python script to run aside from just printing hello world
                    # submit_cmd = f'go find me in a commit i hate dockstore'
                    # print("Running: " + submit_cmd)
                    print("Running nothing yet")
                    # os.popen(submit_cmd).read()
                    # print("Copying results out of staging bucket...")
                    # staging_cmd = f'gsutil cp -r gs://{cluster_staging_bucket}/{cluster_name}/~{prefix}.vcf.bgz ~{prefix}.vcf.bgz'
                    # print(staging_cmd)
                    # os.popen(staging_cmd).read()
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            time.sleep(300)
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Goodbye cluster"
    >>>

    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " SSD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}