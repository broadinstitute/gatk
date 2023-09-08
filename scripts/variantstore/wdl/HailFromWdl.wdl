version 1.0
# Largely "borrowing" from Lee's work
# https://github.com/broadinstitute/aou-ancestry/blob/a57bbab3ccee4d06317fecb8ca109424bca373b7/script/wdl/hail_in_wdl/filter_VDS_and_shard_by_contig.wdl

#
# Given a VDS and a bed file, render a VCF (sharded by chromosome).
#   All bed files referenced in this WDL are UCSC bed files (as opposed to PLINK bed files).
#
# This has not been tested on any reference other than hg38.
# Inputs:
#
#         ## ANALYSIS PARAMETERS
#        #  ie, parameters that go to the Hail python code (submission_script below)
#        String vds_url
#
#        # Genomic region for the output VCFs to cover
#        File bed_url
#
#        # VCF Header that will be used in the output
#        File vcf_header_url
#
#        # Contigs of interest.  If a contig is present in the bed file, but not in this list, the contig will be ignored.
#        #   In other words, this is a contig level intersection with the bed file.
#        #     This list of contigs that must be present in the reference.  Each contig will be processed separately (shard)
#        # This list should be ordered.  Eg, ["chr21", "chr22"]
#        Array[String] contigs
#
#        # String used in construction of output filename
#        #  Cannot contain any special characters, ie, characters must be alphanumeric or "_"
#        String prefix
#
#        ## CLUSTER PARAMETERS
#        # Number of workers (per shard) to use in the Hail cluster.
#        Int num_workers
#
#        # Set to 'subnetwork' if running in Terra Cromwell
#        String gcs_subnetwork_name='subnetwork'
#
#        # The script that is run on the cluster
#        #  See filter_VDS_and_shard_by_contig.py for an example.
#        File submission_script
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

import "GvsUtils.wdl" as Utils

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filter_vds_to_VCF_by_chr {
    ### Change here:  You will need to specify all parameters (both analysis and runtime) that need to go to the
    # cluster, VM spinning up the cluster, and the script being run on the cluster.
    input {

        ## ANALYSIS PARAMETERS
        #  ie, parameters that go to the Hail python code (submission_script below)
        String vds_url

        # Genomic region for the output VCFs to cover
        File bed_url

        # VCF Header that will be used in the output
        File vcf_header_url

        # Contigs of interest.  If a contig is present in the bed file, but not in this list, the contig will be ignored.
        #   In other words, this is a contig level intersection with the bed file.
        #     This list of contigs that must be present in the reference.  Each contig will be processed separately (shard)
        # This list should be ordered.  Eg, ["chr21", "chr22"]
        Array[String] contigs

        # String used in construction of output filename
        #  Cannot contain any special characters, ie, characters must be alphanumeric or "_"
        String prefix

        ## CLUSTER PARAMETERS
        # Number of workers (per shard) to use in the Hail cluster.
        Int num_workers

        # Set to 'subnetwork' if running in Terra Cromwell
        String gcs_subnetwork_name='subnetwork'

        # The script that is run on the cluster
        #  See filter_VDS_and_shard_by_contig.py for an example.
        File submission_script

        # Set to "us-central1" if running in Terra Cromwell
        String region = "us-central1"
    }

    call Utils.GetToolVersions

    scatter (contig in contigs) {
        call filter_vds_and_export_as_vcf {
            input:
                vds_url = vds_url,
                bed_url = bed_url,
                contig = contig,
                prefix = prefix,
                gcs_project = GetToolVersions.google_project,
                num_workers = num_workers,
                gcs_subnetwork_name = gcs_subnetwork_name,
                vcf_header_url = vcf_header_url,
                submission_script = submission_script,
                hail_docker = GetToolVersions.hail_docker,
                region = region,
        }
    }

    output {
        Array[File] vcfs = filter_vds_and_export_as_vcf.vcf
    }
}

task filter_vds_and_export_as_vcf {
    input {
        # You must treat a VDS as a String, since it is a directory and not a single file
        String vds_url
        File bed_url
        File vcf_header_url

        # Cannot make this localization_optional
        File submission_script

        # contig must be in the reference
        String contig
        String prefix
        String gcs_project
        String region = "us-central1"
        Int num_workers
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String hail_docker
    }

    parameter_meta {
        bed_url: {
                     localization_optional: true
                 }
        vcf_header_url: {
                            localization_optional: true
                        }
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
            python3 -c "print('hello world')"
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
        from google.cloud import dataproc_v1 as dataproc

        # Must match pattern (?:[a-z](?:[-a-z0-9]{0,49}[a-z0-9])?)
        cluster_name = f'~{prefix}-~{contig}-hail-{str(uuid.uuid4())[0:13]}'

        # Must be local filepath
        script_path = "~{submission_script}"

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

            for cluster in cluster_client.list_clusters(request={"project_id": "~{gcs_project}", "region": "~{region}"}):
                if cluster.cluster_name == cluster_name:
                    cluster_staging_bucket = cluster.config.temp_bucket

                    #### THIS IS WHERE YOU CALL YOUR SCRIPT AND COPY THE OUTPUT LOCALLY (so that it can get back into WDL-space)
                    submit_cmd = f'gcloud dataproc jobs submit pyspark {script_path} --cluster={cluster_name} --project ~{gcs_project}  --region=~{region} --account {account} --driver-log-levels root=WARN -- --vds_url ~{vds_url} --bed_url ~{bed_url} --vcf_header_url ~{vcf_header_url} --contig ~{contig} --output_gs_url gs://{cluster_staging_bucket}/{cluster_name}/~{prefix}.~{contig}.vcf.bgz'
                    print("Running: " + submit_cmd)
                    os.popen(submit_cmd).read()
                    print("Copying results out of staging bucket...")
                    staging_cmd = f'gsutil cp -r gs://{cluster_staging_bucket}/{cluster_name}/~{prefix}.~{contig}.vcf.bgz ~{prefix}.~{contig}.vcf.bgz'
                    print(staging_cmd)
                    os.popen(staging_cmd).read()
                    ###########

                    break

        except Exception as e:
            print(e)
            raise
        finally:
            print(f'Stopping cluster: {cluster_name}')
            os.popen("gcloud dataproc clusters delete --project {} --region {} --account {} {}".format("~{gcs_project}", "~{region}", account, cluster_name)).read()

        EOF

        echo "Complete"
    >>>

    output {
        File vcf = "~{prefix}.~{contig}.vcf.bgz"
    }

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
