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
#        String bed_url
#
#        # VCF Header that will be used in the output
#        String vcf_header_url
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

        String? git_branch_or_tag
        String? hail_version
        String? worker_machine_type

        # Genomic region for the output VCFs to cover
        String bed_url = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"

        # VCF Header that will be used in the output
        String vcf_header_url = "gs://gvs_quickstart_storage/hail_from_wdl/vcf_header.txt"

        # Contigs of interest.  If a contig is present in the bed file, but not in this list, the contig will be ignored.
        #   In other words, this is a contig level intersection with the bed file.
        #     This list of contigs that must be present in the reference.  Each contig will be processed separately (shard)
        # This list should be ordered.  Eg, ["chr21", "chr22"]
        Array[String] contigs = ["chr20"]

        # String used in construction of output filename
        #  Cannot contain any special characters, ie, characters must be alphanumeric or "-"
        String prefix = "hail-from-wdl"

        ## CLUSTER PARAMETERS
        # Number of workers (per shard) to use in the Hail cluster.
        Int num_workers = 10

        # Set to 'subnetwork' if running in Terra Cromwell
        String gcs_subnetwork_name = 'subnetwork'

        # The script that is run on the cluster
        #  See filter_VDS_and_shard_by_contig.py for an example.
        File? submission_script

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
                git_branch_or_tag = git_branch_or_tag,
                hail_version = hail_version,
                worker_machine_type = worker_machine_type,
                submission_script = submission_script,
                variants_docker = GetToolVersions.variants_docker,
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
        String bed_url
        String vcf_header_url

        String? git_branch_or_tag
        File? submission_script
        String? hail_version
        String? worker_machine_type

        # contig must be in the reference
        String contig
        String prefix
        String gcs_project
        String region = "us-central1"
        Int num_workers
        RuntimeAttr? runtime_attr_override
        String gcs_subnetwork_name

        String variants_docker
    }

    RuntimeAttr runtime_default = object {
                                      mem_gb: 30,
                                      disk_gb: 100,
                                      cpu_cores: 1,
                                      preemptible_tries: 0,
                                      max_retries: 0,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    String default_script_filename = "filter_VDS_and_shard_by_contig.py"

    command <<<
        # Prepend date, time and pwd to xtrace log entries.
        PS4='\D{+%F %T} \w $ '
        set -o errexit -o nounset -o pipefail -o xtrace

        account_name=$(gcloud config list account --format "value(core.account)")

        pip3 install --upgrade pip
        pip3 install hail~{'==' + hail_version}
        pip3 install --upgrade google-cloud-dataproc

        if [[ -z "~{git_branch_or_tag}" && -z "~{submission_script}" ]] || [[ ! -z "~{git_branch_or_tag}" && ! -z "~{submission_script}" ]]
        then
            echo "Must specify git_branch_or_tag XOR submission_script"
            exit 1
        elif [[ ! -z "~{git_branch_or_tag}" ]]
        then
            script_url="https://raw.githubusercontent.com/broadinstitute/gatk/~{git_branch_or_tag}/scripts/variantstore/wdl/extract/~{default_script_filename}"
            curl --silent --location --remote-name "${script_url}"
        fi

        if [[ ! -z "~{submission_script}" ]]
        then
            script_path="~{submission_script}"
        else
            script_path="~{default_script_filename}"
        fi

        python3 /app/run_in_hail_cluster.py \
            --script-path ${script_path} \
            --account ${account_name} \
            --num-workers ~{num_workers} \
            ~{'--worker-machine-type' + worker_machine_type} \
            --region ~{region} \
            --gcs-project ~{gcs_project} \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --vds-url ~{vds_url} \
            --vcf-header-url ~{vcf_header_url} \
            --bed-url ~{bed_url}

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
        # `slim` here to be able to use Java
        docker: variants_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}
