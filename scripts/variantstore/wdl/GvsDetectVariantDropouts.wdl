version 1.0

# This WDL detects variant-call dropout regions across the entire genome for
# every superpartition in a GVS Hail VDS.  It runs genome_wide_dropouts.py on
# a Hail Dataproc cluster following the same pattern as GvsValidateVDS.wdl.
import "GvsUtils.wdl" as Utils

workflow GvsDetectVariantDropouts {
    input {
        # ---------------------------------------------------------------
        # Required
        # ---------------------------------------------------------------
        String vds_path
        # GCS path to the TSV produced by the BigQuery superpartition-sampling
        # query (columns: sample_name, sample_id, superpartition).
        String samples_path
        # GCS path where the output TSV will be written.
        String output_path

        # ---------------------------------------------------------------
        # Optional tuning
        # ---------------------------------------------------------------
        # Genomic bin size in base-pairs.
        Int    bin_size              = 50000
        # Bins below this fraction of the per-superpartition/contig median are
        # flagged as dropouts.
        Float  dropout_fraction      = 0.5

        # ---------------------------------------------------------------
        # Cluster / infra
        # ---------------------------------------------------------------
        String  cluster_prefix             = "dropouts-cluster"
        String? hail_temp_path
        String  region                     = "us-central1"
        Int?    cluster_max_idle_minutes
        Int?    cluster_max_age_minutes
        Boolean leave_cluster_running_at_end = false
        Float?  master_memory_fraction

        # ---------------------------------------------------------------
        # Docker / version overrides (all optional; resolved via
        # GetToolVersions when not supplied)
        # ---------------------------------------------------------------
        String? git_branch_or_tag
        String? hail_version
        File?   hail_wheel
        String? basic_docker
        String? variants_docker
        String? workspace_bucket
        String? workspace_project
        String? cloud_sdk_slim_docker
    }

    parameter_meta {
        vds_path: {
            help: "GCS path to the Hail VDS to inspect for dropouts."
        }
        samples_path: {
            help: "GCS path to the superpartition samples TSV produced by the BigQuery sampling query (columns: sample_name, sample_id, superpartition)."
        }
        output_path: {
            help: "GCS path where the output dropout TSV will be written (columns: contig, bin_start, bin_end, superpartition, n_variants, median_bin_count, dropout_flag)."
        }
        bin_size: {
            help: "Genomic bin size in base-pairs (default: 50000)."
        }
        dropout_fraction: {
            help: "Bins whose variant count is below this fraction of the per-superpartition/contig median are flagged as dropouts (default: 0.5)."
        }
        cluster_prefix: {
            help: "Prefix for the Dataproc cluster name."
        }
        hail_temp_path: {
            help: "Optional GCS path to use as Hail's temporary directory."
        }
        region: {
            help: "GCP region in which to create the Dataproc cluster (default: us-central1)."
        }
        hail_version: {
            help: "Optional Hail version to install. Cannot be specified together with hail_wheel."
        }
        hail_wheel: {
            help: "Optional Hail wheel file. Cannot be specified together with hail_version."
        }
    }

    if (!defined(variants_docker) || !defined(basic_docker) || !defined(cloud_sdk_slim_docker) || !defined(workspace_bucket) || !defined(workspace_project) || !defined(hail_version)) {
        call Utils.GetToolVersions {
            input:
                git_branch_or_tag = git_branch_or_tag,
        }
    }

    String effective_workspace_bucket    = select_first([workspace_bucket,    GetToolVersions.workspace_bucket])
    String effective_google_project      = select_first([workspace_project,   GetToolVersions.google_project])
    String effective_basic_docker        = select_first([basic_docker,         GetToolVersions.basic_docker])
    String effective_variants_docker     = select_first([variants_docker,      GetToolVersions.variants_docker])
    String effective_cloud_sdk_slim_docker = select_first([cloud_sdk_slim_docker, GetToolVersions.cloud_sdk_slim_docker])

    if (defined(hail_version) && defined(hail_wheel)) {
        call Utils.TerminateWorkflow as BothHailVersionAndHailWheelDefined {
            input:
                message      = "Cannot define both `hail_version` and `hail_wheel`, exiting.",
                basic_docker = effective_basic_docker,
        }
    }

    String effective_hail_version = select_first([hail_version, GetToolVersions.hail_version])

    call Utils.GetHailScripts {
        input:
            variants_docker = effective_variants_docker,
    }

    call DetectDropouts {
        input:
            prefix                        = cluster_prefix,
            vds_path                      = vds_path,
            samples_path                  = samples_path,
            output_path                   = output_path,
            bin_size                      = bin_size,
            dropout_fraction              = dropout_fraction,
            hail_version                  = effective_hail_version,
            hail_wheel                    = hail_wheel,
            hail_temp_path                = hail_temp_path,
            run_in_hail_cluster_script    = GetHailScripts.run_in_hail_cluster_script,
            genome_wide_dropouts_script   = GetHailScripts.genome_wide_dropouts_script,
            workspace_project             = effective_google_project,
            region                        = region,
            workspace_bucket              = effective_workspace_bucket,
            leave_cluster_running_at_end  = leave_cluster_running_at_end,
            cluster_max_idle_minutes      = cluster_max_idle_minutes,
            cluster_max_age_minutes       = cluster_max_age_minutes,
            master_memory_fraction        = master_memory_fraction,
            cloud_sdk_slim_docker         = effective_cloud_sdk_slim_docker,
    }

    output {
        String cluster_name   = DetectDropouts.cluster_name
        String dropout_output = output_path
        Boolean done          = true
    }
}

task DetectDropouts {
    input {
        String  prefix
        String  vds_path
        String  samples_path
        String  output_path
        Int     bin_size
        Float   dropout_fraction

        String? hail_version
        File?   hail_wheel
        String? hail_temp_path

        File    run_in_hail_cluster_script
        File    genome_wide_dropouts_script

        String  workspace_project
        String  workspace_bucket
        String  region

        Boolean leave_cluster_running_at_end
        Int?    cluster_max_idle_minutes
        Int?    cluster_max_age_minutes
        Float?  master_memory_fraction

        String  cloud_sdk_slim_docker
    }

    meta {
        # Always re-run: output depends on live GCS data.
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

        # Generate a UUIDish random hex string of <8 hex chars>-<4 hex chars>
        hex="$(head -c4 < /dev/urandom | od -h -An | tr -d '[:space:]')-$(head -c2 < /dev/urandom | od -h -An | tr -d '[:space:]')"

        cluster_name="~{prefix}-${hex}"
        echo ${cluster_name} > cluster_name.txt

        if [[ -z "~{hail_temp_path}" ]]
        then
            hail_temp_path="~{workspace_bucket}/hail-temp/hail-temp-${hex}"
        else
            hail_temp_path="~{hail_temp_path}"
        fi

        # Set up the autoscaling policy.
        # No secondary (preemptible) workers: per-chromosome jobs complete in
        # ~2 minutes each, so we don't need the extra capacity, and preemptible
        # nodes being hard-killed mid-shuffle cause FetchFailedException.
        cat > auto-scale-policy.yaml <<FIN
        workerConfig:
            minInstances: 2
            maxInstances: 2
        secondaryWorkerConfig:
            maxInstances: 0
        basicAlgorithm:
            cooldownPeriod: 120s
            yarnConfig:
                scaleUpFactor: 1.0
                scaleDownFactor: 1.0
                gracefulDecommissionTimeout: 120s
        FIN
        gcloud dataproc autoscaling-policies import gvs-autoscaling-policy \
            --project=~{workspace_project} \
            --source=auto-scale-policy.yaml \
            --region=~{region} \
            --quiet

        # Build the JSON argument file for genome_wide_dropouts.py
        cat > script-arguments.json <<FIN
        {
            "vds-path":          "~{vds_path}",
            "temp-path":         "${hail_temp_path}",
            "samples-path":      "~{samples_path}",
            "output-path":       "~{output_path}",
            "bin-size":          ~{bin_size},
            "dropout-fraction":  ~{dropout_fraction}
        }
        FIN

        python3 ~{run_in_hail_cluster_script} \
            --script-path ~{genome_wide_dropouts_script} \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            --autoscaling-policy gvs-autoscaling-policy \
            --region ~{region} \
            --gcs-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes '  + cluster_max_age_minutes} \
            ~{'--master-memory-fraction '   + master_memory_fraction} \
            --job-max-failures-total 1 \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end}
    >>>

    runtime {
        memory:          "6.5 GB"
        disks:           "local-disk 100 SSD"
        cpu:             1
        preemptible:     0
        docker:          cloud_sdk_slim_docker
        bootDiskSizeGb:  10
    }

    output {
        String cluster_name = read_string("cluster_name.txt")
        Boolean done = true
    }
}

