version 1.0

# This WDL counts loci present in a "new" VDS but absent from an "old" VDS, running
# the comparison in Hail on an autoscaling Dataproc cluster.
import "GvsUtils.wdl" as Utils

workflow GvsCountVdsNovelLoci {
    input {
        # Paths to the two VDSes to compare.
        String new_vds_path
        String old_vds_path

        String cluster_prefix = "vds-compare-cluster"
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
        new_vds_path: {
            help: "Path to the newer VDS whose novel loci will be counted (e.g. the v9 VDS)."
        }
        old_vds_path: {
            help: "Path to the older VDS used as the reference set of loci (e.g. the v8 VDS)."
        }
        cluster_prefix: {
            help: "Prefix of the Dataproc cluster name."
        }
        hail_temp_path: {
            help: "Optional GCS path for Hail temporary files."
        }
        region: {
            help: "GCP region in which to create the Dataproc cluster, e.g. us-central1."
        }
        hail_version: {
            help: "Optional Hail version. Cannot define both this parameter and `hail_wheel`."
        }
        hail_wheel: {
            help: "Optional Hail wheel file. Cannot define both this parameter and `hail_version`."
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

    call CountVdsNovelLoci {
        input:
            prefix = cluster_prefix,
            new_vds_path = new_vds_path,
            old_vds_path = old_vds_path,
            hail_version = effective_hail_version,
            hail_wheel = hail_wheel,
            hail_temp_path = hail_temp_path,
            run_in_hail_cluster_script = GetHailScripts.run_in_hail_cluster_script,
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
        String cluster_name = CountVdsNovelLoci.cluster_name
        Int novel_loci_count = CountVdsNovelLoci.novel_loci_count
        Boolean done = true
    }
}

task CountVdsNovelLoci {
    input {
        String prefix
        String new_vds_path
        String old_vds_path
        Boolean leave_cluster_running_at_end
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

    meta {
        # Always run: results may differ as VDSes are updated.
        volatile: true
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

        # Write the Hail analysis script inline.
        cat > compare_vds_loci.py << 'PYTHON_SCRIPT'
        import argparse
        import hail as hl


        if __name__ == '__main__':
            parser = argparse.ArgumentParser(allow_abbrev=False,
                                             description='Count loci present in a "new" VDS but absent from an "old" VDS.')
            parser.add_argument('--new-vds-path', type=str, required=True,
                                help='Path to the newer VDS whose novel loci will be counted.')
            parser.add_argument('--old-vds-path', type=str, required=True,
                                help='Path to the older VDS used as the reference set of loci.')
            parser.add_argument('--temp-path', type=str, required=True,
                                help='Path to a GCS directory for Hail temporary files.')
            parser.add_argument('--output-count-path', type=str, required=True,
                                help='Path to a GCS file to write the novel locus count to.')

            args = parser.parse_args()

            hl.init(tmp_dir=f'{args.temp_path}/hail_tmp_general')
            hl.default_reference('GRCh38')

            new_vds = hl.vds.read_vds(args.new_vds_path)
            old_vds = hl.vds.read_vds(args.old_vds_path)

            new_loci = new_vds.variant_data.rows().select().key_by('locus')
            old_loci = old_vds.variant_data.rows().select().key_by('locus')

            # Note that all GVS-produced VDSes will be unsplit, with all alleles for a locus on a single row. Thus we
            # can determine which loci are new with this anti-join.
            new_only = new_loci.anti_join(old_loci)
            count = new_only.count()

            print(f'Loci in new VDS ({args.new_vds_path}) but not in old VDS ({args.old_vds_path}): {count}')

            with hl.hadoop_open(args.output_count_path, 'w') as out:
                out.write(str(count) + '\n')
        PYTHON_SCRIPT

        # Construct a JSON of arguments for the Python script to be run in the Hail cluster.
        cat > script-arguments.json <<FIN
        {
            "new-vds-path": "~{new_vds_path}",
            "old-vds-path": "~{old_vds_path}",
            "temp-path": "${hail_temp_path}",
            "output-count-path": "${hail_temp_path}/novel_count.txt"
        }
        FIN

        # Run the Hail Python script to compare loci between the two VDSes.
        # The locus count is written to the path specified by `output-count-path`.
        python3 ~{run_in_hail_cluster_script} \
            --script-path compare_vds_loci.py \
            --script-arguments-json-path script-arguments.json \
            --account ${account_name} \
            --region ~{region} \
            --workspace-project ~{workspace_project} \
            --cluster-name ${cluster_name} \
            ~{'--cluster-max-idle-minutes ' + cluster_max_idle_minutes} \
            ~{'--cluster-max-age-minutes ' + cluster_max_age_minutes} \
            ~{'--master-memory-fraction ' + master_memory_fraction} \
            ~{true='--leave-cluster-running-at-end' false='' leave_cluster_running_at_end}

        gsutil cp "${hail_temp_path}/novel_count.txt" .
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
        Int novel_loci_count = read_int("novel_count.txt")
        Boolean done = true
    }
}

