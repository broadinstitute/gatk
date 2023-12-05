import argparse
import os
import uuid
from google.cloud import dataproc_v1 as dataproc
from logging import info


def configure_logging():
    import logging
    import sys
    # https://stackoverflow.com/a/14058475
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root.addHandler(handler)


def unwrap(string):
    import re
    return re.sub("\\s{2,}", " ", string).strip()


def run_in_cluster(cluster_name, account, worker_machine_type, master_machine_type, region, gcs_project,
                   autoscaling_policy, script_path, secondary_script_path, use_classic_vqsr, vds_path, temp_path,
                   avro_path, leave_cluster_running_at_end, cluster_max_idle_minutes, cluster_max_age_minutes,
                   master_memory_fraction):

    use_classic_vqsr_flag = ""
    if use_classic_vqsr:
        use_classic_vqsr_flag = "--use-classic-vqsr"

    cluster_max_idle_arg = f"--max-idle {cluster_max_idle_minutes}m" if cluster_max_idle_minutes else ""
    cluster_max_age_arg = f"--max-age {cluster_max_age_minutes}m" if cluster_max_age_minutes else ""

    try:
        cluster_start_cmd = unwrap(f"""
        
        hailctl dataproc start 
         --autoscaling-policy={autoscaling_policy}
         --worker-machine-type {worker_machine_type}
         --master-machine-type {master_machine_type}
         --master-memory-fraction {master_memory_fraction}
         --region {region}
         --project {gcs_project}
         --service-account {account}
         {cluster_max_idle_arg}
         {cluster_max_age_arg}
         --num-master-local-ssds 1
         --num-worker-local-ssds 1 
         --subnet=projects/{gcs_project}/regions/{region}/subnetworks/subnetwork
         --properties=dataproc:dataproc.monitoring.stackdriver.enable=true,dataproc:dataproc.logging.stackdriver.enable=true,core:fs.gs.outputstream.sync.min.interval=5
         {cluster_name}
         
        """)

        info(f"Starting cluster '{cluster_name}'...")
        info(cluster_start_cmd)
        pipe = os.popen(cluster_start_cmd)
        info(pipe.read())
        wait_status = pipe.close()
        if wait_status:
            exit_code = os.waitstatus_to_exitcode(wait_status)
            raise RuntimeError(f"Unexpected exit code from cluster creation: {exit_code}")

        cluster_client = dataproc.ClusterControllerClient(
            client_options={"api_endpoint": f"{region}-dataproc.googleapis.com:443"}
        )

        for cluster in cluster_client.list_clusters(request={"project_id": gcs_project, "region": region}):
            if cluster.cluster_name == cluster_name:
                submit_cmd = unwrap(f"""

                gcloud dataproc jobs submit pyspark {script_path}
                 {'--py-files=' + secondary_script_path if secondary_script_path else ''}
                 --cluster={cluster_name}
                 --project {gcs_project}
                 --region={region}
                 --account {account}
                 --driver-log-levels root=WARN
                 --
                 --vds-path {vds_path}
                 --temp-path {temp_path}
                 {'--avro-path ' + avro_path if avro_path else ''}
                 {use_classic_vqsr_flag}
                """)

                info("Running: " + submit_cmd)
                pipe = os.popen(submit_cmd)
                pipe.read()
                wait_status = pipe.close()
                if wait_status:
                    exit_code = os.waitstatus_to_exitcode(wait_status)
                    raise RuntimeError(f"Unexpected exit code running submitted job: {exit_code}")
                break
    except Exception as e:
        info(e)
        raise
    finally:
        if leave_cluster_running_at_end:
            info(f"Leaving cluster {cluster_name} running as `leave_cluster_running_at_end` option is True.")
        else:
            info(f'Stopping cluster: {cluster_name}')
            delete_cmd = unwrap(f"""

                gcloud dataproc clusters delete
                  --project {gcs_project}
                  --region {region}
                  --account {account}
                  --quiet
                  {cluster_name}

            """)

            pipe = os.popen(delete_cmd)
            pipe.read()
            wait_status = pipe.close()
            if wait_status:
                exit_code = os.waitstatus_to_exitcode(wait_status)
                raise RuntimeError(f"Unexpected exit code deleting cluster: {exit_code}")


if __name__ == "__main__":
    configure_logging()

    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--cluster-name', type=str, required=True, help='Name of the Hail cluster')
    parser.add_argument('--account', type=str, help='GCP account name')
    parser.add_argument('--worker-machine-type', type=str, required=False, default="n1-highmem-8",
                        help='Dataproc cluster worker machine type')
    parser.add_argument('--master-machine-type', type=str, required=False, default="n1-highmem-16",
                        help='Dataproc cluster master machine type')
    parser.add_argument('--master-memory-fraction', type=float, default=0.8, help='Dataproc master memory fraction')
    parser.add_argument('--region', type=str, required=True, help='GCS region')
    parser.add_argument('--gcs-project', type=str, required=True, help='GCS project')
    parser.add_argument('--autoscaling-policy', type=str, help='Name of the autoscaling policy that should get used')
    parser.add_argument('--script-path', type=str, required=True, help='Path to script to run in Hail cluster')
    parser.add_argument('--secondary-script-path', type=str, help='Path to secondary script to run in Hail cluster')
    parser.add_argument("--use-classic-vqsr", action="store_true",
                        help="If set, expect that the input GVS Avro files were generated using VQSR Classic")
    parser.add_argument('--vds-path', type=str, required=True, help='VDS URL')
    parser.add_argument('--avro-path', type=str, help='Avro URL')
    parser.add_argument('--temp-path', type=str, required=True, help='Cruft URL')
    parser.add_argument('--cluster-max-idle-minutes', type=int, help='Maximum idle time of cluster in minutes')
    parser.add_argument('--cluster-max-age-minutes', type=int, help='Maximum age of cluster in minutes')
    parser.add_argument('--leave-cluster-running-at-end', action="store_true", default=False)

    args = parser.parse_args()

    run_in_cluster(cluster_name=args.cluster_name,
                   account=args.account,
                   master_machine_type=args.master_machine_type,
                   worker_machine_type=args.worker_machine_type,
                   region=args.region,
                   gcs_project=args.gcs_project,
                   autoscaling_policy=args.autoscaling_policy,
                   script_path=args.script_path,
                   secondary_script_path=args.secondary_script_path,
                   use_classic_vqsr=args.use_classic_vqsr,
                   vds_path=args.vds_path,
                   temp_path=args.temp_path,
                   avro_path=args.avro_path,
                   leave_cluster_running_at_end=args.leave_cluster_running_at_end,
                   cluster_max_idle_minutes=args.cluster_max_idle_minutes,
                   cluster_max_age_minutes=args.cluster_max_age_minutes,
                   master_memory_fraction=args.master_memory_fraction
                   )
