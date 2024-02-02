import argparse
import ijson
import os
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
                   autoscaling_policy, script_path, secondary_script_path_list, script_arguments_json_path, leave_cluster_running_at_end, cluster_max_idle_minutes, cluster_max_age_minutes, master_memory_fraction):

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

        # prepare custom arguments
        secondary_script_path_arg = f'--py-files {" ".join(secondary_script_path_list)}' if secondary_script_path_list else ''
        with open(script_arguments_json_path, 'r') as input_file:
            items = ijson.items(input_file, '', use_float=True)
            arguments = items.__next__();
            custom_script_args = [f"--{key} {arguments.get(key)}" for key in arguments.keys()]

        for cluster in cluster_client.list_clusters(request={"project_id": gcs_project, "region": region}):
            if cluster.cluster_name == cluster_name:
                submit_cmd = unwrap(f"""

                gcloud dataproc jobs submit pyspark {script_path}
                 {secondary_script_path_arg}
                 --cluster={cluster_name}
                 --project {gcs_project}
                 --region={region}
                 --account {account}
                 --driver-log-levels root=WARN
                 --
                 {"\n".join(custom_script_args)}
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
    parser.add_argument('--secondary-script-path-list', type=str, required=False, action="append", default=[],
                        help='List of paths to secondary scripts to run in Hail cluster')
    parser.add_argument('--script-arguments-json-path', type=str, required=True,
                        help='JSON file of arguments for script')
    parser.add_argument('--temp-path', type=str, required=True, help='Cruft URL')
    parser.add_argument('--cluster-max-idle-minutes', type=int, help='Maximum idle time of cluster in minutes')
    parser.add_argument('--cluster-max-age-minutes', type=int, help='Maximum age of cluster in minutes')
    parser.add_argument('--intermediate-resume-point', type=int, required=False,
                        help='Intermediate VDS index at which to resume')

    args = parser.parse_args()

    run_in_cluster(cluster_name=args.cluster_name,
                   account=args.account,
                   master_machine_type=args.master_machine_type,
                   worker_machine_type=args.worker_machine_type,
                   region=args.region,
                   gcs_project=args.gcs_project,
                   autoscaling_policy=args.autoscaling_policy,
                   script_path=args.script_path,
                   secondary_script_path_list=args.secondary_script_path_list,
                   script_arguments_json_path=args.script_arguments_json_path,
                   leave_cluster_running_at_end=args.leave_cluster_running_at_end,
                   cluster_max_idle_minutes=args.cluster_max_idle_minutes,
                   cluster_max_age_minutes=args.cluster_max_age_minutes,
                   master_memory_fraction=args.master_memory_fraction,
                   )
