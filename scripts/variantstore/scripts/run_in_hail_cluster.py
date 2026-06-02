import argparse
import ijson
import os
import tempfile
from google.cloud import dataproc_v1 as dataproc
from logging import info


TINY_NUM_WORKERS = 2
TINY_MAX_SECONDARY = 5
DEFAULT_NUM_WORKERS = 2
DEFAULT_MAX_SECONDARY = 200


def autoscaling_policy_name(num_workers, max_secondary):
    """Return a policy name that encodes its configuration, e.g. gvs-spark-autoscaling-w2-s200."""
    return f"gvs-spark-autoscaling-w{num_workers}-s{max_secondary}"


def build_autoscaling_config(num_workers, max_secondary):
    """Build an autoscaling policy YAML string for the given worker counts."""
    return f"""\
workerConfig:
    minInstances: {num_workers}
    maxInstances: {num_workers}
secondaryWorkerConfig:
    maxInstances: {max_secondary}
basicAlgorithm:
    cooldownPeriod: 120s
    yarnConfig:
        scaleUpFactor: 0.2
        scaleDownFactor: 0.5
        gracefulDecommissionTimeout: 1200s
"""


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
    return " ".join(string.split())


def _autoscaling_policy_exists(policy_name, workspace_project, region):
    """Return True if the named autoscaling policy already exists in Dataproc."""
    describe_cmd = unwrap(f"""
        gcloud dataproc autoscaling-policies describe {policy_name}
          --project={workspace_project}
          --region={region}
          --quiet
    """)
    wait_status = os.popen(f"{describe_cmd} > /dev/null 2>&1").close()
    return wait_status is None  # None means exit code 0 (policy found)


def create_autoscaling_policy(use_tiny_dataproc_cluster, workspace_project, region,
                               num_primary_workers=DEFAULT_NUM_WORKERS,
                               max_secondary_workers=DEFAULT_MAX_SECONDARY):
    """Create a GVS Spark autoscaling policy in Dataproc if it does not already exist.

    When use_tiny_dataproc_cluster is True the policy is fixed at
    TINY_NUM_WORKERS primary / TINY_MAX_SECONDARY secondary workers (suitable
    for integration tests).  Otherwise the caller-supplied num_primary_workers
    and max_secondary_workers values are used.

    Policy names encode their configuration (e.g. gvs-spark-autoscaling-w2-s200)
    so that a policy with a given name is written at most once and never
    overwritten with different content.

    Returns a tuple of (policy_name, num_workers).
    """
    if use_tiny_dataproc_cluster:
        num_workers = TINY_NUM_WORKERS
        max_secondary = TINY_MAX_SECONDARY
    else:
        num_workers = num_primary_workers
        max_secondary = max_secondary_workers

    policy_name = autoscaling_policy_name(num_workers, max_secondary)

    if _autoscaling_policy_exists(policy_name, workspace_project, region):
        info(f"Autoscaling policy '{policy_name}' already exists, skipping creation.")
        return policy_name, num_workers

    info(f"Creating autoscaling policy '{policy_name}' "
         f"(primary workers: {num_workers}, max secondary workers: {max_secondary})...")

    config = build_autoscaling_config(num_workers, max_secondary)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        f.write(config)
        yaml_path = f.name

    try:
        import_cmd = unwrap(f"""
            gcloud dataproc autoscaling-policies import {policy_name}
              --project={workspace_project}
              --source={yaml_path}
              --region={region}
              --quiet
        """)
        info(import_cmd)
        pipe = os.popen(import_cmd)
        info(pipe.read())
        wait_status = pipe.close()
        if wait_status:
            exit_code = os.waitstatus_to_exitcode(wait_status)
            raise RuntimeError(f"Unexpected exit code importing autoscaling policy: {exit_code}")
    finally:
        os.unlink(yaml_path)

    return policy_name, num_workers


def run_in_cluster(cluster_name, account, worker_machine_type, master_machine_type, region, use_tiny_dataproc_cluster, workspace_project,
                   script_path, secondary_script_path_list, script_arguments_json_path, leave_cluster_running_at_end, cluster_max_idle_minutes, cluster_max_age_minutes, master_memory_fraction,
                   num_primary_workers=DEFAULT_NUM_WORKERS, max_secondary_workers=DEFAULT_MAX_SECONDARY):

    cluster_max_idle_arg = f"--max-idle {cluster_max_idle_minutes}m" if cluster_max_idle_minutes else ""
    cluster_max_age_arg = f"--max-age {cluster_max_age_minutes}m" if cluster_max_age_minutes else ""

    try:
        autoscaling_policy, num_workers = create_autoscaling_policy(use_tiny_dataproc_cluster, workspace_project, region,
                                                                     num_primary_workers=num_primary_workers,
                                                                     max_secondary_workers=max_secondary_workers)

        cluster_start_cmd = unwrap(f"""
        
        hailctl dataproc start 
         --autoscaling-policy={autoscaling_policy}
         --num-workers {num_workers}
         --worker-machine-type {worker_machine_type}
         --master-machine-type {master_machine_type}
         --master-memory-fraction {master_memory_fraction}
         --enable-component-gateway
         --region {region}
         --project {workspace_project}
         --service-account {account}
         {cluster_max_idle_arg}
         {cluster_max_age_arg}
         --num-master-local-ssds 1
         --num-worker-local-ssds 1 
         --master-boot-disk-type=pd-ssd
         --worker-boot-disk-type=pd-ssd
         --subnet=projects/{workspace_project}/regions/{region}/subnetworks/subnetwork
         --properties=dataproc:dataproc.monitoring.stackdriver.enable=true,dataproc:dataproc.logging.stackdriver.enable=true,core:fs.gs.outputstream.sync.min.interval=5
         --packages=python-snappy
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

        run_in_existing_cluster(cluster_name, account, region, workspace_project,
                        script_path, secondary_script_path_list, script_arguments_json_path, leave_cluster_running_at_end)

    except Exception:
        raise

def run_in_existing_cluster(cluster_name, account, region, workspace_project,
                   script_path, secondary_script_path_list, script_arguments_json_path, leave_cluster_running_at_end):

    try:
        cluster_client = dataproc.ClusterControllerClient(
            client_options={"api_endpoint": f"{region}-dataproc.googleapis.com:443"}
        )

        # prepare custom arguments
        # the following says `--py-files` is supposed to be a comma separated list
        # https://fig.io/manual/gcloud/dataproc/jobs/submit/pyspark
        secondary_script_path_arg = f'--py-files {",".join(secondary_script_path_list)}' if secondary_script_path_list else ''
        with open(script_arguments_json_path, 'r') as input_file:
            items = ijson.items(input_file, '', use_float=True)
            arguments = items.__next__();
            custom_script_args = [f"--{key} {arguments.get(key)}" for key in arguments.keys()]

        found_cluster = False
        for cluster in cluster_client.list_clusters(request={"project_id": workspace_project, "region": region}):
            if cluster.cluster_name == cluster_name:
                found_cluster = True
                info("Found cluster: " + cluster_name)

                submit_cmd = unwrap(f"""

                gcloud dataproc jobs submit pyspark {script_path}
                 {secondary_script_path_arg}
                 --cluster={cluster_name}
                 --project {workspace_project}
                 --region={region}
                 --account {account}
                 --driver-log-levels root=WARN
                 --
                 {' '.join(custom_script_args)}
                """)

                info("Running: " + submit_cmd)
                pipe = os.popen(submit_cmd)
                pipe.read()
                wait_status = pipe.close()
                if wait_status:
                    exit_code = os.waitstatus_to_exitcode(wait_status)
                    raise RuntimeError(f"Unexpected exit code running submitted job: {exit_code}")
                break
        if not found_cluster:
            raise RuntimeError(f"Unable to find cluster: {cluster_name}")
    except Exception:
        raise

    finally:
        if leave_cluster_running_at_end:
            info(f"Leaving cluster {cluster_name} running as `leave_cluster_running_at_end` option is True.")
        else:
            info(f'Stopping cluster: {cluster_name}')
            delete_cmd = unwrap(f"""

                gcloud dataproc clusters delete
                  --project {workspace_project}
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
    parser.add_argument('--master-machine-type', type=str, required=False, default="n1-highmem-32",
                        help='Dataproc cluster master machine type')
    parser.add_argument('--master-memory-fraction', type=float, default=0.8, help='Dataproc master memory fraction')
    parser.add_argument('--region', type=str, required=True, help='GCS region')
    parser.add_argument('--use-tiny-dataproc-cluster', action='store_true', default=False,
                        help='Use a small autoscaling configuration suited for integration tests rather than '
                             'the default large configuration suited for production callsets.')
    parser.add_argument('--workspace-project', type=str, required=True, help='GCP project for the Dataproc cluster')
    parser.add_argument('--script-path', type=str, required=True, help='Path to script to run in Hail cluster')
    parser.add_argument('--secondary-script-path-list', type=str, required=False, action="append", default=[],
                        help='List of paths to secondary scripts to run in Hail cluster')
    parser.add_argument('--script-arguments-json-path', type=str, required=True,
                        help='JSON file of arguments for script')
    parser.add_argument('--leave-cluster-running-at-end', action="store_true", default=False)
    parser.add_argument('--cluster-max-idle-minutes', type=int, help='Maximum idle time of cluster in minutes')
    parser.add_argument('--cluster-max-age-minutes', type=int, help='Maximum age of cluster in minutes')
    parser.add_argument('--num-primary-workers', type=int, default=DEFAULT_NUM_WORKERS,
                        help=f'Number of primary workers for the non-tiny autoscaling policy (default: {DEFAULT_NUM_WORKERS}). '
                             f'Ignored when --use-tiny-dataproc-cluster is set.')
    parser.add_argument('--max-secondary-workers', type=int, default=DEFAULT_MAX_SECONDARY,
                        help=f'Maximum number of secondary workers for the non-tiny autoscaling policy (default: {DEFAULT_MAX_SECONDARY}). '
                             f'Ignored when --use-tiny-dataproc-cluster is set.')

    args = parser.parse_args()

    run_in_cluster(cluster_name=args.cluster_name,
                   account=args.account,
                   master_machine_type=args.master_machine_type,
                   worker_machine_type=args.worker_machine_type,
                   region=args.region,
                   use_tiny_dataproc_cluster=args.use_tiny_dataproc_cluster,
                   workspace_project=args.workspace_project,
                   script_path=args.script_path,
                   secondary_script_path_list=args.secondary_script_path_list,
                   script_arguments_json_path=args.script_arguments_json_path,
                   leave_cluster_running_at_end=args.leave_cluster_running_at_end,
                   cluster_max_idle_minutes=args.cluster_max_idle_minutes,
                   cluster_max_age_minutes=args.cluster_max_age_minutes,
                   master_memory_fraction=args.master_memory_fraction,
                   num_primary_workers=args.num_primary_workers,
                   max_secondary_workers=args.max_secondary_workers,
                   )
