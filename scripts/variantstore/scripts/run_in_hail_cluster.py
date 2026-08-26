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

# Substrings identifying a Compute Engine capacity shortage. These are properties of a
# single zone, not of the region: the same request commonly succeeds in a sibling zone. By
# default Dataproc places a cluster with Auto Zone placement, which picks one zone and
# fails outright rather than moving on, so a stockout there is a hard workflow failure.
# Deliberately narrow -- the enclosing error code is UNAVAILABLE, which also covers
# transient conditions that retrying in a different zone would not help.
STOCKOUT_MARKERS = (
    'STOCKOUT',
    'does not have enough resources',
    'ZONE_RESOURCE_POOL_EXHAUSTED',
)

# Capacity shortages that do not report themselves as such. When a zone can provide some
# but not all of the requested workers, Dataproc does not emit a STOCKOUT; the create
# simply times out waiting for nodes that never arrive, and the message blames firewall
# rules. Observed in us-central1-b immediately after a genuine STOCKOUT in us-central1-a,
# with only 1 of 4 requested workers provisioned.
#
# Retrying these in another zone is right when VM-to-VM networking is known to work, which
# it is for a Terra project that routinely runs Dataproc. If networking really is broken
# every zone will fail the same way, costing one attempt per zone before the error
# surfaces -- so the log distinguishes this case from an outright stockout.
PARTIAL_CAPACITY_MARKERS = (
    'Timed out waiting for',
    'minimum required datanodes',
    'Cannot start master',
)


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


def looks_like_stockout(output):
    """Whether command output indicates a zone ran out of capacity."""
    return any(marker in output for marker in STOCKOUT_MARKERS)


def looks_like_partial_capacity(output):
    """Whether output indicates only some requested nodes could be provisioned."""
    return any(marker in output for marker in PARTIAL_CAPACITY_MARKERS)


def retry_reason(output):
    """Why this failure is worth retrying in another zone, or None to fail fast.

    Returned as text rather than a boolean so the log records which of the two very
    different-looking capacity failures was seen.
    """
    if looks_like_stockout(output):
        return 'the zone is out of capacity'
    if looks_like_partial_capacity(output):
        return ('the zone could not provide all requested nodes (Dataproc reports this as '
                'a node timeout and blames firewall rules; if every zone fails this way, '
                'check VM-to-VM firewall rules for real)')
    return None


def zones_for_region(region, workspace_project):
    """List the zones of a region, for `--zones auto`.

    Derived from the region rather than hardcoded so this stays correct outside
    us-central1. On failure the caller falls back to Dataproc's own zone placement, since
    an inability to enumerate zones is not a reason to fail the workflow.
    """
    # The format and filter expressions are quoted because this runs through /bin/sh, where
    # the parentheses in value(name) are metacharacters. Unquoted they produce
    # `Syntax error: "(" unexpected`, the command never runs, and the caller silently falls
    # back to Dataproc's single-zone placement -- which is exactly the behavior --zones
    # exists to avoid.
    list_cmd = unwrap(f"""
        gcloud compute zones list
         --project={workspace_project}
         --filter='region:{region}'
         --format='value(name)'
         --quiet
    """)
    # stderr is folded in rather than discarded so a failure here is visible in the log.
    pipe = os.popen(list_cmd + " 2>&1")
    output = pipe.read()
    if pipe.close():
        info(f"Could not list zones in {region}; falling back to automatic zone placement. "
             f"Output: {output.strip()}")
        return []
    if not output.strip():
        info(f"No zones reported for {region}; falling back to automatic zone placement.")
        return []
    return sorted(zone.strip() for zone in output.split() if zone.strip())


def resolve_zones(zones, region, workspace_project):
    """Turn the --zones argument into an ordered list of zones to try.

    Empty means today's behavior: no --zone flag, and Dataproc places the cluster itself.
    """
    if not zones:
        return []
    if zones.strip().lower() == 'auto':
        return zones_for_region(region, workspace_project)
    return [zone.strip() for zone in zones.split(',') if zone.strip()]


def delete_failed_cluster(cluster_name, region, workspace_project):
    """Best-effort teardown so the cluster name is free for the next attempt.

    A create that fails on capacity usually leaves the cluster behind in ERROR state, and
    without this the retry would fail with ALREADY_EXISTS instead of surfacing the real
    problem.
    """
    delete_cmd = unwrap(f"""
        gcloud dataproc clusters delete {cluster_name}
         --region={region}
         --project={workspace_project}
         --quiet
    """)
    info(f"Removing failed cluster '{cluster_name}'.")
    pipe = os.popen(delete_cmd + " 2>&1")
    output = pipe.read()
    if pipe.close():
        info(f"Nothing to remove, or removal failed (continuing anyway): {output.strip()}")


def run_in_cluster(cluster_name, account, worker_machine_type, master_machine_type, region, use_tiny_dataproc_cluster, workspace_project,
                   script_path, secondary_script_path_list, script_arguments_json_path, leave_cluster_running_at_end, cluster_max_idle_minutes, cluster_max_age_minutes, master_memory_fraction,
                   num_primary_workers=DEFAULT_NUM_WORKERS, max_secondary_workers=DEFAULT_MAX_SECONDARY,
                   zones=None, num_local_ssds=1):

    cluster_max_idle_arg = f"--max-idle {cluster_max_idle_minutes}m" if cluster_max_idle_minutes else ""
    cluster_max_age_arg = f"--max-age {cluster_max_age_minutes}m" if cluster_max_age_minutes else ""

    try:
        autoscaling_policy, num_workers = create_autoscaling_policy(use_tiny_dataproc_cluster, workspace_project, region,
                                                                     num_primary_workers=num_primary_workers,
                                                                     max_secondary_workers=max_secondary_workers)

        # An empty list means one attempt with no --zone, i.e. Dataproc's own placement.
        candidate_zones = resolve_zones(zones, region, workspace_project) or [None]
        if len(candidate_zones) > 1:
            info(f"Will try up to {len(candidate_zones)} zone(s) in order: "
                 f"{', '.join(candidate_zones)}")

        for index, zone in enumerate(candidate_zones):
            zone_arg = f"--zone {zone}" if zone else ""
            placement = f"zone {zone}" if zone else "an automatically placed zone"

            cluster_start_cmd = unwrap(f"""

            hailctl dataproc start 
             --autoscaling-policy={autoscaling_policy}
             --num-workers {num_workers}
             --worker-machine-type {worker_machine_type}
             --master-machine-type {master_machine_type}
             --master-memory-fraction {master_memory_fraction}
             --enable-component-gateway
             --region {region}
             {zone_arg}
             --project {workspace_project}
             --service-account {account}
             {cluster_max_idle_arg}
             {cluster_max_age_arg}
             --num-master-local-ssds {num_local_ssds}
             --num-worker-local-ssds {num_local_ssds} 
             --master-boot-disk-type=pd-ssd
             --worker-boot-disk-type=pd-ssd
             --subnet=projects/{workspace_project}/regions/{region}/subnetworks/subnetwork
             --properties=dataproc:dataproc.monitoring.stackdriver.enable=true,dataproc:dataproc.logging.stackdriver.enable=true,core:fs.gs.outputstream.sync.min.interval=5
             --packages=python-snappy
             {cluster_name}
             
            """)

            info(f"Starting cluster '{cluster_name}' in {placement}...")
            info(cluster_start_cmd)
            # stderr is folded into stdout because gcloud reports capacity failures there,
            # and the retry decision depends on reading them.
            pipe = os.popen(cluster_start_cmd + " 2>&1")
            output = pipe.read()
            info(output)
            wait_status = pipe.close()
            if not wait_status:
                break

            exit_code = os.waitstatus_to_exitcode(wait_status)
            more_zones_to_try = index + 1 < len(candidate_zones)
            reason = retry_reason(output)
            retrying = more_zones_to_try and reason is not None

            # Always tear down after a failed create, whether or not another zone will be
            # tried. A cluster that fails to create is left behind in ERROR state with its
            # VMs running: it bills, and it holds regional CPU quota that the next attempt
            # -- or an unrelated cluster in the same project -- then cannot get.
            delete_failed_cluster(cluster_name, region, workspace_project)

            if retrying:
                info(f"Cluster creation failed in {placement}: {reason}. "
                     f"Retrying in {candidate_zones[index + 1]}.")
                continue
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
    parser.add_argument('--zones', type=str, required=False, default=None,
                        help='Comma-separated zones to try in order when creating the '
                             'cluster, or "auto" to use every zone in --region. A Compute '
                             'Engine stockout is a per-zone condition, so a create that '
                             'fails for capacity is retried in the next zone. Omit to keep '
                             "Dataproc's own placement, which picks one zone and does not "
                             'retry.')
    parser.add_argument('--num-local-ssds', type=int, required=False, default=1,
                        help='Local SSDs attached to the master and to each worker. Local '
                             'SSD availability is zone-specific and a common stockout '
                             'cause, so a workload that does not spill shuffle to disk can '
                             'set 0 to widen the usable zones. Defaults to 1.')

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
            zones=args.zones,
            num_local_ssds=args.num_local_ssds,
                   )
