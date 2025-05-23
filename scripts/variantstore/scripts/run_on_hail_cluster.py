import argparse
from logging import info

from run_in_hail_cluster import *

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


if __name__ == "__main__":
    configure_logging()

    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Get workspace information')

    parser.add_argument('--cluster-name', type=str, required=True, help='Name of the Hail cluster')
    parser.add_argument('--account', type=str, help='GCP account name')
    parser.add_argument('--region', type=str, required=True, help='GCS region')
    parser.add_argument('--gcs-project', type=str, required=True, help='GCS project')
    parser.add_argument('--script-path', type=str, required=True, help='Path to script to run in Hail cluster')
    parser.add_argument('--secondary-script-path-list', type=str, required=False, action="append", default=[],
                        help='List of paths to secondary scripts to run in Hail cluster')
    parser.add_argument('--script-arguments-json-path', type=str, required=True,
                        help='JSON file of arguments for script')

    args = parser.parse_args()

    run_on_existing_cluster(cluster_name=args.cluster_name,
                   account=args.account,
                   region=args.region,
                   gcs_project=args.gcs_project,
                   script_path=args.script_path,
                   secondary_script_path_list=args.secondary_script_path_list,
                   script_arguments_json_path=args.script_arguments_json_path,
                   leave_cluster_running_at_end=args.leave_cluster_running_at_end,
                   )
