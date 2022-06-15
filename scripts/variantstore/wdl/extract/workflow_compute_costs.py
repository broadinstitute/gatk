from logging import info, warning


def fapi_list_submissions(workspace_namespace: str, workspace_name: str):
    from firecloud import api as fapi
    return fapi.list_submissions(workspace_namespace, workspace_name).json()


def fapi_get_submission(workspace_namespace: str, workspace_name: str, submission_id: str):
    from firecloud import api as fapi
    return fapi.get_submission(workspace_namespace, workspace_name, submission_id).json()


def fapi_get_workflow_metadata(workspace_namespace: str, workspace_name: str, submission_id: str, workflow_id: str):
    from firecloud import api as fapi
    return fapi.get_workflow_metadata(workspace_namespace, workspace_name, submission_id, workflow_id).json()


def compute_costs(workspace_namespace, workspace_name, excluded_submission_ids,
                  list_submissions=fapi_list_submissions,
                  get_submission=fapi_get_submission,
                  get_workflow_metadata=fapi_get_workflow_metadata):

    submissions = list_submissions(workspace_namespace, workspace_name)
    submission_ids = [s['submissionId'] for s in submissions]

    # Two-level default dictionary for workflow name -> workflow id -> cost
    from collections import defaultdict
    workflow_costs = defaultdict(lambda: defaultdict(dict))

    for submission_id in submission_ids:
        if excluded_submission_ids and submission_id in excluded_submission_ids:
            info(f"Submission id '{submission_id}' in exclude list, skipping cost calculation.")
            continue
        submission = get_submission(workspace_namespace, workspace_name, submission_id)
        workflows = submission['workflows']
        if not workflows:
            print(f"Submission '{submission_id}' has no workflows, skipping cost calculation.")
            continue
        workflow_ids = [w['workflowId'] for w in workflows]
        for workflow_id in workflow_ids:
            workflow = get_workflow_metadata(workspace_namespace, workspace_name, submission_id, workflow_id)
            workflow_name = workflow['workflowName']
            if not workflow_name:
                print(f"Workflow {workflow_id} has no workflow name, skipping cost calculation.")
                continue
            # info(f'Submission {submission_id}: workflow id {workflow_id}, name {workflow_name}')

            if workflow_name.startswith('Gvs'):
                if len(workflow_ids) == 1:
                    # If this run is < 1 day old the cost data may not yet be available.
                    if 'cost' not in submission:
                        print(f"No cost found for workflow {workflow_id} in submission {submission_id}; possibly too recent.")
                        continue
                    workflow_costs[workflow_name][workflow_id] = submission['cost']
                else:
                    print(f"Unexpectedly found GVS workflow '{workflow_name}' among {len(workflow_ids)} workflows in submission {submission_id}. Expecting only 1 workflow in submission, skipping cost calculation.")
                    continue
            else:
                print(f"Workflow '{workflow_name}' in submission {submission_id} not recognized as a GVS workflow, skipping cost calculation.")

    return workflow_costs


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Populate an alt_allele table for the BigQuery Variant Store')

    parser.add_argument('--workspace_namespace', type=str, help='Namespace for workspace', required=True)
    parser.add_argument('--workspace_name', type=str, help='Name of workspace', required=True)
    parser.add_argument('--exclude', type=str, help='Submission ID to exclude from cost calculations',
                        required=False, action="append")

    return parser.parse_args()


def configure_logging():
    import logging
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # create formatter and add it to the handler
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    # add the handler to the logger
    logger.addHandler(ch)


if __name__ == '__main__':
    configure_logging()
    args = parse_args()
    costs = compute_costs(args.workspace_namespace, args.workspace_name, args.exclude)

    wdls = [
        'GvsAssignIds',
        'GvsImportGenomes',
        'GvsCreateAltAllele',
        'GvsCreateFilterSet',
        'GvsPrepareCallset',
        'GvsExtractCallset',
        'GvsUnified',
        'GvsQuickstartIntegration'
    ]

    print("workflow_name\tworkflow_id\tcost")
    for wdl in wdls:
        if wdl in costs:
            for workflow_id, cost in costs[wdl].items():
                print(f"{wdl}\t{workflow_id}\t{cost}")
