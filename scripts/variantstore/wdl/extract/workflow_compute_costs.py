
def fapi_list_submissions(workspace_namespace: str, workspace_name: str, submission_id: str):
    from firecloud import api as fapi
    return fapi.mock_list_submissions(workspace_namespace, workspace_name, submission_id).json()


def fapi_get_submission(workspace_namespace: str, workspace_name: str, submission_id: str):
    from firecloud import api as fapi
    return fapi.mock_get_submission(workspace_namespace, workspace_name, submission_id).json()


def fapi_get_workflow_metadata(workspace_namespace: str, workspace_name: str, submission_id: str, workflow_id: str):
    from firecloud import api as fapi
    return fapi.mock_get_workflow_metadata(workspace_namespace, workspace_name, submission_id, workflow_id).json()


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
        if submission_id in excluded_submission_ids:
            print(f"Submission id '{submission_id}' in exclude list, skipping cost calculation.")
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
            print(f'Submission {submission_id}: workflow id {workflow_id}, name {workflow_name}')

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


if __name__ == '__main__':
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

    for wdl in wdls:
        if wdl in costs:
            cost = sum(costs[wdl].values())
            print(f"{wdl}: {len(costs[wdl])} workflows, total compute cost ${cost:.2f}")
