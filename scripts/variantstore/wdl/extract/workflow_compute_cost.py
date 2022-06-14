from collections import defaultdict
from firecloud import api as fapi
# import os
# workspace_namespace = os.environ['WORKSPACE_NAMESPACE']
# workspace_name = os.environ['WORKSPACE_NAME']

workspace_namespace = 'broad-firecloud-dsde'
workspace_name = 'VS-415 GVS Quickstart Default Extract Scatter'

submissions = fapi.list_submissions(workspace_namespace, workspace_name).json()
excluded_submission_ids = []

submission_ids = [s['submissionId'] for s in submissions]
filtered_submission_ids = [s for s in submission_ids if s not in excluded_submission_ids]

# Two-level default dictionary for workflow name -> workflow id -> cost
workflow_costs = defaultdict(lambda: defaultdict(dict))

for submission_id in filtered_submission_ids:
    submission = fapi.get_submission(workspace_namespace, workspace_name, submission_id).json()
    workflows = submission['workflows']
    if not workflows:
        print(f"Submission '{submission_id}' has no workflows, skipping cost calculation.")
        continue
    workflow_ids = [w['workflowId'] for w in workflows]
    for workflow_id in workflow_ids:
        workflow = fapi.get_workflow_metadata(workspace_namespace, workspace_name, submission_id, workflow_id).json()
        workflow_name = workflow['workflowName']
        if not workflow_name:
            print(f"Workflow {workflow_id} has no workflow name, skipping cost calculation.")
            continue
        print(f'Submission {submission_id}: workflow id {workflow_id}, name {workflow_name}')

        if workflow_name.startswith('Gvs'):
            if len(workflow_ids) == 1:
                cost = fapi.get_workflow_cost(workspace_namespace, workspace_name, submission_id, workflow_id).json()
                # If this run is < 1 day old the cost data may not yet be available.
                if 'cost' not in cost:
                    print(f"No cost found for workflow {workflow_id} in submission {submission_id}; possibly too recent.")
                    continue
                workflow_costs[workflow_name][workflow_id] = cost['cost']
                print(workflow_costs)
            else:
                print(f"Unexpectedly found GVS workflow '{workflow_name}' among {len(workflow_ids)} workflows in submission {submission_id}. Expecting only 1 workflow in submission, skipping cost calculation.")
                continue
        else:
            print(f"Workflow '{workflow_name}' in submission {submission_id} not recognized as a GVS workflow, skipping cost calculation.")
