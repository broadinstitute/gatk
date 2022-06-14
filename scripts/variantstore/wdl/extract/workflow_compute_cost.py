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

core_workflows = [
    'GvsAssignIds',
    'GvsImportGenomes',
    'GvsCreateAltAllele',
    'GvsCreateFilterSet',
    'GvsPrepareCallset',
    'GvsExtractCallset'
]

core_workflow_wrappers = [
    'GvsUnified',
    'GvsQuickstartIntegration'
]

workflow_costs = defaultdict(lambda: defaultdict(dict))


def assign_workflow_costs_from_parent_workflow(parent):
    # print(parent)
    for fully_qualified_name, subworkflows in parent['calls'].items():
        subworkflow_id = subworkflows[0]['subWorkflowId']
        subworkflow_name = fully_qualified_name.split('.')[-1]
        subworkflow_cost = fapi.get_workflow_cost(workspace_namespace, workspace_name, submission_id, subworkflow_id).json()
        print(f"cost is {subworkflow_cost}")
        workflow_costs[subworkflow_name][subworkflow_id] = subworkflow_cost
        print(f"Assigned cost of ${subworkflow_cost} to '{subworkflow_name}' id {subworkflow_id}")


for submission_id in filtered_submission_ids:
    submission = fapi.get_submission(workspace_namespace, workspace_name, submission_id).json()
    workflows = submission['workflows']
    if not workflows:
        continue
    workflow_ids = [w['workflowId'] for w in workflows]
    for workflow_id in workflow_ids:
        workflow = fapi.get_workflow_metadata(workspace_namespace, workspace_name, submission_id, workflow_id).json()
        workflow_name = workflow['workflowName']
        if not workflow_name:
            continue
        print(f'Submission {submission_id}: workflow id {workflow_id}, name {workflow_name}')

        if workflow_name in core_workflows:
            cost = fapi.get_workflow_cost(workspace_namespace, workspace_name, submission_id, workflow_id)
            workflow_costs[workflow_name][workflow_id] = cost
        elif workflow_name == 'GvsUnified':
            assign_workflow_costs_from_parent_workflow(workflow)
        elif workflow_name == 'GvsQuickstartIntegration':
            unified_id = workflow['calls']['GvsQuickstartIntegration.GvsUnified'][0]['subWorkflowId']
            unified = fapi.get_workflow_metadata(workspace_namespace, workspace_name, submission_id, unified_id).json()
            assign_workflow_costs_from_parent_workflow(unified)
        else:
            print(f"Workflow '{workflow_name}' not recognized as a GVS workflow, skipping cost calculations.")

    print(workflow_costs)
