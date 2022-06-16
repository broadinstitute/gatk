import unittest

from workflow_compute_costs import compute_costs
from workflow_compute_costs_test_helper import *

WORKSPACE_NAMESPACE = 'test_workspace_namespace'
WORKSPACE_NAME = 'test_workspace'


def mock_list_submissions(_workspace_namespace, _workspace_name):
    return [
        FAILED_INTEGRATION_SUBMISSION,
        SUCCEEDED_ASSIGN_IDS_SUBMISSION,
        SUCCEEDED_INTEGRATION_SUBMISSION
    ]


def mock_get_submission(_workspace_namespace, _workspace_name, submission_id):
    if submission_id == FAILED_INTEGRATION_SUBMISSION_ID:
        return FAILED_INTEGRATION_SUBMISSION
    elif submission_id == SUCCEEDED_INTEGRATION_SUBMISSION_ID:
        return SUCCEEDED_INTEGRATION_SUBMISSION
    elif submission_id == SUCCEEDED_ASSIGN_IDS_SUBMISSION_ID:
        return SUCCEEDED_ASSIGN_IDS_SUBMISSION
    else:
        raise ValueError(f"Unrecognized submission id '{submission_id}")


def mock_get_workflow_metadata(_workspace_namespace, _workspace_name, submission_id, workflow_id):
    if submission_id == FAILED_INTEGRATION_SUBMISSION_ID:
        expected_workflow = FAILED_INTEGRATION_WORKFLOW
    elif submission_id == SUCCEEDED_INTEGRATION_SUBMISSION_ID:
        expected_workflow = SUCCEEDED_INTEGRATION_WORKFLOW
    elif submission_id == SUCCEEDED_ASSIGN_IDS_SUBMISSION_ID:
        expected_workflow = SUCCEEDED_ASSIGN_IDS_WORKFLOW
    else:
        raise ValueError(f"Unrecognized submission id '{submission_id}")

    expected_workflow_id = expected_workflow['id']
    if workflow_id != expected_workflow_id:
        raise ValueError(
            f"Unexpected: workflow id '{workflow_id}' was expected to be '{expected_workflow_id}' for submission '{submission_id}'")

    return expected_workflow


class TestWorkflowComputeCosts(unittest.TestCase):
    def test_empty_submission_list(self):
        costs = compute_costs(workspace_namespace=WORKSPACE_NAMESPACE,
                              workspace_name=WORKSPACE_NAME,
                              excluded_submission_ids=[],
                              list_submissions=lambda _ns, _n: [],  # return an empty list of submissions
                              get_submission=None,
                              get_workflow_metadata=None
                              )
        # Asserts both empty costs and the non-invocation of `get_submission` and `get_workflow_metadata`.
        self.assertEqual(len(costs), 0, msg="Costs should be empty when there are no submissions in the workspace.")

    def check_totals(self, costs):
        expected_links = [w['link'] for w in costs['workflows']]
        self.assertEqual(costs['total']['links'], expected_links)
        self.assertAlmostEqual(sum([w['cost'] for w in costs['workflows']]), costs['total']['cost'])

    def test_all_submissions(self):
        self.maxDiff = None
        costs = compute_costs(workspace_namespace=WORKSPACE_NAMESPACE,
                              workspace_name=WORKSPACE_NAME,
                              excluded_submission_ids=[],
                              list_submissions=mock_list_submissions,
                              get_submission=mock_get_submission,
                              get_workflow_metadata=mock_get_workflow_metadata
                              )

        self.assertEqual(len(costs), 2, msg="Expecting exactly two workflow names in costs dictionary")
        integration = costs['GvsQuickstartIntegration']

        expected = [{'cost': 0.885638,
                     'link': 'https://app.terra.bio/#workspaces/test_workspace_namespace/test_workspace/job_history/73ac71db-0488-46be-a8e8-7f00e795edec',
                     'submission_id': '73ac71db-0488-46be-a8e8-7f00e795edec',
                     'submission_timestamp': '2022-06-10T10:22:34.320Z',
                     'workflow_id': 'a7e9cf65-f3e6-4ede-b64e-38018ca560f3'},
                    {'cost': 0.104732,
                     'link': 'https://app.terra.bio/#workspaces/test_workspace_namespace/test_workspace/job_history/136781a2-64b8-4ede-8652-973136030aaf',
                     'submission_id': '136781a2-64b8-4ede-8652-973136030aaf',
                     'submission_timestamp': '2022-06-09T20:28:57.797Z',
                     'workflow_id': '902df806-5852-43d0-9816-ff46bf7e1716'},
                    ]
        self.assertEqual(integration['workflows'], expected)
        self.check_totals(integration)

        expected = [{'cost': 0.018258,
                     'link': 'https://app.terra.bio/#workspaces/test_workspace_namespace/test_workspace/job_history/7b74f3e2-82ff-4289-92a5-5a801025a0f5',
                     'submission_id': '7b74f3e2-82ff-4289-92a5-5a801025a0f5',
                     'submission_timestamp': '2022-06-13T14:43:57.830Z',
                     'workflow_id': 'e3120691-cf45-471d-85d8-a10a8ef51a07'}]
        assign_ids = costs['GvsAssignIds']
        self.assertEqual(assign_ids['workflows'], expected)
        self.check_totals(assign_ids)

    def test_apply_exclusion(self):
        self.maxDiff = None
        excluded_submission_id = '136781a2-64b8-4ede-8652-973136030aaf'
        costs = compute_costs(workspace_namespace=WORKSPACE_NAMESPACE,
                              workspace_name=WORKSPACE_NAME,
                              excluded_submission_ids=[excluded_submission_id],
                              list_submissions=mock_list_submissions,
                              get_submission=mock_get_submission,
                              get_workflow_metadata=mock_get_workflow_metadata
                              )

        self.assertEqual(len(costs), 2, msg="Expecting exactly two kinds of workflows in costs dictionary")

        integration = costs['GvsQuickstartIntegration']
        expected = [{'cost': 0.885638,
                     'link': 'https://app.terra.bio/#workspaces/test_workspace_namespace/test_workspace/job_history/73ac71db-0488-46be-a8e8-7f00e795edec',
                     'submission_id': '73ac71db-0488-46be-a8e8-7f00e795edec',
                     'submission_timestamp': '2022-06-10T10:22:34.320Z',
                     'workflow_id': 'a7e9cf65-f3e6-4ede-b64e-38018ca560f3'}]
        self.assertEqual(integration['workflows'], expected)
        self.check_totals(integration)

        expected = [{'cost': 0.018258,
                     'link': 'https://app.terra.bio/#workspaces/test_workspace_namespace/test_workspace/job_history/7b74f3e2-82ff-4289-92a5-5a801025a0f5',
                     'submission_id': '7b74f3e2-82ff-4289-92a5-5a801025a0f5',
                     'submission_timestamp': '2022-06-13T14:43:57.830Z',
                     'workflow_id': 'e3120691-cf45-471d-85d8-a10a8ef51a07'}]
        assign_ids = costs['GvsAssignIds']
        self.assertEqual(assign_ids['workflows'], expected)
        self.check_totals(assign_ids)
