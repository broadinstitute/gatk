import unittest

from workflow_compute_costs import compute_costs
from workflow_compute_costs_test_helper import *

WORKSPACE_NAMESPACE='test_workspace_namespace'
WORKSPACE_NAME='test_workspace'


def list_submissions(_workspace_namespace, _workspace_name):
    return [
        FAILED_INTEGRATION_SUBMISSION,
        SUCCEEDED_ASSIGN_IDS_SUBMISSION,
        SUCCEEDED_INTEGRATION_SUBMISSION
    ]


def get_submission(_workspace_namespace, _workspace_name, submission_id):
    if submission_id == FAILED_INTEGRATION_SUBMISSION_ID:
        return FAILED_INTEGRATION_SUBMISSION
    elif submission_id == SUCCEEDED_INTEGRATION_SUBMISSION_ID:
        return SUCCEEDED_INTEGRATION_SUBMISSION
    elif submission_id == SUCCEEDED_ASSIGN_IDS_SUBMISSION_ID:
        return SUCCEEDED_ASSIGN_IDS_SUBMISSION
    else:
        raise ValueError(f"Unrecognized submission id '{submission_id}")


def get_workflow_metadata(_workspace_namespace, _workspace_name, submission_id, workflow_id):
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
        raise ValueError(f"Unexpected: workflow id '{workflow_id}' was expected to be '{expected_workflow_id}' for submission '{submission_id}'")

    return expected_workflow


class TestWorkflowComputeCosts(unittest.TestCase):
    def test_empty_submission_list(self):
        costs = compute_costs(workspace_namespace=WORKSPACE_NAMESPACE,
                              workspace_name=WORKSPACE_NAME,
                              excluded_submission_ids=[],
                              list_submissions=lambda _ns, _n: [],
                              get_submission=None,
                              get_workflow_metadata=None
                              )
        # Asserts both empty costs and the non-invocation of `get_submission` and `get_workflow_metadata`.
        self.assertEqual(len(costs), 0)

    def test_all_submissions(self):
        costs = compute_costs(workspace_namespace=WORKSPACE_NAMESPACE,
                              workspace_name=WORKSPACE_NAME,
                              excluded_submission_ids=[],
                              list_submissions=list_submissions,
                              get_submission=get_submission,
                              get_workflow_metadata=get_workflow_metadata
                              )
        # Asserts both empty costs and the non-invocation of `get_submission` and `get_workflow_metadata`.
        self.assertEqual(len(costs), 2)

    def test_apply_exclusion(self):
        pass



    # def test_downsampled_scaling(self):
    #     import filecmp
    #     import tempfile
    #
    #     for x, y in [(10, 10), (1, 10), (10, 1)]:
    #         with tempfile.NamedTemporaryFile() as actual_output_bed:
    #             scale_xy_bed_values('scale_xy_bed_values_test_files/intervals_downsampled_5.bed',
    #                                 actual_output_bed.name,
    #                                 x,
    #                                 y)
    #
    #             expected_output_bed = f'scale_xy_bed_values_test_files/intervals_downsampled_5_scaled_{x}_{y}.bed'
    #             self.assertTrue(filecmp.cmp(expected_output_bed,
    #                                         actual_output_bed.name,
    #                                         shallow=False), f'fail on X scaling {x} and Y scaling {y}')
