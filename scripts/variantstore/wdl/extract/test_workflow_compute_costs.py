import unittest

from workflow_compute_costs import compute_costs
from workflow_compute_costs_test_helper import *

WORKSPACE_NAMESPACE='test_workspace_namespace'
WORKSPACE_NAME='test_workspace'

FAILED_SUBMISSION_ID='136781a2-64b8-4ede-8652-973136030aaf'
EXCLUDED_SUBMISSION_ID='3b4ca454-bf05-4df0-806d-20880e1c7cfe'
SUCCEEDED_SUBMISSION_ID=''

SUBMISSION_IDS_RESPONSE={

}


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
        self.assertTrue(len(costs) == 0)

    def test_return_all_costs(self):
        pass

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
