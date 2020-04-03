import sys
import mock
import pytest

from ml4cvd.arguments import parse_args, TMAPS
from ml4cvd.test_utils import TMAPS as MOCK_TMAPS
from ml4cvd.test_utils import build_hdf5s
from ml4cvd.recipes import train_multimodal_multitask, compare_multimodal_multitask_models
from ml4cvd.recipes import infer_multimodal_multitask, infer_hidden_layer_multimodal_multitask
from ml4cvd.recipes import compare_multimodal_scalar_task_models, _find_learning_rate
# Imports with test in their name
from ml4cvd.recipes import test_multimodal_multitask as tst_multimodal_multitask
from ml4cvd.recipes import test_multimodal_scalar_tasks as tst_multimodal_scalar_tasks


@pytest.fixture(scope='class')
@mock.patch.dict(TMAPS, MOCK_TMAPS)
def default_arguments(tmpdir_factory):
    temp_dir = tmpdir_factory.mktemp('data')
    build_hdf5s(temp_dir, MOCK_TMAPS.values(), n=50)
    hdf5_dir = str(temp_dir)
    inp_key = '3d_cont'
    out_key = '1d_cat'
    sys.argv = [
        '',
        '--output_folder', hdf5_dir,
        '--input_tensors', inp_key,
        '--output_tensors', out_key,
        '--tensors', hdf5_dir,
        '--pool_x', '1',
        '--pool_y', '1',
        '--pool_z', '1',
        '--training_steps', '2',
        '--test_steps', '3',
        '--validation_steps', '2',
        '--epochs', '2',
        '--num_workers', '0',
        '--batch_size', '2',
    ]
    args = parse_args()
    return args


class TestRecipes:
    """Smoke tests"""

    def test_train(self, default_arguments):
        train_multimodal_multitask(default_arguments)

    def test_test(self, default_arguments):
        tst_multimodal_multitask(default_arguments)

    def test_test_scalar(self, default_arguments):
        tst_multimodal_scalar_tasks(default_arguments)

    def test_infer(self, default_arguments):
        infer_multimodal_multitask(default_arguments)

    def test_infer_hidden(self, default_arguments):
        infer_hidden_layer_multimodal_multitask(default_arguments)

    def test_find_learning_rate(self, default_arguments):
        _find_learning_rate(default_arguments)
