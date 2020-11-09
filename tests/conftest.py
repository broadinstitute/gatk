import sys
import pytest

from ml4h.arguments import parse_args
from ml4h.test_utils import TMAPS as MOCK_TMAPS
from ml4h.test_utils import build_hdf5s


def pytest_configure(config):
    pytest.N_TENSORS = 50
    config.addinivalue_line("markers", "slow: mark tests as slow")


@pytest.fixture(scope='class')
def default_arguments(tmpdir_factory):
    temp_dir = tmpdir_factory.mktemp('data')
    build_hdf5s(temp_dir, MOCK_TMAPS.values(), n=pytest.N_TENSORS)
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
