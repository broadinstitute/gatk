import sys
import pytest
from ml4h.arguments import parse_args
from ml4h.test_utils import TMAPS as MOCK_TMAPS


class TestUConnect:

    def test_no_u(self, tmpdir):
        sys.argv = [
            'train',
            '--output_folder', str(tmpdir),
        ]
        args = parse_args()
        assert len(args.u_connect) == 0

    def test_simple_u(self, tmpdir):
        inp_key = '3d_cont'
        sys.argv = [
            'train',
            '--output_folder', str(tmpdir),
            '--input_tensors', inp_key,
            '--output_tensors', inp_key,
            '--u_connect', inp_key, inp_key,
        ]
        args = parse_args()
        assert len(args.u_connect) == 1
        inp, out = list(args.u_connect.items())[0]
        tmap = MOCK_TMAPS[inp_key]
        assert inp == tmap
        assert out == {tmap, }

    def test_many_to_one(self, tmpdir):
        inp_key1 = '3d_cont'
        inp_key2 = '3d_cat'
        sys.argv = [
            'train',
            '--output_folder', str(tmpdir),
            '--input_tensors', inp_key1, inp_key2,
            '--output_tensors', inp_key1,
            '--u_connect', inp_key1, inp_key1,
            '--u_connect', inp_key2, inp_key1,
        ]
        args = parse_args()
        assert len(args.u_connect) == 2
        assert args.u_connect[MOCK_TMAPS[inp_key1]] == {MOCK_TMAPS[inp_key1]}
        assert args.u_connect[MOCK_TMAPS[inp_key2]] == {MOCK_TMAPS[inp_key1]}

    def test_one_to_many(self, tmpdir):
        key1 = '3d_cont'
        key2 = '3d_cat'
        sys.argv = [
            'train',
            '--output_folder', str(tmpdir),
            '--input_tensors', key1, key2,
            '--output_tensors', key1, key2,
            '--u_connect', key1, key1,
            '--u_connect', key1, key2,
        ]
        args = parse_args()
        assert len(args.u_connect) == 1
        assert args.u_connect[MOCK_TMAPS[key1]] == {MOCK_TMAPS[key1], MOCK_TMAPS[key2]}

    def test_multi_u(self, tmpdir):
        key1 = '3d_cont'
        key2 = '3d_cat'
        sys.argv = [
            'train',
            '--output_folder', str(tmpdir),
            '--input_tensors', key1, key2,
            '--output_tensors', key1, key2,
            '--u_connect', key1, key1,
            '--u_connect', key2, key2,
        ]
        args = parse_args()
        assert len(args.u_connect) == 2
        assert args.u_connect[MOCK_TMAPS[key1]] == {MOCK_TMAPS[key1], }
        assert args.u_connect[MOCK_TMAPS[key2]] == {MOCK_TMAPS[key2], }
