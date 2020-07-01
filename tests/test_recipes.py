import os
import pytest
import pandas as pd
import numpy as np

from ml4cvd.recipes import inference_file_name, hidden_inference_file_name
from ml4cvd.recipes import train_multimodal_multitask, compare_multimodal_multitask_models
from ml4cvd.recipes import infer_multimodal_multitask, infer_hidden_layer_multimodal_multitask
from ml4cvd.recipes import compare_multimodal_scalar_task_models, _find_learning_rate
from ml4cvd.explorations import _continuous_explore_header, _categorical_explore_header, _should_error_detect, explore
# Imports with test in their name
from ml4cvd.recipes import test_multimodal_multitask as tst_multimodal_multitask
from ml4cvd.recipes import test_multimodal_scalar_tasks as tst_multimodal_scalar_tasks
from ml4cvd.test_utils import TMAPS_UP_TO_4D
from ml4cvd.test_utils import build_hdf5s
from ml4cvd.TensorMap import TensorMap, Interpretation


class TestRecipes:
    def test_train(self, default_arguments):
        train_multimodal_multitask(default_arguments)

    def test_test(self, default_arguments):
        tst_multimodal_multitask(default_arguments)

    def test_test_scalar(self, default_arguments):
        tst_multimodal_scalar_tasks(default_arguments)

    def test_infer(self, default_arguments):
        infer_multimodal_multitask(default_arguments)
        tsv = inference_file_name(default_arguments.output_folder, default_arguments.id)
        inferred = pd.read_csv(tsv, sep='\t')
        assert len(set(inferred['sample_id'])) == pytest.N_TENSORS

    def test_infer_genetics(self, default_arguments):
        default_arguments.tsv_style = 'genetics'
        infer_multimodal_multitask(default_arguments)
        default_arguments.tsv_style = 'standard'
        tsv = inference_file_name(default_arguments.output_folder, default_arguments.id)
        inferred = pd.read_csv(tsv, sep='\t')
        assert len(set(inferred['FID'])) == pytest.N_TENSORS

    def test_infer_hidden(self, default_arguments):
        infer_hidden_layer_multimodal_multitask(default_arguments)
        tsv = hidden_inference_file_name(default_arguments.output_folder, default_arguments.id)
        inferred = pd.read_csv(tsv, sep='\t')
        assert len(set(inferred['sample_id'])) == pytest.N_TENSORS

    def test_infer_hidden_genetics(self, default_arguments):
        default_arguments.tsv_style = 'genetics'
        infer_hidden_layer_multimodal_multitask(default_arguments)
        default_arguments.tsv_style = 'standard'
        tsv = hidden_inference_file_name(default_arguments.output_folder, default_arguments.id)
        inferred = pd.read_csv(tsv, sep='\t')
        assert len(set(inferred['FID'])) == pytest.N_TENSORS

    def test_find_learning_rate(self, default_arguments):
        _find_learning_rate(default_arguments)

    def test_explore(self, default_arguments, tmpdir_factory):
        temp_dir = tmpdir_factory.mktemp('explore_tensors')
        default_arguments.tensors = str(temp_dir)
        tmaps = TMAPS_UP_TO_4D[:]
        tmaps.append(TensorMap(f'scalar', shape=(1,), interpretation=Interpretation.CONTINUOUS))
        explore_expected = build_hdf5s(temp_dir, tmaps, n=pytest.N_TENSORS)
        default_arguments.num_workers = 3
        default_arguments.tensor_maps_in = tmaps
        explore(default_arguments)
        csv_path = os.path.join(
            default_arguments.output_folder, default_arguments.id, 'tensors_all_union.csv'
        )
        explore_result = pd.read_csv(csv_path)
        for row in explore_result.iterrows():
            row = row[1]
            for tm in tmaps:
                row_expected = explore_expected[(row['fpath'], tm)]
                if _should_error_detect(tm):
                    actual = getattr(row, _continuous_explore_header(tm))
                    assert not np.isnan(actual)
                    continue
                if tm.is_continuous():
                    actual = getattr(row, _continuous_explore_header(tm))
                    assert actual == row_expected
                    continue
                if tm.is_categorical():
                    for channel, idx in tm.channel_map.items():
                        channel_val = getattr(row, _categorical_explore_header(tm, channel))
                        assert channel_val == row_expected[idx]
