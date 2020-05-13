import pytest
import pandas as pd

from ml4cvd.recipes import inference_file_name, hidden_inference_file_name
from ml4cvd.recipes import train_multimodal_multitask, compare_multimodal_multitask_models
from ml4cvd.recipes import infer_multimodal_multitask, infer_hidden_layer_multimodal_multitask
from ml4cvd.recipes import compare_multimodal_scalar_task_models, _find_learning_rate
# Imports with test in their name
from ml4cvd.recipes import test_multimodal_multitask as tst_multimodal_multitask
from ml4cvd.recipes import test_multimodal_scalar_tasks as tst_multimodal_scalar_tasks


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
