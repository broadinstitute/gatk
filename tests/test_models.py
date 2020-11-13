import os
import pytest
import numpy as np
import tensorflow as tf
from itertools import cycle
from collections import defaultdict
from typing import List, Optional, Dict, Tuple, Iterator

from ml4h.TensorMap import TensorMap
from ml4h.models import make_multimodal_multitask_model, parent_sort, BottleneckType, ACTIVATION_FUNCTIONS, MODEL_EXT, train_model_from_generators, \
    check_no_bottleneck, make_paired_autoencoder_model
from ml4h.test_utils import TMAPS_UP_TO_4D, MULTIMODAL_UP_TO_4D, CATEGORICAL_TMAPS, CONTINUOUS_TMAPS, SEGMENT_IN, SEGMENT_OUT, PARENT_TMAPS, CYCLE_PARENTS
from ml4h.test_utils import LANGUAGE_TMAP_1HOT_WINDOW, LANGUAGE_TMAP_1HOT_SOFTMAX


MEAN_PRECISION_EPS = .02  # how much mean precision degradation is acceptable
DEFAULT_PARAMS = {
    'activation': 'relu',
    'dense_layers': [4, 2],
    'dense_blocks': [5, 3],
    'block_size': 3,
    'learning_rate': 1e-3,
    'optimizer': 'adam',
    'conv_type': 'conv',
    'conv_layers': [6, 5, 3],
    'conv_width': [71]*5,
    'conv_x': [3]*5,
    'conv_y': [3]*5,
    'conv_z': [2]*5,
    'padding': 'same',
    'max_pools': [],
    'pool_type': 'max',
    'pool_x': 1,
    'pool_y': 1,
    'pool_z': 1,
    'conv_regularize': 'spatial_dropout',
    'conv_regularize_rate': .1,
    'conv_normalize': 'batch_norm',
    'dense_regularize': 'dropout',
    'dense_regularize_rate': .1,
    'dense_normalize': 'batch_norm',
    'bottleneck_type': BottleneckType.FlattenRestructure,
    'pair_loss': 'cosine',
    'training_steps': 12,
    'learning_rate': 0.00001,
    'epochs': 6,
    'optimizer': 'adam',
    'learning_rate_schedule': None,
    'model_layers': None,
    'model_file': None,
    'hidden_layer': 'embed',
    'u_connect': defaultdict(dict),
}


TrainType = Dict[str, np.ndarray]  # TODO: better name


def make_training_data(input_tmaps: List[TensorMap], output_tmaps: List[TensorMap]) -> Iterator[Tuple[TrainType, TrainType, List[None]]]:
    return cycle([
        (
            {tm.input_name(): tf.random.normal((2,) + tm.shape) for tm in input_tmaps},
            {tm.output_name(): tf.zeros((2,) + tm.shape) for tm in output_tmaps},
            [None] * len(output_tmaps),
        ), ])


def assert_model_trains(input_tmaps: List[TensorMap], output_tmaps: List[TensorMap], m: Optional[tf.keras.Model] = None, skip_shape_check: bool = False):
    if m is None:
        m = make_multimodal_multitask_model(
            input_tmaps,
            output_tmaps,
            **DEFAULT_PARAMS,
        )
    if not skip_shape_check:
        for tmap, tensor in zip(input_tmaps, m.inputs):
            assert tensor.shape[1:] == tmap.shape
            assert tensor.shape[1:] == tmap.shape
        for tmap, tensor in zip(parent_sort(output_tmaps), m.outputs):
            assert tensor.shape[1:] == tmap.shape
            assert tensor.shape[1:] == tmap.shape
    data = make_training_data(input_tmaps, output_tmaps)
    history = m.fit(data, steps_per_epoch=2, epochs=2, validation_data=data, validation_steps=2)
    for tmap in output_tmaps:
        for metric in tmap.metrics:
            metric_name = metric if type(metric) == str else metric.__name__
            name = f'{tmap.output_name()}_{metric_name}' if len(output_tmaps) > 1 else metric_name
            assert name in history.history


def _rotate(a: List, n: int):
    return a[-n:] + a[:-n]


class TestMakeMultimodalMultitaskModel:
    @pytest.mark.parametrize(
        'input_output_tmaps',
        [
            (CONTINUOUS_TMAPS[:1], CONTINUOUS_TMAPS[1:2]), (CONTINUOUS_TMAPS[1:2], CONTINUOUS_TMAPS[:1]),
            (CONTINUOUS_TMAPS[:2], CONTINUOUS_TMAPS[:2]),
        ],
    )
    def test_multimodal_multitask_quickly(self, input_output_tmaps):
        """
        Tests 1d->2d, 2d->1d, (1d,2d)->(1d,2d)
        """
        assert_model_trains(input_output_tmaps[0], input_output_tmaps[1])

    @pytest.mark.slow
    @pytest.mark.parametrize(
        'input_tmaps',
        MULTIMODAL_UP_TO_4D,
    )
    @pytest.mark.parametrize(
        'output_tmaps',
        MULTIMODAL_UP_TO_4D,
    )
    def test_multimodal(self, input_tmaps: List[TensorMap], output_tmaps: List[TensorMap]):
        assert_model_trains(input_tmaps, output_tmaps)

    @pytest.mark.slow
    @pytest.mark.parametrize(
        'input_tmap',
        CONTINUOUS_TMAPS[:-1],
    )
    @pytest.mark.parametrize(
        'output_tmap',
        TMAPS_UP_TO_4D,
    )
    def test_unimodal_md_to_nd(self, input_tmap: TensorMap, output_tmap: TensorMap):
        assert_model_trains([input_tmap], [output_tmap])

    @pytest.mark.slow
    @pytest.mark.parametrize(
        'input_tmap',
        TMAPS_UP_TO_4D,
    )
    @pytest.mark.parametrize(
        'output_tmap',
        TMAPS_UP_TO_4D,
    )
    def test_load_unimodal(self, tmpdir, input_tmap, output_tmap):
        m = make_multimodal_multitask_model(
            [input_tmap],
            [output_tmap],
            **DEFAULT_PARAMS,
        )
        path = os.path.join(tmpdir, f'm{MODEL_EXT}')
        m.save(path)
        make_multimodal_multitask_model(
            [input_tmap],
            [output_tmap],
            model_file=path,
            **DEFAULT_PARAMS,
        )

    @pytest.mark.slow
    @pytest.mark.parametrize(
        'activation',
        ACTIVATION_FUNCTIONS.keys(),
    )
    def test_load_custom_activations(self, tmpdir, activation):
        inp, out = CONTINUOUS_TMAPS[:2], CATEGORICAL_TMAPS[:2]
        params = DEFAULT_PARAMS.copy()
        params['activation'] = activation
        m = make_multimodal_multitask_model(
            inp,
            out,
            **params,
        )
        path = os.path.join(tmpdir, f'm{MODEL_EXT}')
        m.save(path)
        make_multimodal_multitask_model(
            inp,
            out,
            model_file=path,
            **params,
        )

    @pytest.mark.slow
    @pytest.mark.parametrize(
        'input_tmaps',
        MULTIMODAL_UP_TO_4D,
    )
    @pytest.mark.parametrize(
        'output_tmaps',
        MULTIMODAL_UP_TO_4D,
    )
    def test_load_multimodal(self, tmpdir, input_tmaps: List[TensorMap], output_tmaps: List[TensorMap]):
        m = make_multimodal_multitask_model(
            input_tmaps,
            output_tmaps,
            **DEFAULT_PARAMS,
        )
        path = os.path.join(tmpdir, f'm{MODEL_EXT}')
        m.save(path)
        make_multimodal_multitask_model(
            input_tmaps,
            output_tmaps,
            model_file=path,
            **DEFAULT_PARAMS,
        )

    def test_u_connect_auto_encode(self):
        params = DEFAULT_PARAMS.copy()
        params['pool_x'] = params['pool_y'] = 2
        params['conv_layers'] = [8, 8]
        params['dense_blocks'] = [4, 4, 2]
        m = make_multimodal_multitask_model(
            [SEGMENT_IN],
            [SEGMENT_IN],
            u_connect=defaultdict(set, {SEGMENT_IN: {SEGMENT_IN}}),
            **params,
        )
        assert_model_trains([SEGMENT_IN], [SEGMENT_IN], m)

    def test_u_connect_segment(self):
        params = DEFAULT_PARAMS.copy()
        params['pool_x'] = params['pool_y'] = 2
        m = make_multimodal_multitask_model(
            [SEGMENT_IN],
            [SEGMENT_OUT],
            u_connect=defaultdict(set, {SEGMENT_IN: {SEGMENT_OUT}}),
            **params,
        )
        assert_model_trains([SEGMENT_IN], [SEGMENT_OUT], m)

    @pytest.mark.parametrize(
        'input_output_tmaps',
        [
            (CONTINUOUS_TMAPS[:1], [SEGMENT_IN]), ([SEGMENT_IN], CONTINUOUS_TMAPS[:1]),
            ([SEGMENT_IN], [SEGMENT_IN]),
        ],
    )
    def test_multimodal_multitask_variational(self, input_output_tmaps, tmpdir):
        """
        Tests 1d->2d, 2d->1d, (1d,2d)->(1d,2d)
        """
        params = DEFAULT_PARAMS.copy()
        params['bottleneck_type'] = BottleneckType.Variational
        params['pool_x'] = params['pool_y'] = 2
        m = make_multimodal_multitask_model(
            input_output_tmaps[0],
            input_output_tmaps[1],
            **params
        )
        assert_model_trains(input_output_tmaps[0], input_output_tmaps[1], m)
        m.save(os.path.join(tmpdir, 'vae.h5'))
        path = os.path.join(tmpdir, f'm{MODEL_EXT}')
        m.save(path)
        make_multimodal_multitask_model(
            input_output_tmaps[0],
            input_output_tmaps[1],
            model_file=path,
            **DEFAULT_PARAMS,
        )

    def test_u_connect_adaptive_normalization(self):
        params = DEFAULT_PARAMS.copy()
        params['pool_x'] = params['pool_y'] = 2
        params['bottleneck_type'] = BottleneckType.GlobalAveragePoolStructured
        m = make_multimodal_multitask_model(
            [SEGMENT_IN, TMAPS_UP_TO_4D[0]],
            [SEGMENT_OUT],
            u_connect=defaultdict(set, {SEGMENT_IN: {SEGMENT_OUT}}),
            **params,
        )
        assert_model_trains([SEGMENT_IN, TMAPS_UP_TO_4D[0]], [SEGMENT_OUT], m)

    def test_u_connect_no_bottleneck(self):
        params = DEFAULT_PARAMS.copy()
        params['pool_x'] = params['pool_y'] = 2
        params['bottleneck_type'] = BottleneckType.NoBottleNeck
        m = make_multimodal_multitask_model(
            [SEGMENT_IN, TMAPS_UP_TO_4D[0]],
            [SEGMENT_OUT],
            u_connect=defaultdict(set, {SEGMENT_IN: {SEGMENT_OUT}}),
            **params,
        )
        assert_model_trains([SEGMENT_IN, TMAPS_UP_TO_4D[0]], [SEGMENT_OUT], m)

    def test_no_dense_layers(self):
        params = DEFAULT_PARAMS.copy()
        params['dense_layers'] = []
        inp, out = CONTINUOUS_TMAPS[:2], CATEGORICAL_TMAPS[:2]
        m = make_multimodal_multitask_model(
            inp,
            out,
            **DEFAULT_PARAMS,
        )
        assert_model_trains(inp, out, m)

    @pytest.mark.parametrize(
        'output_tmaps',
        [_rotate(PARENT_TMAPS, i) for i in range(len(PARENT_TMAPS))],
    )
    def test_parents(self, output_tmaps):
        assert_model_trains([TMAPS_UP_TO_4D[-1]], output_tmaps)

    @pytest.mark.parametrize(
        'input_output_tmaps',
        [
            (LANGUAGE_TMAP_1HOT_WINDOW, LANGUAGE_TMAP_1HOT_SOFTMAX),
        ],
    )
    def test_language_models(self, input_output_tmaps, tmpdir):
        params = DEFAULT_PARAMS.copy()
        m = make_multimodal_multitask_model(
            tensor_maps_in=input_output_tmaps[0],
            tensor_maps_out=input_output_tmaps[1],
            **params
        )
        assert_model_trains(input_output_tmaps[0], input_output_tmaps[1], m)
        m.save(os.path.join(tmpdir, 'lstm.h5'))
        path = os.path.join(tmpdir, f'm{MODEL_EXT}')
        m.save(path)
        make_multimodal_multitask_model(
            input_output_tmaps[0],
            input_output_tmaps[1],
            model_file=path,
            **DEFAULT_PARAMS,
        )

    @pytest.mark.parametrize(
        'pairs',
        [
            [(CONTINUOUS_TMAPS[2], CONTINUOUS_TMAPS[1])],
            [(CATEGORICAL_TMAPS[2], CATEGORICAL_TMAPS[1])],
            [(CONTINUOUS_TMAPS[2], CONTINUOUS_TMAPS[1]), (CONTINUOUS_TMAPS[2], CATEGORICAL_TMAPS[3])]
        ],
    )
    def test_paired_models(self, pairs, tmpdir):
        params = DEFAULT_PARAMS.copy()
        pair_list = list(set([p[0] for p in pairs] + [p[1] for p in pairs]))
        params['u_connect'] = {tm: [] for tm in pair_list}
        m, encoders, decoders = make_paired_autoencoder_model(
            pairs=pairs,
            tensor_maps_in=pair_list,
            tensor_maps_out=pair_list,
            **params
        )
        assert_model_trains(pair_list, pair_list, m, skip_shape_check=True)
        m.save(os.path.join(tmpdir, 'paired_ae.h5'))
        path = os.path.join(tmpdir, f'm{MODEL_EXT}')
        m.save(path)
        make_paired_autoencoder_model(
            pairs=pairs,
            tensor_maps_in=pair_list,
            tensor_maps_out=pair_list,
            **params
        )

    @pytest.mark.parametrize(
        'pairs',
        [
            [(CONTINUOUS_TMAPS[2], CONTINUOUS_TMAPS[1])],
            [(CATEGORICAL_TMAPS[2], CATEGORICAL_TMAPS[1])],
            [(CONTINUOUS_TMAPS[2], CONTINUOUS_TMAPS[1]), (CONTINUOUS_TMAPS[2], CATEGORICAL_TMAPS[3])]
        ],
    )
    @pytest.mark.parametrize(
        'output_tmaps',
        [
            [CONTINUOUS_TMAPS[0]],
            [CATEGORICAL_TMAPS[0]],
            [CONTINUOUS_TMAPS[0], CATEGORICAL_TMAPS[0]],
        ],
    )
    def test_semi_supervised_paired_models(self, pairs, output_tmaps, tmpdir):
        params = DEFAULT_PARAMS.copy()
        pair_list = list(set([p[0] for p in pairs] + [p[1] for p in pairs]))
        params['u_connect'] = {tm: [] for tm in pair_list}
        m, encoders, decoders = make_paired_autoencoder_model(
            pairs=pairs,
            tensor_maps_in=pair_list,
            tensor_maps_out=pair_list+output_tmaps,
            **params
        )
        assert_model_trains(pair_list, pair_list+output_tmaps, m, skip_shape_check=True)
        m.save(os.path.join(tmpdir, 'paired_ae.h5'))
        path = os.path.join(tmpdir, f'm{MODEL_EXT}')
        m.save(path)
        make_paired_autoencoder_model(
            pairs=pairs,
            tensor_maps_in=pair_list,
            tensor_maps_out=pair_list+output_tmaps,
            **params
        )

@pytest.mark.parametrize(
    'tmaps',
    [_rotate(PARENT_TMAPS, i) for i in range(len(PARENT_TMAPS))],
)
def test_parent_sort(tmaps):
    assert parent_sort(tmaps) == PARENT_TMAPS


@pytest.mark.parametrize(
    'tmaps',
    [_rotate(CYCLE_PARENTS, i) for i in range(len(CYCLE_PARENTS))],
)
def test_parent_sort_cycle(tmaps):
    with pytest.raises(ValueError):
        parent_sort(tmaps)


@pytest.mark.parametrize(
    'tmaps',
    [_rotate(PARENT_TMAPS + TMAPS_UP_TO_4D, i) for i in range(len(PARENT_TMAPS))],
)
def test_parent_sort_idempotent(tmaps):
    assert parent_sort(tmaps) == parent_sort(parent_sort(tmaps)) == parent_sort(parent_sort(parent_sort(tmaps)))


@pytest.mark.parametrize(
    'tmap_out',
    TMAPS_UP_TO_4D,
)
@pytest.mark.parametrize(
    'u_connect_out',
    TMAPS_UP_TO_4D,
)
def test_check_no_bottleneck(tmap_out, u_connect_out):
    u_connect = defaultdict(set, {tmap_out: {u_connect_out}})
    assert check_no_bottleneck(u_connect, [tmap_out]) == (u_connect_out == tmap_out)


class TestModelPerformance:
    @pytest.mark.slow
    def test_brain_seg(self, tmpdir):
        tensor_path = '/mnt/disks/brains-all-together/2020-02-11/'
        if not os.path.exists(tensor_path):
            pytest.skip('To test brain segmentation performance, attach disk brains-all-together')

        from ml4h.tensor_from_file import TMAPS
        from ml4h.tensor_generators import test_train_valid_tensor_generators, big_batch_from_minibatch_generator
        from multiprocessing import cpu_count
        from sklearn.metrics import average_precision_score

        tmaps_in = [TMAPS['t1_30_slices_4d']]
        tmaps_out = [TMAPS['t1_seg_30_slices']]
        m = make_multimodal_multitask_model(
            tensor_maps_in=tmaps_in, tensor_maps_out=tmaps_out,
            activation='relu',
            learning_rate=1e-3,
            bottleneck_type=BottleneckType.GlobalAveragePoolStructured,
            optimizer='radam',
            dense_layers=[16, 64],
            conv_layers=[32],
            dense_blocks=[32, 24, 16],
            block_size=3,
            conv_type='conv',
            conv_x=[3], conv_y=[3], conv_z=[2],
            pool_x=2, pool_y=2, pool_z=1,
            pool_type='max',
            u_connect=defaultdict(set, {tmaps_in[0]: {tmaps_out[0]}}),
        )
        batch_size = 2
        generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(
            tmaps_in, tmaps_out,
            tensors=tensor_path,
            batch_size=batch_size,
            valid_ratio=.2,
            test_ratio=.2,
            num_workers=cpu_count(),
            cache_size=1e9 / cpu_count(),
            balance_csvs=[],
            training_steps=64,
            validation_steps=18,
            test_modulo=0,
        )
        try:
            m = train_model_from_generators(
                model=m,
                generate_train=generate_train, generate_valid=generate_valid,
                training_steps=64, validation_steps=18, epochs=24, patience=22, batch_size=batch_size,
                output_folder=str(tmpdir), run_id='brain_seg_test',
                inspect_model=True, inspect_show_labels=True,
            )
            test_data, test_labels, test_paths = big_batch_from_minibatch_generator(
                generate_test, 12,
            )
        finally:
            generate_train.kill_workers()
            generate_test.kill_workers()
            generate_valid.kill_workers()
        y_prediction = m.predict(test_data, batch_size=batch_size)
        y_truth = np.array(test_labels[tmaps_out[0].output_name()])
        expected_precisions = {
            'not_brain_tissue': 1.,
            'csf': .921,
            'grey': .963,
            'white': .989,
        }
        actual_precisions = {}
        for name, idx in tmaps_out[0].channel_map.items():
            average_precision = average_precision_score(
                y_truth[..., idx].flatten(), y_prediction[..., idx].flatten(),
            )
            actual_precisions[name] = average_precision
        for name in expected_precisions:
            assert actual_precisions[name] >= expected_precisions[name] - MEAN_PRECISION_EPS
