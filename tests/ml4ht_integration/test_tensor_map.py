from itertools import cycle
from datetime import datetime, timedelta

import pytest
import numpy as np
from torch.utils.data import DataLoader

from ml4h.test_utils import TMAPS_UP_TO_4D, build_hdf5s
from ml4h.ml4ht_integration.tensor_map import TensorMapSampleGetter, tensor_map_from_data_description
from ml4h.models.legacy_models import make_multimodal_multitask_model, BottleneckType
from ml4h.TensorMap import Interpretation
from ml4ht.data.data_loader import numpy_collate_fn, SampleGetterIterableDataset
from ml4ht.data.data_description import DataDescription
from ml4ht.data.util.date_selector import (
    DateRangeOptionPicker,
    first_dt,
    DATE_OPTION_KEY,
)
from ml4ht.data.sample_getter import DataDescriptionSampleGetter


# Tests
def test_tensor_map_from_data_description():
    tmap_in = tensor_map_from_data_description(
        DD1, Interpretation.CONTINUOUS, (1,),
    )
    tmap_out = tensor_map_from_data_description(
        DD2, Interpretation.CONTINUOUS, (1,),
        loss='logcosh',
        metrics=['mae', 'mse'],
    )
    model = make_multimodal_multitask_model(
        [tmap_in], [tmap_out],
    )
    data_set = SampleGetterIterableDataset(
        sample_ids=list(RAW_DATA),
        sample_getter=SAMPLE_GETTER,
    )
    data_loader = DataLoader(
        data_set, batch_size=4, num_workers=2,
        collate_fn=numpy_collate_fn,
    )
    # model can train?
    history = model.fit(iter(data_loader)).history
    # metrics recorded?
    assert 'mae' in history
    assert 'mse' in history


class TestTensorMapSampleGetter:
    def test_expected_batches(self, expected_data):
        paths = [k[0] for k in expected_data]
        sample_getter = TensorMapSampleGetter(
            TMAPS_UP_TO_4D,
            TMAPS_UP_TO_4D,
        )
        for path in paths:
            in_batch, out_batch = sample_getter(path)
            for tmap in TMAPS_UP_TO_4D:
                np.testing.assert_allclose(
                    in_batch[tmap.input_name()],
                    expected_data[path, tmap]
                )
                np.testing.assert_allclose(
                    out_batch[tmap.output_name()],
                    expected_data[path, tmap]
                )

    def test_train_model(self, expected_data, model):
        paths = [k[0] for k in expected_data]
        sample_getter = TensorMapSampleGetter(
            TMAPS_UP_TO_4D,
            TMAPS_UP_TO_4D,
        )
        dataset = SampleGetterIterableDataset(
            sample_getter=sample_getter,
            sample_ids=paths,
        )
        data_loader = DataLoader(
            dataset, batch_size=2, num_workers=3,
            collate_fn=numpy_collate_fn,
        )
        val_loader = DataLoader(
            dataset, batch_size=2, num_workers=3,
            collate_fn=numpy_collate_fn,
        )
        history = model.fit(
            cycle(data_loader), validation_data=cycle(val_loader),
            steps_per_epoch=10, validation_steps=2, epochs=3,
        ).history
        for tmap in TMAPS_UP_TO_4D:
            for metric in tmap.metrics:
                metric_name = metric if type(metric) == str else metric.__name__
                name = f'{tmap.output_name()}_{metric_name}'
                assert name in history
                assert f'val_{name}' in history


# Setup
@pytest.fixture(scope="module")
def expected_data(tmpdir_factory):
    temp_dir = tmpdir_factory.mktemp('data')
    return build_hdf5s(temp_dir, TMAPS_UP_TO_4D, n=pytest.N_TENSORS)


@pytest.fixture(scope="function")
def model():
    return make_multimodal_multitask_model(
        tensor_maps_in=TMAPS_UP_TO_4D,
        tensor_maps_out=TMAPS_UP_TO_4D,
        conv_x=[3], conv_y=[3], conv_z=[1],
        pool_x=1, pool_y=1, pool_z=1,
        dense_blocks=[4], dense_layers=[4],
        block_size=1,
        activation='relu', learning_rate=1e-5,
        bottleneck_type=BottleneckType.FlattenRestructure,
        optimizer='sgd', conv_layers=[8], conv_type='conv',
    )


RAW_DATA = {
    0: {
        datetime(year=2000, month=3, day=1): np.array([10]),
    },
    1: {
        datetime(year=2000, month=3, day=1): np.array([-1]),
    },
    2: {
        datetime(year=2000, month=3, day=1): np.array([10]),
    },
    3: {
        datetime(year=2000, month=3, day=1): np.array([10]),
    },
}


class DictionaryDataDescription(DataDescription):
    def __init__(self, data, fail_idx, name):
        self.data = data
        self.fail_idx = fail_idx
        self._name = name

    def get_loading_options(self, sample_id):
        return [{DATE_OPTION_KEY: dt} for dt in self.data[sample_id]]

    @property
    def name(self) -> str:
        return self._name

    def get_raw_data(self, sample_id, loading_option):
        if sample_id == self.fail_idx:
            raise ValueError("Bad idx")
        dt = loading_option[DATE_OPTION_KEY]
        return self.data[sample_id][dt]


DD1 = DictionaryDataDescription(RAW_DATA, 1, 'dd1')
DD2 = DictionaryDataDescription(RAW_DATA, 2, 'dd2')
DATE_OPTION_PICKER = DateRangeOptionPicker(
    reference_data_description=DD1,
    reference_date_chooser=first_dt,
    time_before=timedelta(days=0),
    time_after=timedelta(days=5),
)
SAMPLE_GETTER = DataDescriptionSampleGetter(
    input_data_descriptions=[DD1],
    output_data_descriptions=[DD2],
    option_picker=DATE_OPTION_PICKER,
)


