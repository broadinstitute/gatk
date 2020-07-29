import os
import h5py
import numpy as np
from typing import List, Tuple, Dict
from itertools import product

from ml4cvd.defines import TENSOR_EXT
from ml4cvd.TensorMap import TensorMap, Interpretation


CONTINUOUS_TMAPS = [
    TensorMap(f'{n}d_cont', shape=tuple(range(2, n + 2)), interpretation=Interpretation.CONTINUOUS)
    for n in range(1, 6)
]
CATEGORICAL_TMAPS = [
    TensorMap(
        f'{n}d_cat', shape=tuple(range(2, n + 2)),
        interpretation=Interpretation.CATEGORICAL,
        channel_map={f'c_{i}': i for i in range(n + 1)},
    )
    for n in range(1, 6)
]
LANGUAGE_TMAP_1HOT_WINDOW = [
    TensorMap(
        f'language_1hot_window', shape=(32, 26),
        interpretation=Interpretation.LANGUAGE,
        channel_map={f'c_{i}': i for i in range(26)},
    ),
]
LANGUAGE_TMAP_1HOT_SOFTMAX = [
    TensorMap(
        f'language_1hot_out', shape=(26,),
        interpretation=Interpretation.LANGUAGE,
        channel_map={f'c_{i}': i for i in range(26)},
    ),
]

TMAPS_UP_TO_4D = CONTINUOUS_TMAPS[:-1] + CATEGORICAL_TMAPS[:-1]
TMAPS_5D = CONTINUOUS_TMAPS[-1:] + CATEGORICAL_TMAPS[-1:]
MULTIMODAL_UP_TO_4D = [list(x) for x in product(CONTINUOUS_TMAPS[:-1], CATEGORICAL_TMAPS[:-1])]
SEGMENT_IN = TensorMap(f'2d_for_segment_in', shape=(32, 32, 1), interpretation=Interpretation.CONTINUOUS, metrics=['mse'])
SEGMENT_OUT = TensorMap(f'2d_for_segment_out', shape=(32, 32, 2), interpretation=Interpretation.CATEGORICAL, channel_map={'yes': 0, 'no': 1})


TMAPS = {
    tmap.name: tmap
    for tmap in CONTINUOUS_TMAPS + CATEGORICAL_TMAPS
}
PARENT_TMAPS = [
    TensorMap(f'parent_test_{i}', shape=(1,), interpretation=Interpretation.CONTINUOUS)
    for i in range(3)
]
for i in range(len(PARENT_TMAPS)):
    PARENT_TMAPS[i].parents = PARENT_TMAPS[:i]
CYCLE_PARENTS = [
    TensorMap(f'parent_test_cycle_{i}', shape=(1,), interpretation=Interpretation.CONTINUOUS)
    for i in range(3)
]
for i in range(len(CYCLE_PARENTS)):
    CYCLE_PARENTS[i].parents = [CYCLE_PARENTS[i - 1]]  # 0th tmap will be child of last


def build_hdf5s(path: str, tensor_maps: List[TensorMap], n=5) -> Dict[Tuple[str, TensorMap], np.ndarray]:
    """
    Builds hdf5s at path given TensorMaps. Only works for Continuous and Categorical TensorMaps.
    """
    out = {}
    for i in range(n):
        hd5_path = os.path.join(path, f'{i}{TENSOR_EXT}')
        with h5py.File(hd5_path, 'w') as hd5:
            for tm in tensor_maps:
                if tm.is_continuous():
                    value = np.full(tm.shape, fill_value=i, dtype=np.float32)
                elif tm.is_categorical():
                    value = np.zeros(tm.shape, dtype=np.float32)
                    value[..., i % tm.shape[-1]] = 1
                else:
                    raise NotImplementedError(f'Cannot automatically build hdf5 from interpretation "{tm.interpretation}"')
                hd5.create_dataset(tm.hd5_key_guess(), data=value)
                out[(hd5_path, tm)] = value
    return out
