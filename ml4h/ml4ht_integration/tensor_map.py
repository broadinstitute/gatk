from typing import Callable, List, Tuple

import h5py

from ml4h.TensorMap import TensorMap, Interpretation
from ml4ht.data.data_description import DataDescription
from ml4ht.data.defines import SampleID, LoadingOption, Tensor, Batch


class TensorMapSampleGetter:
    def __init__(
            self,
            tensor_maps_in: List[TensorMap],
            tensor_maps_out: List[TensorMap],
            augment: bool = False,
    ):
        self.tensor_maps_in = tensor_maps_in
        self.tensor_maps_out = tensor_maps_out
        self.augment = augment

    def __call__(self, path: str) -> Batch:
        dependents = {}
        with h5py.File(path, 'r') as hd5:
            in_batch = {}
            for tm in self.tensor_maps_in:
                in_batch[tm.input_name()] = tm.postprocess_tensor(
                    tm.tensor_from_file(tm, hd5, dependents),
                    augment=self.augment, hd5=hd5,
                )
            out_batch = {}
            for tm in self.tensor_maps_out:
                out_batch[tm.output_name()] = tm.postprocess_tensor(
                    tm.tensor_from_file(tm, hd5, dependents),
                    augment=self.augment, hd5=hd5,
                )
        return in_batch, out_batch


def _not_implemented_tensor_from_file(_, __, ___=None):
    """Used to make sure TensorMap is never used to load data"""
    raise NotImplementedError('This TensorMap cannot load data.')


def tensor_map_from_data_description(
        data_description: DataDescription,
        interpretation: Interpretation,
        shape: Tuple[int, ...],
        **tensor_map_kwargs,
) -> TensorMap:
    """
    Allows a DataDescription to be used in the model factory
    """
    tmap = TensorMap(
        name=data_description.name,
        interpretation=interpretation,
        shape=shape,
        tensor_from_file=_not_implemented_tensor_from_file,
        **tensor_map_kwargs,
    )
    # hacky way to make sure the model factory will work with data description data loaders
    tmap.input_name = lambda: data_description.name
    tmap.output_name = lambda: data_description.name
    return tmap
