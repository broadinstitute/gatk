import logging
from typing import List

from torch.utils.data import DataLoader
from ml4ht.data.data_loader import SampleGetterIterableDataset, numpy_collate_fn

from ml4h.TensorMap import TensorMap
from ml4h.defines import TensorGeneratorABC
from ml4h.ml4ht_integration.tensor_map import TensorMapSampleGetter


class TensorMapDataLoader(TensorGeneratorABC):
    def __init__(
        self, batch_size: int, input_maps: List[TensorMap], output_maps: List[TensorMap],
        paths: List[str], num_workers: int,
        keep_paths: bool = False,
        drop_last: bool = True,
        augment: bool = False,
        **kwargs,
    ):
        self.paths = paths
        self.input_maps = input_maps
        self.output_maps = output_maps
        self.keep_paths = keep_paths
        self.sample_getter = TensorMapSampleGetter(
            input_maps, output_maps, augment,
            return_path=keep_paths,
        )
        self.dset = SampleGetterIterableDataset(
            paths, self.sample_getter,
            get_epoch=SampleGetterIterableDataset.shuffle_get_epoch,
        )
        self.data_loader = DataLoader(
            self.dset, batch_size=batch_size, num_workers=num_workers,
            collate_fn=self._collate_fn, drop_last=drop_last,
        )
        self.iter_loader = iter(self.data_loader)


    def _collate_fn(self, batches):
        if self.keep_paths:
            return numpy_collate_fn([batch[:2] for batch in batches]) + ([batch[2] for batch in batches],)
        return numpy_collate_fn(batches)

    @staticmethod
    def can_apply(paths, weights, mixup, siamese, **kwargs):
        """Can you substitute this TensorGenerator for the ml4h legacy TensorGenerator"""
        if isinstance(paths[0], list):
            raise NotImplementedError(
                "TensorMapDataLoader cannot sample from multiple lists of paths. Pass a list of paths for the paths argument"
            )
        if weights is not None:
            raise NotImplementedError(
                "TensorMapDataLoader cannot sample from multiple lists of paths. Do not pass 'weights' argument."
            )
        if mixup:
            raise NotImplementedError("Mixup not implemented for TensorMapDataLoader")
        if siamese:
            raise NotImplementedError("Siamese not implemented for TensorMapDataLoader")

    def __iter__(self):
        return self

    def __next__(self):
        """Infinite iterator over data loader"""
        try:
            return next(self.iter_loader)
        except StopIteration:
            self.iter_loader = iter(self.data_loader)
            logging.info("Completed one epoch.")
            return next(self.iter_loader)

    def __call__(self):
        try:
            next(self.iter_loader)
        except StopIteration:
            self.iter_loader = iter(self.data_loader)
        return self

    def kill_workers(self):
        """necessary for legacy compatibility"""
        pass
