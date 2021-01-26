import numpy as np
from abc import ABC, abstractmethod

from ml4h.defines import EPS
from tensorflow.python.keras.applications import imagenet_utils


class Normalizer(ABC):
    @abstractmethod
    def normalize(self, tensor: np.ndarray) -> np.ndarray:
        """Shape preserving transformation"""
        pass

    def un_normalize(self, tensor: np.ndarray) -> np.ndarray:
        """The inverse of normalize if possible. Otherwise identity."""
        return tensor


class Standardize(Normalizer):
    def __init__(self, mean: float, std: float):
        self.mean, self.std = mean, std

    def normalize(self, tensor: np.ndarray) -> np.ndarray:
        return (tensor - self.mean) / self.std

    def un_normalize(self, tensor: np.ndarray) -> np.ndarray:
        return tensor * self.std + self.mean


class ZeroMeanStd1(Normalizer):
    def normalize(self, tensor: np.ndarray) -> np.ndarray:
        tensor -= np.mean(tensor)
        tensor /= np.std(tensor) + EPS
        return tensor


class NonZeroNormalize(Normalizer):
    def normalize(self, tensor: np.ndarray) -> np.ndarray:
        nonzero = tensor > 0
        tensor[nonzero] = (tensor[nonzero] - tensor[nonzero].mean() + 1e-9) / (
            tensor[nonzero].std() + 1e-9
        )
        return tensor


class TopKNormalize(Normalizer):
    def __init__(self, n_top: int = 50):
        self.n_top = n_top

    def normalize(self, tensor: np.ndarray) -> np.ndarray:
        """Find top K itensity voxels are set upper range to the mean of those"""
        upper = np.mean(sorted(np.max(tensor, axis=-1).flatten())[::-1][0:self.n_top])
        tensor = np.where(tensor >= upper, upper, tensor)
        tensor /= tensor.max()
        return tensor


class ImagenetNormalizeTorch(Normalizer):
    def normalize(self, tensor: np.ndarray) -> np.ndarray:
        # This is equivalent to:
        # x /= 255.
        # mean = [0.485, 0.456, 0.406]
        # std = [0.229, 0.224, 0.225]
        # when mode is torch
        return imagenet_utils.preprocess_input(tensor, data_format=None, mode="torch")
