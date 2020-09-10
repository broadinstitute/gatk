import numpy as np
from abc import ABC, abstractmethod

from ml4h.defines import EPS


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
