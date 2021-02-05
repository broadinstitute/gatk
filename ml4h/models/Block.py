from abc import ABC, abstractmethod
from typing import Dict, List, Callable

import tensorflow as tf

from ml4h.TensorMap import TensorMap

Tensor = tf.Tensor
Block = Callable[[Tensor, Dict[TensorMap, List[Tensor]]], Tensor]


class Block(ABC):
    @abstractmethod
    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        pass

    def can_apply(self):
        return True
