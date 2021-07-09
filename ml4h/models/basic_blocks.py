from typing import Dict, List

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Dense, LSTM, Concatenate

from ml4h.models.Block import Block
from ml4h.TensorMap import TensorMap
from ml4h.models.layer_wrappers import _activation_layer, _regularization_layer, _normalization_layer

Tensor = tf.Tensor


class DenseBlock(Block):
    def __init__(
            self,
            *,
            widths: List[int],
            activation: str,
            normalization: str,
            regularization: str,
            regularization_rate: float,
            name: str = None,
            **kwargs
    ):
        final_dense = Dense(units=widths[-1], name=name) if name else Dense(units=widths[-1])
        self.denses = [Dense(units=width) for width in widths[:-1]] + [final_dense]
        self.activations = [_activation_layer(activation) for _ in widths]
        self.regularizations = [_regularization_layer(1, regularization, regularization_rate) for _ in widths]
        self.norms = [_normalization_layer(normalization) for _ in widths]

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]]) -> Tensor:
        for dense, normalize, activate, regularize in zip(self.denses, self.norms, self.activations, self.regularizations):
            x = normalize(regularize(activate(dense(x))))
        return x


class DenseEncoder(Block):
    def __init__(
            self,
            *,
            tensor_map: TensorMap,
            dense_layers: List[int],
            activation: str,
            dense_normalize: str,
            dense_regularize: str,
            dense_regularize_rate: float,
            **kwargs,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        self.fully_connected = DenseBlock(
            widths=dense_layers,
            activation=activation,
            normalization=dense_normalize,
            regularization=dense_regularize,
            regularization_rate=dense_regularize_rate,
            name=tensor_map.embed_name(),
        )

    def can_apply(self):
        return self.tensor_map.axes() == 1 and not self.tensor_map.is_embedding()

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]]) -> Tensor:
        if not self.can_apply():
            return x
        y = self.fully_connected(x, intermediates)
        intermediates[self.tensor_map].append(y)
        return y


class DenseDecoder(Block):
    def __init__(
            self,
            tensor_map: TensorMap,
            activation: str,
            parents: List[TensorMap] = None,
            **kwargs,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        self.parents = parents
        self.activation = _activation_layer(activation)
        self.dense = Dense(units=tensor_map.shape[0], name=tensor_map.output_name(), activation=tensor_map.activation)
        self.units = tensor_map.annotation_units

    def can_apply(self):
        return self.tensor_map.axes() == 1

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]]) -> Tensor:
        if not self.can_apply():
            return x
        if self.parents:
            x = Concatenate()([x] + [intermediates[parent][-1] for parent in self.parents])
            x = Dense(units=self.units)(x)
            x = self.dense(self.activation(x))
        x = self.dense(x)
        intermediates[self.tensor_map].append(x)
        return x


class ModelAsBlock(Block):
    """Takes a serialized model and applies it, can be used to encode or decode Tensors"""
    def __init__(
            self,
            *,
            tensor_map: TensorMap,
            model: Model,
            **kwargs,
    ):
        self.tensor_map = tensor_map
        self.model = model

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]]) -> Tensor:
        x = self.model(x)
        intermediates[self.tensor_map].append(x)
        return x


class LSTMEncoderBlock(Block):
    def __init__(
            self,
            tensor_map,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        self.lstm = LSTM(tensor_map.annotation_units)

    def can_apply(self):
        return self.tensor_map.is_language()

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]]) -> Tensor:
        if not self.can_apply():
            return x
        return self.lstm(x)


class LanguageDecoderBlock(Block):
    def __init__(
            self,
            tensor_map: TensorMap,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        self.dense = Dense(tensor_map.shape[-1], activation=tensor_map.activation, name=tensor_map.output_name())

    def can_apply(self):
        return self.tensor_map.is_language()

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]]) -> Tensor:
        if not self.can_apply():
            return x
        return self.dense(x)
