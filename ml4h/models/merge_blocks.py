from typing import Dict, List, Tuple

import tensorflow as tf
import tensorflow.keras.backend as K
import tensorflow_probability as tfp
from tensorflow.keras.layers import concatenate, Flatten, Average, Layer

from ml4h.models.Block import Block
from ml4h.TensorMap import TensorMap
from ml4h.models.basic_blocks import DenseBlock
from ml4h.models.layer_wrappers import global_average_pool

Tensor = tf.Tensor
tfd = tfp.distributions


class FlatConcatBlock(Block):
    """
    Flattens then concatenates all inputs
    """
    def __init__(self,  **kwargs):
        pass

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        y = [Flatten()(x[-1]) for tm, x in intermediates.items() if not tm.is_embedding()]
        y = concatenate(y) if len(y) > 1 else y[0]
        return y


class FlatConcatDenseBlock(Block):
    """
    Flattens then concatenates all inputs, applies a dense layer
    """
    def __init__(
            self,
            activation: str,
            dense_layers: List[int],
            dense_normalize: str,
            dense_regularize: str,
            dense_regularize_rate: float,
            **kwargs,
    ):
        self.fully_connected = DenseBlock(
            widths=dense_layers,
            activation=activation,
            normalization=dense_normalize,
            regularization=dense_regularize,
            regularization_rate=dense_regularize_rate,
            name='embed',
        ) if dense_layers else None

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        y = [Flatten()(x[-1]) for tm, x in intermediates.items() if not tm.is_embedding()]
        y = concatenate(y) if len(y) > 1 else y[0]
        y = self.fully_connected(y, intermediates) if self.fully_connected else y
        return y


class GlobalAveragePoolBlock(Block):
    """
    GAPs then concatenates all inputs, applies a dense layer
    """
    def __init__(
            self,
            activation: str,
            dense_layers: List[int],
            dense_normalize: str,
            dense_regularize: str,
            dense_regularize_rate: float,
            **kwargs,
    ):
        self.fully_connected = DenseBlock(
            widths=dense_layers,
            activation=activation,
            normalization=dense_normalize,
            regularization=dense_regularize,
            regularization_rate=dense_regularize_rate,
            name='embed',
        ) if dense_layers else None

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        y = [Flatten()(x[-1]) for tm, x in intermediates.items() if tm.axes() == 1]  # Flat tensors
        y += [global_average_pool(x[-2 if self.fully_connected else -1]) for tm, x in intermediates.items() if tm.axes() > 1]  # Structured tensors
        y = concatenate(y) if len(y) > 1 else y[0]
        y = self.fully_connected(y, intermediates) if self.fully_connected else y
        return y


class AverageBlock(Block):
    """
    Average the last tensors in intermediates dictionary
    """
    def __init__(self, **kwargs):
        pass

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        return Average()([Flatten()(x[-1]) for tm, x in intermediates.items()])


class EncodeIdentityBlock(Block):
    """
    Adds the input tensor to the intermediates dictionary, useful for TensorMaps with pretrained embeddings
    """
    def __init__(self, tensor_map, **kwargs):
        self.tensor_map = tensor_map

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        if self.tensor_map.is_embedding():
            intermediates[self.tensor_map].append(x)
        return x


class PairLossBlock(Block):
    """
    Flattens or GAPs then concatenates all inputs, applies a dense layer, then restructures to provided shapes
    """
    def __init__(
            self,
            pairs: List[Tuple[TensorMap, TensorMap]],
            pair_loss: str = 'cosine',
            pair_loss_weight: float = 1.0,
            **kwargs,
    ):
        self.pairs = pairs
        if pair_loss == 'cosine':
            self.loss_layer = CosineLossLayer(pair_loss_weight)
        elif pair_loss == 'euclid':
            self.loss_layer = L2LossLayer(pair_loss_weight)

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        for left, right in self.pairs:
            x = self.loss_layer([intermediates[left][-1], intermediates[right][-1]])
            intermediates[left].extend(x)
            intermediates[right].extend(x)
        return intermediates[left][-1]


def l2_norm(x, axis=None):
    """
    takes an input tensor and returns the l2 norm along specified axis
    """

    square_sum = K.sum(K.square(x), axis=axis, keepdims=True)
    norm = K.sqrt(K.maximum(square_sum, K.epsilon()))

    return norm


def pairwise_cosine_difference(t1, t2):
    t1_norm = t1 / l2_norm(t1, axis=-1)
    t2_norm = t2 / l2_norm(t2, axis=-1)
    dot = K.clip(K.batch_dot(t1_norm, t2_norm), -1, 1)
    return K.mean(tf.acos(dot))


class CosineLossLayer(Layer):
    """Layer that creates an Cosine loss."""

    def __init__(self, weight, **kwargs):
        super(CosineLossLayer, self).__init__(**kwargs)
        self.weight = weight

    def get_config(self):
        config = super().get_config().copy()
        config.update({'weight': self.weight})
        return config

    def call(self, inputs):
        # We use `add_loss` to create a regularization loss
        # that depends on the inputs.
        self.add_loss(self.weight * pairwise_cosine_difference(inputs[0], inputs[1]))
        return inputs


class L2LossLayer(Layer):
    """Layer that creates an L2 loss."""

    def __init__(self, weight, **kwargs):
        super(L2LossLayer, self).__init__(**kwargs)
        self.weight = weight

    def get_config(self):
        config = super().get_config().copy()
        config.update({'weight': self.weight})
        return config

    def call(self, inputs):
        self.add_loss(self.weight * tf.reduce_sum(tf.square(inputs[0] - inputs[1])))
        return inputs


class VariationalDiagNormal(Layer):
    def __init__(
            self,
            latent_size: int,
            kl_divergence_weight: float = 1.,
            **kwargs
    ):
        self.latent_size = latent_size
        self.kl_divergence_weight = kl_divergence_weight
        super(VariationalDiagNormal, self).__init__(**kwargs)
        self.prior = tfd.MultivariateNormalDiag(loc=tf.zeros([latent_size]), scale_identity_multiplier=1.0)

    def call(self, mu: Tensor, log_sigma: Tensor, **kwargs):
        """mu and sigma must be shape (None, latent_size)"""
        approx_posterior = tfd.MultivariateNormalDiag(loc=mu, scale_diag=tf.math.exp(log_sigma))
        kl = tf.reduce_mean(tfd.kl_divergence(approx_posterior, self.prior))
        self.add_loss(kl * self.kl_divergence_weight)
        self.add_metric(kl, 'mean', name='KL_divergence')
        return approx_posterior.sample()

    def get_config(self):
        return {'latent_size': self.latent_size, 'kl_divergence_weight': self.kl_divergence_weight}


