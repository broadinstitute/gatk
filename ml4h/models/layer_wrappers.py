
# Imports
import os
import time
import logging
import numpy as np
from enum import Enum, auto
from itertools import chain
from abc import ABC, abstractmethod
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Iterable, Union, Optional, Set, Sequence, Callable, DefaultDict, Any

import tensorflow as tf
import tensorflow_addons as tfa
import tensorflow.keras.backend as K
from tensorflow.keras.callbacks import History
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.utils import model_to_dot
from tensorflow.keras.layers import LeakyReLU, PReLU, ELU, ThresholdedReLU, Lambda, Reshape, LayerNormalization
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau, Callback
from tensorflow.keras.layers import SpatialDropout1D, SpatialDropout2D, SpatialDropout3D, add, concatenate
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization, Activation, Flatten, LSTM, RepeatVector
from tensorflow.keras.layers import Conv1D, Conv2D, Conv3D, UpSampling1D, UpSampling2D, UpSampling3D, MaxPooling1D
from tensorflow.keras.layers import MaxPooling2D, MaxPooling3D, Average, AveragePooling1D, AveragePooling2D, AveragePooling3D, Layer
from tensorflow.keras.layers import SeparableConv1D, SeparableConv2D, DepthwiseConv2D, Concatenate, Add
from tensorflow.keras.layers import GlobalAveragePooling1D, GlobalAveragePooling2D, GlobalAveragePooling3D


Tensor = tf.Tensor

ACTIVATION_CLASSES = {
    'leaky': LeakyReLU(),
    'prelu': PReLU(),
    'elu': ELU(),
    'thresh_relu': ThresholdedReLU,
}
ACTIVATION_FUNCTIONS = {
    'swish': tf.nn.swish,
    'gelu': tfa.activations.gelu,
    'lisht': tfa.activations.lisht,
    'mish': tfa.activations.mish,
}
NORMALIZATION_CLASSES = {
    'batch_norm': BatchNormalization,
    'layer_norm': LayerNormalization,
    'instance_norm': tfa.layers.InstanceNormalization,
    'poincare_norm': tfa.layers.PoincareNormalize,
}

CONV_REGULARIZATION_CLASSES = {
    # class name -> (dimension -> class)
    'spatial_dropout': {2: SpatialDropout1D, 3: SpatialDropout2D, 4: SpatialDropout3D},
    'dropout': defaultdict(lambda _: Dropout),
}
DENSE_REGULARIZATION_CLASSES = {
    'dropout': Dropout,  # TODO: add l1, l2
}


def _one_by_n_kernel(dimension):
    return tuple([1] * (dimension - 1))


def _conv_layer_from_kind_and_dimension(
        dimension: int, conv_layer_type: str, conv_x: List[int], conv_y: List[int], conv_z: List[int],
) -> Tuple[Layer, List[Tuple[int, ...]]]:
    if dimension == 4 and conv_layer_type == 'conv':
        conv_layer = Conv3D
        kernel = zip(conv_x, conv_y, conv_z)
    elif dimension == 3 and conv_layer_type == 'conv':
        conv_layer = Conv2D
        kernel = zip(conv_x, conv_y)
    elif dimension == 2 and conv_layer_type == 'conv':
        conv_layer = Conv1D
        kernel = zip(conv_x)
    elif dimension == 3 and conv_layer_type == 'separable':
        conv_layer = SeparableConv2D
        kernel = zip(conv_x, conv_y)
    elif dimension == 2 and conv_layer_type == 'separable':
        conv_layer = SeparableConv1D
        kernel = zip(conv_x)
    elif dimension == 3 and conv_layer_type == 'depth':
        conv_layer = DepthwiseConv2D
        kernel = zip(conv_x, conv_y)
    else:
        raise ValueError(f'Unknown convolution type: {conv_layer_type} for dimension: {dimension}')
    return conv_layer, list(kernel)


def _pool_layers_from_kind_and_dimension(dimension, pool_type, pool_number, pool_x, pool_y, pool_z):
    if dimension == 4 and pool_type == 'max':
        return [MaxPooling3D(pool_size=(pool_x, pool_y, pool_z)) for _ in range(pool_number)]
    elif dimension == 3 and pool_type == 'max':
        return [MaxPooling2D(pool_size=(pool_x, pool_y)) for _ in range(pool_number)]
    elif dimension == 2 and pool_type == 'max':
        return [MaxPooling1D(pool_size=pool_x) for _ in range(pool_number)]
    elif dimension == 4 and pool_type == 'average':
        return [AveragePooling3D(pool_size=(pool_x, pool_y, pool_z)) for _ in range(pool_number)]
    elif dimension == 3 and pool_type == 'average':
        return [AveragePooling2D(pool_size=(pool_x, pool_y)) for _ in range(pool_number)]
    elif dimension == 2 and pool_type == 'average':
        return [AveragePooling1D(pool_size=pool_x) for _ in range(pool_number)]
    else:
        raise ValueError(f'Unknown pooling type: {pool_type} for dimension: {dimension}')


def _upsampler(dimension, pool_x, pool_y, pool_z):
    if dimension == 4:
        return UpSampling3D(size=(pool_x, pool_y, pool_z))
    elif dimension == 3:
        return UpSampling2D(size=(pool_x, pool_y))
    elif dimension == 2:
        return UpSampling1D(size=pool_x)


def _activation_layer(activation: str) -> Activation:
    return (
        ACTIVATION_CLASSES.get(activation, None)
        or Activation(ACTIVATION_FUNCTIONS.get(activation, None) or activation)
    )


def _normalization_layer(norm: str) -> Layer:
    if not norm:
        return lambda x: x
    return NORMALIZATION_CLASSES[norm]()


def _regularization_layer(dimension: int, regularization_type: str, rate: float):
    if not regularization_type:
        return lambda x: x
    if regularization_type in DENSE_REGULARIZATION_CLASSES:
        return DENSE_REGULARIZATION_CLASSES[regularization_type](rate)
    return CONV_REGULARIZATION_CLASSES[regularization_type][dimension](rate)


def adaptive_normalize_from_tensor(tensor: Tensor, target: Tensor) -> Tensor:
    """Uses Dense layers to convert `tensor` to a mean and standard deviation to normalize `target`"""
    return adaptive_normalization(mu=Dense(target.shape[-1])(tensor), sigma=Dense(target.shape[-1])(tensor), target=target)


def adaptive_normalization(mu: Tensor, sigma: Tensor, target: Tensor) -> Tensor:
    target = tfa.layers.InstanceNormalization()(target)
    normalizer_shape = (1,) * (len(target.shape) - 2) + (target.shape[-1],)
    mu = Reshape(normalizer_shape)(mu)
    sigma = Reshape(normalizer_shape)(sigma)
    target *= sigma
    target += mu
    return target


def global_average_pool(x: Tensor) -> Tensor:
    if len(x.shape) == 3:
        return GlobalAveragePooling1D()(x)
    elif len(x.shape) == 4:
        return GlobalAveragePooling2D()(x)
    elif len(x.shape) == 5:
        return GlobalAveragePooling3D()(x)
