import logging
from typing import Dict, List, Tuple, Sequence

import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Dense, Flatten, Reshape, LayerNormalization, DepthwiseConv2D, concatenate, Concatenate, Add

from ml4h.models.Block import Block
from ml4h.TensorMap import TensorMap
from ml4h.models.basic_blocks import DenseBlock
from ml4h.models.layer_wrappers import _upsampler, _activation_layer, _regularization_layer, _normalization_layer
from ml4h.models.layer_wrappers import _conv_layer_from_kind_and_dimension, _pool_layers_from_kind_and_dimension, _one_by_n_kernel

Tensor = tf.Tensor


class ConvEncoderBlock(Block):
    def __init__(
            self,
            *,
            tensor_map: TensorMap,
            dense_blocks: List[int],
            dense_layers: List[int],
            dense_normalize: str,
            dense_regularize: str,
            dense_regularize_rate: float,
            conv_layers: List[int],
            conv_type: str,
            conv_width: List[int],
            conv_x: List[int],
            conv_y: List[int],
            conv_z: List[int],
            block_size: int,
            activation: str,
            conv_normalize: str,
            conv_regularize: str,
            conv_regularize_rate: float,
            conv_dilate: bool,
            pool_type: str,
            pool_x: int,
            pool_y: int,
            pool_z: int,
            **kwargs,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        dimension = self.tensor_map.axes()

        # list of filter dimensions should match the total number of convolutional layers
        x_filters = _repeat_dimension(conv_width if dimension == 2 else conv_x, len(conv_layers)+len(dense_blocks))
        y_filters = _repeat_dimension(conv_y, len(conv_layers)+len(dense_blocks))
        z_filters = _repeat_dimension(conv_z, len(conv_layers)+len(dense_blocks))

        #self.preprocess_block = PreprocessBlock(['rotate'], [0.3])
        self.res_block = Residual(
            dimension=dimension, filters_per_conv=conv_layers, conv_layer_type=conv_type, conv_x=x_filters[:len(conv_layers)],
            conv_y=y_filters[:len(conv_layers)], conv_z=z_filters[:len(conv_layers)], activation=activation, normalization=conv_normalize,
            regularization=conv_regularize, regularization_rate=conv_regularize_rate, dilate=conv_dilate,
        )

        self.dense_blocks = [
            DenseConvolutional(
                dimension=dimension, conv_layer_type=conv_type, filters=filters, conv_x=[x] * block_size, conv_y=[y] * block_size,
                conv_z=[z]*block_size, block_size=block_size, activation=activation, normalization=conv_normalize,
                regularization=conv_regularize, regularization_rate=conv_regularize_rate,
            ) for filters, x, y, z in zip(dense_blocks, x_filters[len(conv_layers):], y_filters[len(conv_layers):], z_filters[len(conv_layers):])
        ]
        self.pools = _pool_layers_from_kind_and_dimension(dimension, pool_type, len(dense_blocks) + 1, pool_x, pool_y, pool_z)
        self.fully_connected = DenseBlock(
            widths=dense_layers,
            activation=activation,
            normalization=dense_normalize,
            regularization=dense_regularize,
            regularization_rate=dense_regularize_rate,
            name=self.tensor_map.embed_name(),
        ) if dense_layers else None

    def can_apply(self):
        return self.tensor_map.axes() > 1

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        if not self.can_apply():
            return x
        #x = self.preprocess_block(x)  # TODO: upgrade to tensorflow 2.3+
        x = self.res_block(x)
        intermediates[self.tensor_map].append(x)
        for i, (dense_block, pool) in enumerate(zip(self.dense_blocks, self.pools)):
            x = pool(x)
            x = dense_block(x)
            intermediates[self.tensor_map].append(x)
        if self.fully_connected:
            x = Flatten()(x)
            x = self.fully_connected(x, intermediates)
            intermediates[self.tensor_map].append(x)
        return x


class ConvDecoderBlock(Block):
    def __init__(
            self,
            *,
            tensor_map: TensorMap,
            dense_blocks: List[int],
            conv_type: str,
            conv_width: List[int],
            conv_x: List[int],
            conv_y: List[int],
            conv_z: List[int],
            block_size: int,
            activation: str,
            conv_normalize: str,
            conv_regularize: str,
            conv_regularize_rate: float,
            pool_x: int,
            pool_y: int,
            pool_z: int,
            u_connect_parents: List[TensorMap] = None,
            **kwargs,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        dimension = tensor_map.axes()
        x_filters = _repeat_dimension(conv_width if dimension == 2 else conv_x, len(dense_blocks))
        y_filters = _repeat_dimension(conv_y, len(dense_blocks))
        z_filters = _repeat_dimension(conv_z, len(dense_blocks))
        self.dense_conv_blocks = [
            DenseConvolutional(
                dimension=tensor_map.axes(), conv_layer_type=conv_type, filters=filters, conv_x=[x] * block_size,
                conv_y=[y]*block_size, conv_z=[z]*block_size, block_size=block_size, activation=activation, normalization=conv_normalize,
                regularization=conv_regularize, regularization_rate=conv_regularize_rate,
            )
            for filters, x, y, z in zip(dense_blocks, x_filters, y_filters, z_filters)
        ]
        conv_layer, _ = _conv_layer_from_kind_and_dimension(dimension, 'conv', conv_x, conv_y, conv_z)
        self.conv_label = conv_layer(tensor_map.shape[-1], _one_by_n_kernel(dimension), activation=tensor_map.activation, name=tensor_map.output_name())
        self.upsamples = [_upsampler(dimension, pool_x, pool_y, pool_z) for _ in range(len(dense_blocks) + 1)]
        self.u_connect_parents = u_connect_parents or []
        self.start_shape = _start_shape_before_pooling(num_upsamples=len(dense_blocks), output_shape=tensor_map.shape,
                                                       upsample_rates=[pool_x, pool_y, pool_z], channels=dense_blocks[-1])
        self.reshape = FlatToStructure(output_shape=self.start_shape, activation=activation, normalization=conv_normalize)
        logging.info(f'Built a decoder with: {len(self.dense_conv_blocks)} and reshape {self.start_shape}')

    def can_apply(self):
        return self.tensor_map.axes() > 1

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        if not self.can_apply():
            return x
        if x.shape != self.start_shape:
            x = self.reshape(x)
        for i, (dense_block, upsample) in enumerate(zip(self.dense_conv_blocks, self.upsamples)):
            intermediate = [intermediates[tm][len(self.upsamples)-(i+1)] for tm in self.u_connect_parents]
            x = concatenate(intermediate + [x]) if intermediate else x
            x = upsample(x)
            x = dense_block(x)
        intermediate = [intermediates[tm][0] for tm in self.u_connect_parents]
        x = concatenate(intermediate + [x]) if intermediate else x
        return self.conv_label(x)


class ResidualBlock(Block):
    def __init__(
            self,
            *,
            tensor_map: TensorMap,
            conv_layers: List[int],
            conv_type: str,
            conv_width: List[int],
            conv_x: List[int],
            conv_y: List[int],
            conv_z: List[int],
            activation: str,
            conv_normalize: str,
            conv_regularize: str,
            conv_regularize_rate: float,
            conv_dilate: bool,
            **kwargs,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        block_size = len(conv_layers)
        x_filters, y_filters, z_filters = _get_xyz_filters(block_size, conv_x if self.tensor_map.axes() > 2 else conv_width, conv_y, conv_z)
        conv_layer, kernels = _conv_layer_from_kind_and_dimension(self.tensor_map.axes(), conv_type, x_filters, y_filters, z_filters)
        self.conv_layers = []
        for i, (num_filters, kernel) in enumerate(zip(conv_layers, kernels)):
            if isinstance(conv_layer, DepthwiseConv2D):
                self.conv_layers.append(conv_layer(kernel_size=kernel, padding='same', dilation_rate=2 ** i if conv_dilate else 1))
            else:
                self.conv_layers.append(conv_layer(filters=num_filters, kernel_size=kernel, padding='same', dilation_rate=2**i if conv_dilate else 1))

        self.activations = [_activation_layer(activation) for _ in range(block_size)]
        self.normalizations = [_normalization_layer(conv_normalize) for _ in range(block_size)]
        self.regularizations = [_regularization_layer(self.tensor_map.axes(), conv_regularize, conv_regularize_rate) for _ in range(block_size)]
        residual_conv_layer, _ = _conv_layer_from_kind_and_dimension(self.tensor_map.axes(), 'conv', x_filters, y_filters, z_filters)
        self.residual_convs = [residual_conv_layer(filters=conv_layers[0], kernel_size=_one_by_n_kernel(self.tensor_map.axes())) for _ in range(block_size - 1)]
        logging.info(f'Residual Block Convolutional Layers (num_filters, kernel_size): {list(zip(conv_layers, kernels))}')

    def can_apply(self):
        return self.tensor_map.axes() > 1

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        if not self.can_apply():
            return x
        previous = x
        for convolve, activate, normalize, regularize, one_by_n_convolve in zip(
                self.conv_layers, self.activations, self.normalizations, self.regularizations, [None] + self.residual_convs,
        ):
            x = regularize(normalize(activate(convolve(x))))
            if one_by_n_convolve is not None:  # Do not residual add the input
                x = Add()([one_by_n_convolve(x), previous])
            previous = x
        intermediates[self.tensor_map].append(x)
        return x


class PoolBlock(Block):
    def __init__(
            self,
            *,
            tensor_map: TensorMap,
            pool_type: str,
            pool_x: int,
            pool_y: int,
            pool_z: int,
            **kwargs,
    ):
        self.tensor_map = tensor_map
        if not self.can_apply():
            return
        dimension = self.tensor_map.axes()
        self.pool = _pool_layers_from_kind_and_dimension(dimension, pool_type, 1, pool_x, pool_y, pool_z)[0]

    def can_apply(self):
        return self.tensor_map.axes() > 1

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]] = None) -> Tensor:
        if not self.can_apply():
            return x
        x = self.pool(x)
        intermediates[self.tensor_map].append(x)
        return x


class Residual:
    def __init__(
            self,
            *,
            dimension: int,
            filters_per_conv: List[int],
            conv_layer_type: str,
            conv_x: List[int],
            conv_y: List[int],
            conv_z: List[int],
            activation: str,
            normalization: str,
            regularization: str,
            regularization_rate: float,
            dilate: bool,
    ):
        block_size = len(filters_per_conv)
        assert len(conv_x) == len(conv_y) == len(conv_z) == block_size
        conv_layer, kernels = _conv_layer_from_kind_and_dimension(dimension, conv_layer_type, conv_x, conv_y, conv_z)
        self.conv_layers = []
        for i, (num_filters, kernel) in enumerate(zip(filters_per_conv, kernels)):
            if isinstance(conv_layer, DepthwiseConv2D):
                self.conv_layers.append(conv_layer(kernel_size=kernel, padding='same', dilation_rate=2 ** i if dilate else 1))
            else:
                self.conv_layers.append(conv_layer(filters=num_filters, kernel_size=kernel, padding='same', dilation_rate=2**i if dilate else 1))

        self.activations = [_activation_layer(activation) for _ in range(block_size)]
        self.normalizations = [_normalization_layer(normalization) for _ in range(block_size)]
        self.regularizations = [_regularization_layer(dimension, regularization, regularization_rate) for _ in range(block_size)]
        residual_conv_layer, _ = _conv_layer_from_kind_and_dimension(dimension, 'conv', conv_x, conv_y, conv_z)
        self.residual_convs = [residual_conv_layer(filters=filters_per_conv[0], kernel_size=_one_by_n_kernel(dimension)) for _ in range(block_size - 1)]
        logging.info(f'Residual Block Convolutional Layers (num_filters, kernel_size): {list(zip(filters_per_conv, kernels))}')

    def __call__(self, x: Tensor) -> Tensor:
        previous = x
        for convolve, activate, normalize, regularize, one_by_n_convolve in zip(
                self.conv_layers, self.activations, self.normalizations, self.regularizations, [None] + self.residual_convs,
        ):
            x = regularize(normalize(activate(convolve(x))))
            if one_by_n_convolve is not None:  # Do not residual add the input
                x = Add()([one_by_n_convolve(x), previous])
            previous = x
        return x


class DenseConvolutional:
    def __init__(
            self,
            *,
            dimension: int,
            block_size: int,
            conv_layer_type: str,
            filters: int,
            conv_x: List[int],
            conv_y: List[int],
            conv_z: List[int],
            activation: str,
            normalization: str,
            regularization: str,
            regularization_rate: float,
    ):
        conv_layer, kernels = _conv_layer_from_kind_and_dimension(dimension, conv_layer_type, conv_x, conv_y, conv_z)
        if isinstance(conv_layer, DepthwiseConv2D):
            self.conv_layers = [conv_layer(kernel_size=kernel, padding='same') for kernel in kernels]
        else:
            self.conv_layers = [conv_layer(filters=filters, kernel_size=kernel, padding='same') for kernel in kernels]
        self.activations = [_activation_layer(activation) for _ in range(block_size)]
        self.normalizations = [_normalization_layer(normalization) for _ in range(block_size)]
        self.regularizations = [_regularization_layer(dimension, regularization, regularization_rate) for _ in range(block_size)]
        logging.info(f'Dense Block Convolutional Layers (num_filters, kernel_size): {list(zip([filters]*len(kernels), kernels))}')

    def __call__(self, x: Tensor) -> Tensor:
        dense_connections = [x]
        for i, (convolve, activate, normalize, regularize) in enumerate(
            zip(
                    self.conv_layers, self.activations, self.normalizations, self.regularizations,
            ),
        ):
            x = normalize(regularize(activate(convolve(x))))
            if i < len(self.conv_layers) - 1:  # output of block does not get concatenated to
                dense_connections.append(x)
                x = Concatenate()(dense_connections[:])  # [:] is necessary because of tf weirdness
        return x


class FlatToStructure:
    """Takes a flat input, applies a dense layer, then restructures to output_shape"""
    def __init__(
            self,
            output_shape: Tuple[int, ...],
            activation: str,
            normalization: str,
    ):
        self.input_shapes = output_shape
        self.dense = Dense(units=int(np.prod(output_shape)))
        self.activation = _activation_layer(activation)
        self.reshape = Reshape(output_shape)
        self.norm = _normalization_layer(normalization)

    def __call__(self, x: Tensor) -> Tensor:
        return self.reshape(self.norm(self.activation(self.dense(x))))


def _repeat_dimension(filters: List[int], num_filters_needed: int) -> List[int]:
    if len(filters) < num_filters_needed:
        repeat = num_filters_needed // len(filters) + 1
        filters = (filters * repeat)[:num_filters_needed]
    return filters


def _start_shape_before_pooling(
        num_upsamples: int, output_shape: Tuple[int, ...], upsample_rates: Sequence[int], channels: int,
) -> Tuple[int, ...]:
    """
    Given the number of blocks in the decoder and the upsample rates, return required input shape to get to output shape
    """
    upsample_rates = list(upsample_rates) + [1] * len(output_shape)
    return tuple((shape // rate**num_upsamples for shape, rate in zip(output_shape[:-1], upsample_rates))) + (channels,)


def _get_xyz_filters(num_filters, conv_x, conv_y, conv_z):
    # list of filter dimensions should match the total number of convolutional layers
    x_filters = _repeat_dimension(conv_x, num_filters)
    y_filters = _repeat_dimension(conv_y, num_filters)
    z_filters = _repeat_dimension(conv_z, num_filters)
    return x_filters, y_filters, z_filters
