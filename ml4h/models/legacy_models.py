# legacy_models.py
# This file defines model factories.
# Model factories connect input TensorMaps to output TensorMaps with computational graphs.

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

# Keras imports
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
#from tensorflow.keras.layers.experimental.preprocessing import RandomRotation, RandomZoom, RandomContrast
import tensorflow_probability as tfp

from ml4h.metrics import get_metric_dict
from ml4h.plots import plot_metric_history
from ml4h.TensorMap import TensorMap, Interpretation
from ml4h.optimizers import get_optimizer, NON_KERAS_OPTIMIZERS
from ml4h.defines import JOIN_CHAR, IMAGE_EXT, MODEL_EXT, ECG_CHAR_2_IDX, PARTNERS_CHAR_2_IDX, PARTNERS_READ_TEXT

CHANNEL_AXIS = -1  # Set to 1 for Theano backend
LANGUAGE_MODEL_SUFFIX = '_next_character'
tfd = tfp.distributions


class BottleneckType(Enum):
    FlattenRestructure = auto()  # All decoder outputs are flattened to put into embedding
    GlobalAveragePoolStructured = auto()  # Structured (not flat) decoder outputs are global average pooled
    Variational = auto()  # All decoder outputs are flattened then variationally sampled to put into embedding
    NoBottleNeck = auto()  # only works when everything is u_connected
    PairLoss = auto()  # Distance between paired modalities is added to the loss


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
# PREPROCESS_CLASSES = {
#     'zoom': RandomZoom,
#     'rotate': RandomRotation,
#     'contrast': RandomContrast,
# }
CONV_REGULARIZATION_CLASSES = {
    # class name -> (dimension -> class)
    'spatial_dropout': {2: SpatialDropout1D, 3: SpatialDropout2D, 4: SpatialDropout3D},
    'dropout': defaultdict(lambda _: Dropout),
}
DENSE_REGULARIZATION_CLASSES = {
    'dropout': Dropout,  # TODO: add l1, l2
}


def make_shallow_model(
    tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap],
    learning_rate: float, model_file: str = None, model_layers: str = None,
) -> Model:
    """Make a shallow model (e.g. linear or logistic regression)

    Input and output tensor maps are set from the command line.
    Model summary printed to output

    :param tensor_maps_in: List of input TensorMaps, only 1 input TensorMap is currently supported,
                            otherwise there are layer name collisions.
    :param tensor_maps_out: List of output TensorMaps
    :param learning_rate: Size of learning steps in SGD optimization
    :param model_file: Optional HD5 model file to load and return.
    :param model_layers: Optional HD5 model file whose weights will be loaded into this model when layer names match.
    :return: a compiled keras model
    """
    if model_file is not None:
        m = load_model(model_file, custom_objects=get_metric_dict(tensor_maps_out))
        m.summary()
        logging.info("Loaded model file from: {}".format(model_file))
        return m

    losses = []
    outputs = []
    my_metrics = {}
    loss_weights = []
    input_tensors = [Input(shape=tm.shape, name=tm.input_name()) for tm in tensor_maps_in]
    it = concatenate(input_tensors) if len(input_tensors) > 1 else input_tensors[0]

    for ot in tensor_maps_out:
        losses.append(ot.loss)
        loss_weights.append(ot.loss_weight)
        my_metrics[ot.output_name()] = ot.metrics
        outputs.append(Dense(units=len(ot.channel_map), activation=ot.activation, name=ot.output_name())(it))

    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    m = Model(inputs=input_tensors, outputs=outputs)
    m.compile(optimizer=opt, loss=losses, loss_weights=loss_weights, metrics=my_metrics)
    m.summary()

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        logging.info('Loaded model weights from:{}'.format(model_layers))

    return m


def make_waveform_model_unet(
    tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], learning_rate: float,
    model_file: str = None, model_layers: str = None,
) -> Model:
    """Make a waveform predicting model

    Input and output tensor maps are set from the command line.
    Model summary printed to output

    :param tensor_maps_in: List of input TensorMaps, only 1 input TensorMap is currently supported,
                            otherwise there are layer name collisions.
    :param tensor_maps_out: List of output TensorMaps
    :param learning_rate: Size of learning steps in SGD optimization
    :param model_file: Optional HD5 model file to load and return.
    :param model_layers: Optional HD5 model file whose weights will be loaded into this model when layer names match.
    :return: a compiled keras model
    """
    if model_file is not None:
        m = load_model(model_file, custom_objects=get_metric_dict(tensor_maps_out))
        m.summary()
        logging.info("Loaded model file from: {}".format(model_file))
        return m

    neurons = 24
    input_tensor = residual = Input(shape=tensor_maps_in[0].shape, name=tensor_maps_in[0].input_name())
    x = c600 = Conv1D(filters=neurons, kernel_size=11, activation='relu', padding='same')(input_tensor)
    x = Conv1D(filters=neurons, kernel_size=51, activation='relu', padding='same')(x)
    x = MaxPooling1D(2)(x)
    x = c300 = Conv1D(filters=neurons, kernel_size=111, activation='relu', padding='same')(x)
    x = Conv1D(filters=neurons, kernel_size=201, activation='relu', padding='same')(x)
    x = MaxPooling1D(2)(x)
    x = Conv1D(filters=neurons, kernel_size=301, activation='relu', padding='same')(x)
    x = Conv1D(filters=neurons, kernel_size=301, activation='relu', padding='same')(x)
    x = UpSampling1D(2)(x)
    x = concatenate([x, c300])
    x = Conv1D(filters=neurons, kernel_size=201, activation='relu', padding='same')(x)
    x = UpSampling1D(2)(x)
    x = Conv1D(filters=neurons, kernel_size=111, activation='relu', padding='same')(x)
    x = concatenate([x, c600])
    x = Conv1D(filters=neurons, kernel_size=51, activation='relu', padding='same')(x)
    x = concatenate([x, residual])
    conv_label = Conv1D(filters=tensor_maps_out[0].shape[CHANNEL_AXIS], kernel_size=1, activation="linear")(x)
    output_y = Activation(tensor_maps_out[0].activation, name=tensor_maps_out[0].output_name())(conv_label)
    m = Model(inputs=[input_tensor], outputs=[output_y])
    m.summary()
    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0)
    m.compile(optimizer=opt, loss=tensor_maps_out[0].loss, metrics=tensor_maps_out[0].metrics)

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        logging.info('Loaded model weights from:{}'.format(model_layers))

    return m


def make_character_model_plus(
    tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], learning_rate: float, base_model: Model, language_layer: str,
    language_prefix: str, model_layers: str = None,
) -> Tuple[Model, Model]:
    """Make a ECG captioning model from an ECG embedding model

    The base_model must have an embedding layer, but besides that can have any number of other predicition TensorMaps.
    Input and output tensor maps are set from the command line.
    Model summary printed to output

    :param tensor_maps_in: List of input TensorMaps, only 1 input TensorMap is currently supported,
                            otherwise there are layer name collisions.
    :param tensor_maps_out: List of output TensorMaps
    :param learning_rate: Size of learning steps in SGD optimization
    :param base_model: The model the computes the ECG embedding
    :param language_layer: The name of TensorMap for the language string to learn
    :param language_prefix: The path prefix of the TensorMap of the language to learn
    :param model_layers: Optional HD5 model file whose weights will be loaded into this model when layer names match.
    :return: a tuple of the compiled keras model and the character emitting sub-model
    """
    char_maps_in, char_maps_out = _get_tensor_maps_for_characters(tensor_maps_in, base_model, language_layer, language_prefix)
    tensor_maps_in.extend(char_maps_in)
    tensor_maps_out.extend(char_maps_out)
    char_model = make_character_model(tensor_maps_in, tensor_maps_out, learning_rate, language_layer, model_layers=model_layers)
    losses = []
    my_metrics = {}
    loss_weights = []
    output_layers = []
    for tm in tensor_maps_out:
        losses.append(tm.loss)
        loss_weights.append(tm.loss_weight)
        my_metrics[tm.output_name()] = tm.metrics
        if tm.name == f'{language_layer}{LANGUAGE_MODEL_SUFFIX}':
            output_layers.append(char_model.get_layer(tm.output_name()))
        else:
            output_layers.append(base_model.get_layer(tm.output_name()))

    m = Model(inputs=base_model.inputs + char_model.inputs, outputs=base_model.outputs + char_model.outputs)
    m.summary()
    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0)
    m.compile(optimizer=opt, loss=losses, loss_weights=loss_weights, metrics=my_metrics)

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        _plot_dot_model_in_color(model_to_dot(m, show_shapes=True, expand_nested=True), model_layers.replace(MODEL_EXT, IMAGE_EXT), True)
        logging.info(f'Loaded and plotted model weights from:{model_layers}')

    return m, char_model


def make_character_model(
    tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], learning_rate: float,
    language_layer: str, model_file: str = None, model_layers: str = None,
) -> Model:
    """Make a ECG captioning model

    Input and output tensor maps are set from the command line. Model summary is logged

    :param tensor_maps_in: List of input TensorMaps, only 1 input TensorMap is currently supported, otherwise there are layer name collisions.
    :param tensor_maps_out: List of output TensorMaps
    :param learning_rate: Size of learning steps in SGD optimization
    :param language_layer: The name of TensorMap for the language string to learn
    :param model_file: Optional HD5 model file to load and return.
    :param model_layers: Optional HD5 model file whose weights will be loaded into this model when layer names match.
    :return: a compiled keras model
    """
    if model_file is not None:
        m = load_model(model_file, custom_objects=get_metric_dict(tensor_maps_out))
        m.summary()
        logging.info(f'Loaded model file from: {model_file}')
        return m

    input_layers = []
    for it in tensor_maps_in:
        if it.is_embedding():
            embed_in = Input(shape=it.shape, name=it.input_name())
            input_layers.append(embed_in)
        elif it.is_language():
            burn_in = Input(shape=it.shape, name=it.input_name())
            input_layers.append(burn_in)
            repeater = RepeatVector(it.shape[0])
        else:
            logging.warning(f"character model can not handle input TensorMap:{it.name} with interpretation:{it.interpretation}")

    logging.info(f"inputs: {[il.name for il in input_layers]}")
    wave_embeds = repeater(embed_in)
    lstm_in = concatenate([burn_in, wave_embeds], name='concat_embed_and_text')
    lstm_out = LSTM(128)(lstm_in)  # TODO this should be argument

    output_layers = []
    for ot in tensor_maps_out:
        if ot.name == f'{language_layer}{LANGUAGE_MODEL_SUFFIX}':
            output_layers.append(Dense(ot.shape[-1], activation=ot.activation, name=ot.output_name())(lstm_out))

    m = Model(inputs=input_layers, outputs=output_layers)
    m.summary()
    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0)
    m.compile(optimizer=opt, loss='categorical_crossentropy')

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        logging.info(f'Loaded model weights from:{model_layers}')

    return m


def make_siamese_model(
    base_model: Model,
    tensor_maps_in: List[TensorMap],
    hidden_layer: str,
    learning_rate: float = None,
    optimizer: str = 'adam',
    **kwargs
) -> Model:
    in_left = [Input(shape=tm.shape, name=tm.input_name() + '_left') for tm in tensor_maps_in]
    in_right = [Input(shape=tm.shape, name=tm.input_name() + '_right') for tm in tensor_maps_in]
    encode_model = make_hidden_layer_model(base_model, tensor_maps_in, hidden_layer)
    h_left = encode_model(in_left)
    h_right = encode_model(in_right)

    # Compute the L1 distance
    l1_layer = Lambda(lambda tensors: K.abs(tensors[0] - tensors[1]))
    l1_distance = l1_layer([h_left, h_right])

    # Add a dense layer with a sigmoid unit to generate the similarity score
    prediction = Dense(1, activation='sigmoid', name='output_siamese')(l1_distance)

    m = Model(inputs=in_left + in_right, outputs=prediction)
    opt = get_optimizer(optimizer, learning_rate, kwargs.get('optimizer_kwargs'))
    m.compile(optimizer=opt, loss='binary_crossentropy')

    if kwargs['model_layers'] is not None:
        m.load_weights(kwargs['model_layers'], by_name=True)
        logging.info(f"Loaded model weights from:{kwargs['model_layers']}")

    return m


def make_hidden_layer_model_from_file(parent_file: str, tensor_maps_in: List[TensorMap], output_layer_name: str, tensor_maps_out: List[TensorMap]) -> Model:
    parent_model = load_model(parent_file, custom_objects=get_metric_dict(tensor_maps_out))
    return make_hidden_layer_model(parent_model, tensor_maps_in, output_layer_name)


def make_hidden_layer_model(parent_model: Model, tensor_maps_in: List[TensorMap], output_layer_name: str) -> Model:
    target_layer = None
    # TODO: handle more nested models?
    for layer in parent_model.layers:
        if isinstance(layer, Model):
            try:
                target_layer = layer.get_layer(output_layer_name)
                parent_model = layer
                break
            except ValueError:
                continue
    else:
        target_layer = parent_model.get_layer(output_layer_name)
    parent_inputs = [parent_model.get_layer(tm.input_name()).input for tm in tensor_maps_in]
    dummy_input = {tm.input_name(): np.zeros((1,) + parent_model.get_layer(tm.input_name()).input_shape[0][1:]) for tm in tensor_maps_in}
    intermediate_layer_model = Model(inputs=parent_inputs, outputs=target_layer.output)
    # If we do not predict here then the graph is disconnected, I do not know why?!
    intermediate_layer_model.predict(dummy_input)
    return intermediate_layer_model


Tensor = tf.Tensor
Encoder = Callable[[Tensor], Tuple[Tensor, List[Tensor]]]
Decoder = Callable[[Tensor, Dict[TensorMap, List[Tensor]], Dict[TensorMap, Tensor]], Tensor]
BottleNeck = Callable[[Dict[TensorMap, Tensor]], Dict[TensorMap, Tensor]]


class Preprocess:
    def __init__(
            self,
            *,
            augmentations: List[str],
            factors: List[float],
    ):
        self.preprocess = [PREPROCESS_CLASSES[augmenter](factor) for augmenter, factor in zip(augmentations, factors)]

    def __call__(self, x: Tensor) -> Tensor:
        for augmentation in self.preprocess:
            x = augmentation(x)
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


class FullyConnected:
    def __init__(
            self,
            *,
            widths: List[int],
            activation: str,
            normalization: str,
            regularization: str,
            regularization_rate: float,
            is_encoder: bool = False,
            name: str = None,
            parents: List[TensorMap] = None,
    ):
        final_dense = Dense(units=widths[-1], name=name) if name else Dense(units=widths[-1])
        self.denses = [Dense(units=width) for width in widths[:-1]] + [final_dense]
        self.activations = [_activation_layer(activation) for _ in widths]
        self.regularizations = [_regularization_layer(1, regularization, regularization_rate) for _ in widths]
        self.norms = [_normalization_layer(normalization) for _ in widths]
        self.is_encoder = is_encoder
        self.parents = parents or []

    def __call__(self, x: Tensor) -> Union[Tensor, Tuple[Tensor, List[Tensor]]]:
        for dense, normalize, activate, regularize in zip(self.denses, self.norms, self.activations, self.regularizations):
            x = normalize(regularize(activate(dense(x))))
        if self.is_encoder:
            return x, []
        return x


class LSTMEncoder:
    def __init__(
            self,
            tensor_map_in,
    ):
        self.lstm = LSTM(tensor_map_in.annotation_units)

    def __call__(self, x: Tensor) -> Union[Tensor, Tuple[Tensor, List[Tensor]]]:
        return self.lstm(x), []


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


def check_no_bottleneck(u_connect: DefaultDict[TensorMap, Set[TensorMap]], tensor_maps_out: List[TensorMap]) -> bool:
    """Checks if every output tensor is u-connected to"""
    return all(any(tm in ucon_out for ucon_out in u_connect.values()) for tm in tensor_maps_out)


class UConnectBottleNeck:

    def __init__(
            self,
            u_connect: DefaultDict[TensorMap, Set[TensorMap]],
    ):
        self.u_connect = u_connect

    def __call__(self, encoder_outputs: Dict[TensorMap, Tensor]) -> Dict[TensorMap, Tensor]:
        out = {}
        for tmap_in, tensor in encoder_outputs.items():
            out = {
                **out,
                **{tmap_out: tensor for tmap_out in self.u_connect[tmap_in]},
            }
        return out


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


class VariationalBottleNeck:
    def __init__(
            self,
            activation: str,
            normalization: str,
            fully_connected_widths: List[int],
            latent_size: int,
            regularization: str,
            regularization_rate: float,
            pre_decoder_shapes: Dict[TensorMap, Optional[Tuple[int, ...]]],
    ):
        self.fully_connected = FullyConnected(
            widths=fully_connected_widths,
            activation=activation,
            normalization=normalization,
            regularization=regularization,
            regularization_rate=regularization_rate,
        ) if fully_connected_widths else None
        self.restructures = {
            tm: FlatToStructure(output_shape=shape, activation=activation, normalization=normalization)
            for tm, shape in pre_decoder_shapes.items() if shape is not None
        }
        self.latent_size = latent_size
        self.sampler: Layer = VariationalDiagNormal(latent_size, regularization_rate)
        self.no_restructures = [tm for tm, shape in pre_decoder_shapes.items() if shape is None]

    def __call__(self, encoder_outputs: Dict[TensorMap, Tensor]) -> Dict[TensorMap, Tensor]:
        y = [Flatten()(x) for x in encoder_outputs.values()]
        if len(y) > 1:
            y = concatenate(y)
        else:
            y = y[0]
        y = self.fully_connected(y) if self.fully_connected else y
        mu = Dense(self.latent_size, name='embed')(y)
        log_sigma = Dense(self.latent_size, name='log_sigma')(y)
        y = self.sampler(mu, log_sigma)
        return {
            **{tm: restructure(y) for tm, restructure in self.restructures.items()},
            **{tm: y for tm in self.no_restructures},
        }


class ConcatenateRestructure:
    """
    Flattens or GAPs then concatenates all inputs, applies a dense layer, then restructures to provided shapes
    """
    def __init__(
            self,
            pre_decoder_shapes: Dict[TensorMap, Optional[Tuple[int, ...]]],
            activation: str,
            normalization: str,
            widths: List[int],
            regularization: str,
            regularization_rate: float,
            u_connect: DefaultDict[TensorMap, Set[TensorMap]],
            bottleneck_type: BottleneckType,
    ):
        self.fully_connected = FullyConnected(
            widths=widths,
            activation=activation,
            normalization=normalization,
            regularization=regularization,
            regularization_rate=regularization_rate,
            name='embed',
        ) if widths else None
        self.restructures = {
            tm: FlatToStructure(output_shape=shape, activation=activation, normalization=normalization)
            for tm, shape in pre_decoder_shapes.items() if shape is not None
        }
        self.no_restructures = [tm for tm, shape in pre_decoder_shapes.items() if shape is None]
        self.u_connect = u_connect
        self.bottleneck_type = bottleneck_type

    def __call__(self, encoder_outputs: Dict[TensorMap, Tensor]) -> Dict[TensorMap, Tensor]:
        if self.bottleneck_type == BottleneckType.FlattenRestructure:
            y = [Flatten()(x) for x in encoder_outputs.values()]
        elif self.bottleneck_type == BottleneckType.GlobalAveragePoolStructured:
            y = [Flatten()(x) for tm, x in encoder_outputs.items() if len(x.shape) == 2]  # Flat tensors
            y += [global_average_pool(x) for tm, x in encoder_outputs.items() if len(x.shape) > 2]  # Structured tensors
        else:
            raise NotImplementedError(f'bottleneck_type {self.bottleneck_type} does not exist.')
        if len(y) > 1:
            y = concatenate(y)
        else:
            y = y[0]
        y = self.fully_connected(y) if self.fully_connected else y
        outputs: Dict[TensorMap, Tensor] = {}
        for input_tm, output_tms in self.u_connect.items():
            for output_tm in output_tms:
                outputs[output_tm] = adaptive_normalize_from_tensor(y, encoder_outputs[input_tm])
        return {
            **{tm: restructure(y) for tm, restructure in self.restructures.items()},
            **{tm: y for tm in self.no_restructures if tm not in outputs},
            **outputs,
        }


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


class ConvEncoder:

    def __init__(
            self,
            *,
            filters_per_dense_block: List[int],
            dimension: int,
            res_filters: List[int],
            conv_layer_type: str,
            conv_x: List[int],
            conv_y: List[int],
            conv_z: List[int],
            block_size: int,
            activation: str,
            normalization: str,
            regularization: str,
            regularization_rate: float,
            dilate: bool,
            pool_type: str,
            pool_x: int,
            pool_y: int,
            pool_z: int,
    ):
        num_res = len(res_filters)
        res_x, res_y, res_z = conv_x[:num_res], conv_y[:num_res], conv_z[:num_res]
        #self.preprocess_block = PreprocessBlock(['rotate'], [0.3])
        self.res_block = Residual(
            dimension=dimension, filters_per_conv=res_filters, conv_layer_type=conv_layer_type, conv_x=res_x,
            conv_y=res_y, conv_z=res_z, activation=activation, normalization=normalization,
            regularization=regularization, regularization_rate=regularization_rate, dilate=dilate,
        )

        dense_x, dense_y, dense_z = conv_x[num_res:], conv_y[num_res:], conv_z[num_res:]
        self.dense_blocks = [
            DenseConvolutional(
                dimension=dimension, conv_layer_type=conv_layer_type, filters=filters, conv_x=[x]*block_size, conv_y=[y]*block_size,
                conv_z=[z]*block_size, block_size=block_size, activation=activation, normalization=normalization,
                regularization=regularization, regularization_rate=regularization_rate,
            ) for filters, x, y, z in zip(filters_per_dense_block, dense_x, dense_y, dense_z)
        ]
        self.pools = _pool_layers_from_kind_and_dimension(dimension, pool_type, len(filters_per_dense_block) + 1, pool_x, pool_y, pool_z)

    def __call__(self, x: Tensor) -> Tuple[Tensor, List[Tensor]]:
        intermediates = []
        #x = self.preprocess_block(x)
        x = self.res_block(x)
        intermediates.append(x)
        x = self.pools[0](x)
        for i, (dense_block, pool) in enumerate(zip(self.dense_blocks, self.pools[1:])):
            x = dense_block(x)
            intermediates.append(x)
            x = pool(x) if i < len(self.dense_blocks) - 1 else x  # don't pool after final dense block
        return x, intermediates


def _calc_start_shape(
        num_upsamples: int, output_shape: Tuple[int, ...], upsample_rates: Sequence[int], channels: int,
) -> Tuple[int, ...]:
    """
    Given the number of blocks in the decoder and the upsample rates, return required input shape to get to output shape
    """
    upsample_rates = list(upsample_rates) + [1] * len(output_shape)
    return tuple((shape // rate**num_upsamples for shape, rate in zip(output_shape[:-1], upsample_rates))) + (channels,)


class DenseDecoder:
    def __init__(
            self,
            tensor_map_out: TensorMap,
            activation: str,
            parents: List[TensorMap] = None,
    ):
        self.parents = parents
        self.activation = _activation_layer(activation)
        self.dense = Dense(units=tensor_map_out.shape[0], name=tensor_map_out.output_name(), activation=tensor_map_out.activation)
        self.units = tensor_map_out.annotation_units

    def __call__(self, x: Tensor, _, decoder_outputs: Dict[TensorMap, Tensor]) -> Tensor:
        if self.parents:
            x = Concatenate()([x] + [decoder_outputs[parent] for parent in self.parents])
            x = Dense(units=self.units)(x)
            x = self.activation(x)
        return self.dense(x)


class LanguageDecoder:
    def __init__(
            self,
            tensor_map_out: TensorMap,
    ):
        self.dense = Dense(tensor_map_out.shape[-1], activation=tensor_map_out.activation, name=tensor_map_out.output_name())

    def __call__(self, x: Tensor, _, decoder_outputs: Dict[TensorMap, Tensor]) -> Tensor:
        return self.dense(x)


class ConvDecoder:
    def __init__(
            self,
            *,
            tensor_map_out: TensorMap,
            filters_per_dense_block: List[int],
            conv_layer_type: str,
            conv_x: List[int],
            conv_y: List[int],
            conv_z: List[int],
            block_size: int,
            activation: str,
            normalization: str,
            regularization: str,
            regularization_rate: float,
            upsample_x: int,
            upsample_y: int,
            upsample_z: int,
            u_connect_parents: List[TensorMap] = None,
    ):
        dimension = tensor_map_out.axes()
        self.dense_blocks = [
            DenseConvolutional(
                dimension=tensor_map_out.axes(), conv_layer_type=conv_layer_type, filters=filters, conv_x=[x]*block_size,
                conv_y=[y]*block_size, conv_z=[z]*block_size, block_size=block_size, activation=activation, normalization=normalization,
                regularization=regularization, regularization_rate=regularization_rate,
            )
            for filters, x, y, z in zip(filters_per_dense_block, conv_x, conv_y, conv_z)
        ]
        conv_layer, _ = _conv_layer_from_kind_and_dimension(dimension, 'conv', conv_x, conv_y, conv_z)
        self.conv_label = conv_layer(tensor_map_out.shape[-1], _one_by_n_kernel(dimension), activation=tensor_map_out.activation, name=tensor_map_out.output_name())
        self.upsamples = [_upsampler(dimension, upsample_x, upsample_y, upsample_z) for _ in range(len(filters_per_dense_block) + 1)]
        self.u_connect_parents = u_connect_parents or []

    def __call__(self, x: Tensor, intermediates: Dict[TensorMap, List[Tensor]], _) -> Tensor:
        for i, (dense_block, upsample) in enumerate(zip(self.dense_blocks, self.upsamples)):
            intermediate = [intermediates[tm][-(i + 1)] for tm in self.u_connect_parents]
            x = concatenate(intermediate + [x]) if intermediate else x
            x = upsample(x)
            x = dense_block(x)
        intermediate = [intermediates[tm][0] for tm in self.u_connect_parents]
        x = concatenate(intermediate + [x]) if intermediate else x
        return self.conv_label(x)


def parent_sort(tms: List[TensorMap]) -> List[TensorMap]:
    """
    Parents will always appear before their children after sorting. Idempotent and slow.
    """
    to_process = sorted(tms, key=lambda x: str(x))
    final: List[TensorMap] = []
    visited = Counter()
    while to_process:
        tm = to_process.pop()
        visited[tm] += 1
        if visited[tm] > len(tms):
            raise ValueError('Problem detected in parent structure. Could be cycle or missing parent.')
        if not tm.parents or set(tm.parents) <= set(final):
            final.append(tm)
        else:
            to_process.insert(0, tm)
    return final


def _get_custom_objects(tensor_maps_out: List[TensorMap]) -> Dict[str, Any]:
    custom_objects = {
        obj.__name__: obj
        for obj in chain(
            NON_KERAS_OPTIMIZERS.values(), ACTIVATION_FUNCTIONS.values(), NORMALIZATION_CLASSES.values(),
            [VariationalDiagNormal, L2LossLayer, CosineLossLayer],
        )
    }
    return {**custom_objects, **get_metric_dict(tensor_maps_out)}


def _repeat_dimension(filters: List[int], num_filters_needed: int) -> List[int]:
    if len(filters) < num_filters_needed:
        repeat = num_filters_needed // len(filters) + 1
        filters = (filters * repeat)[:num_filters_needed]
    return filters


def make_multimodal_multitask_model(
        tensor_maps_in: List[TensorMap],
        tensor_maps_out: List[TensorMap],
        activation: str,
        learning_rate: float,
        bottleneck_type: BottleneckType,
        optimizer: str,
        dense_layers: List[int] = None,
        dense_normalize: str = None,
        dense_regularize: str = None,
        dense_regularize_rate: float = None,
        conv_layers: List[int] = None,
        dense_blocks: List[int] = None,
        block_size: int = None,
        conv_type: str = None,
        conv_normalize: str = None,
        conv_regularize: str = None,
        conv_regularize_rate: float = None,
        conv_x: List[int] = None,
        conv_y: List[int] = None,
        conv_z: List[int] = None,
        conv_width: List[int] = None,
        conv_dilate: bool = None,
        u_connect: DefaultDict[TensorMap, Set[TensorMap]] = None,
        pool_x: int = None,
        pool_y: int = None,
        pool_z: int = None,
        pool_type: str = None,
        training_steps: int = None,
        learning_rate_schedule: str = None,
        **kwargs,
) -> Model:
    """Make multi-task, multi-modal feed forward neural network for all kinds of prediction

    This model factory can be used to make networks for classification, regression, and segmentation
    The tasks attempted are given by the output TensorMaps.
    The modalities and the first layers in the architecture are determined by the input TensorMaps.

    Hyperparameters are exposed to the command line.
    Model summary printed to output

    :param tensor_maps_in: List of input TensorMaps
    :param tensor_maps_out: List of output TensorMaps
    :param activation: Activation function as a string (e.g. 'relu', 'linear, or 'softmax)
    :param learning_rate: learning rate for optimizer
    :param bottleneck_type: How to merge the representations coming from the different input modalities
    :param dense_layers: List of number of filters in each dense layer.
    :param dense_normalize: How to normalize dense layers (e.g. batchnorm)
    :param dense_regularize: How to regularize dense leayers (e.g. dropout)
    :param dense_regularize_rate: Rate of dense_regularize
    :param conv_layers: List of number of filters in each convolutional layer
    :param dense_blocks: List of number of filters in densenet modules for densenet convolutional models
    :param block_size: Number of layers within each Densenet module for densenet convolutional models
    :param conv_type: Type of convolution to use, e.g. separable
    :param conv_normalize: Type of normalization layer for convolutions, e.g. batch norm
    :param conv_regularize: Type of regularization for convolutions (e.g. dropout)
    :param conv_regularize_rate: Rate of conv_regularize
    :param conv_width: Size of X dimension for 1D convolutional kernels
    :param conv_x: Size of X dimension for 2D and 3D convolutional kernels
    :param conv_y: Size of Y dimension for 2D and 3D convolutional kernels
    :param conv_z: Size of Z dimension for 3D convolutional kernels
    :param conv_dilate: whether to use dilation in conv layers
    :param u_connect: dictionary of input TensorMap -> output TensorMaps to u connect to
    :param pool_x: Pooling in the X dimension for Convolutional models.
    :param pool_y: Pooling in the Y dimension for Convolutional models.
    :param pool_z: Pooling in the Z dimension for 3D Convolutional models.
    :param pool_type: max or average pooling following convolutional blocks
    :param optimizer: which optimizer to use. See optimizers.py.
    :return: a compiled keras model
    :param learning_rate_schedule: learning rate schedule to train with, e.g. triangular
    :param training_steps: How many training steps to train the model. Only needed if learning_rate_schedule given
    :param model_file: HD5 model file to load and return.
    :param model_layers: HD5 model file whose weights will be loaded into this model when layer names match.
    :param freeze_model_layers: Whether to freeze layers from loaded from model_layers
    """
    tensor_maps_out = parent_sort(tensor_maps_out)
    u_connect: DefaultDict[TensorMap, Set[TensorMap]] = u_connect or defaultdict(set)
    custom_dict = _get_custom_objects(tensor_maps_out)
    opt = get_optimizer(
        optimizer, learning_rate, steps_per_epoch=training_steps, learning_rate_schedule=learning_rate_schedule,
        optimizer_kwargs=kwargs.get('optimizer_kwargs'),
    )
    if 'model_file' in kwargs and kwargs['model_file'] is not None:
        logging.info("Attempting to load model file from: {}".format(kwargs['model_file']))
        m = load_model(kwargs['model_file'], custom_objects=custom_dict, compile=False)
        m.compile(optimizer=opt, loss=custom_dict['loss'])
        m.summary()
        logging.info("Loaded model file from: {}".format(kwargs['model_file']))
        return m

    # list of filter dimensions should match the number of convolutional layers = len(dense_blocks) + [ + len(conv_layers) if convolving input tensors]
    num_dense = len(dense_blocks)
    num_res = len(conv_layers) if any(tm.axes() > 1 for tm in tensor_maps_in) else 0
    num_filters_needed = num_res + num_dense
    conv_x = _repeat_dimension(conv_x, num_filters_needed)
    conv_y = _repeat_dimension(conv_y, num_filters_needed)
    conv_z = _repeat_dimension(conv_z, num_filters_needed)

    encoders: Dict[TensorMap: Layer] = {}
    for tm in tensor_maps_in:
        if tm.is_language():
            encoders[tm] = LSTMEncoder(tm)
        elif tm.axes() > 1:
            encoders[tm] = ConvEncoder(
                filters_per_dense_block=dense_blocks,
                dimension=tm.axes(),
                res_filters=conv_layers,
                conv_layer_type=conv_type,
                conv_x=conv_x if tm.axes() > 2 else conv_width,
                conv_y=conv_y,
                conv_z=conv_z,
                block_size=block_size,
                activation=activation,
                normalization=conv_normalize,
                regularization=conv_regularize,
                regularization_rate=conv_regularize_rate,
                dilate=conv_dilate,
                pool_type=pool_type,
                pool_x=pool_x,
                pool_y=pool_y,
                pool_z=pool_z,
            )
        else:
            encoders[tm] = FullyConnected(
                widths=[tm.annotation_units],
                activation=activation,
                normalization=dense_normalize,
                regularization=dense_regularize,
                regularization_rate=dense_regularize_rate,
                is_encoder=True,
            )

    pre_decoder_shapes: Dict[TensorMap, Optional[Tuple[int, ...]]] = {}
    for tm in tensor_maps_out:
        if any([tm in out for out in u_connect.values()]) or tm.axes() == 1 or tm.is_language():
            pre_decoder_shapes[tm] = None
        else:
            pre_decoder_shapes[tm] = _calc_start_shape(
                num_upsamples=len(dense_blocks), output_shape=tm.shape, upsample_rates=[pool_x, pool_y, pool_z],
                channels=dense_blocks[-1],
            )

    if bottleneck_type in {BottleneckType.FlattenRestructure, BottleneckType.GlobalAveragePoolStructured}:
        bottleneck = ConcatenateRestructure(
            widths=dense_layers,
            activation=activation,
            regularization=dense_regularize,
            regularization_rate=dense_regularize_rate,
            normalization=dense_normalize,
            pre_decoder_shapes=pre_decoder_shapes,
            u_connect=u_connect,
            bottleneck_type=bottleneck_type,
        )
    elif bottleneck_type == BottleneckType.Variational:
        bottleneck = VariationalBottleNeck(
            fully_connected_widths=dense_layers[:-1],
            latent_size=dense_layers[-1],
            activation=activation,
            regularization=dense_regularize,
            regularization_rate=dense_regularize_rate,
            normalization=dense_normalize,
            pre_decoder_shapes=pre_decoder_shapes,
        )
    elif bottleneck_type == BottleneckType.NoBottleNeck:
        if not check_no_bottleneck(u_connect, tensor_maps_out):
            bottleneck = None
        else:
            bottleneck = UConnectBottleNeck(u_connect)
    else:
        raise NotImplementedError(f'Unknown BottleneckType {bottleneck_type}.')

    conv_x, conv_y, conv_z = conv_x[num_res:], conv_y[num_res:], conv_z[num_res:]
    decoders: Dict[TensorMap, Layer] = {}
    for tm in tensor_maps_out:
        if tm.is_language():
            decoders[tm] = LanguageDecoder(tensor_map_out=tm)
        elif tm.axes() > 1:
            decoders[tm] = ConvDecoder(
                tensor_map_out=tm,
                filters_per_dense_block=dense_blocks[::-1],
                conv_layer_type=conv_type,
                conv_x=conv_x if tm.axes() > 2 else conv_width,
                conv_y=conv_y,
                conv_z=conv_z,
                block_size=block_size,
                activation=activation,
                normalization=conv_normalize,
                regularization=conv_regularize,
                regularization_rate=conv_regularize_rate,
                upsample_x=pool_x,
                upsample_y=pool_y,
                upsample_z=pool_z,
                u_connect_parents=[tm_in for tm_in in tensor_maps_in if tm in u_connect[tm_in]],
            )
        else:
            decoders[tm] = DenseDecoder(
                tensor_map_out=tm,
                parents=tm.parents,
                activation=activation,
            )

    m = _make_multimodal_multitask_model(encoders, bottleneck, decoders)

    # load layers for transfer learning
    model_layers = kwargs.get('model_layers', False)
    if model_layers:
        loaded = 0
        freeze = kwargs.get('freeze_model_layers', False)
        m.load_weights(model_layers, by_name=True)
        try:
            m_other = load_model(model_layers, custom_objects=custom_dict, compile=False)
            for other_layer in m_other.layers:
                try:
                    target_layer = m.get_layer(other_layer.name)
                    target_layer.set_weights(other_layer.get_weights())
                    loaded += 1
                    if freeze:
                        target_layer.trainable = False
                except (ValueError, KeyError):
                    logging.warning(f'Error loading layer {other_layer.name} from model: {model_layers}. Will still try to load other layers.')
        except ValueError as e:
            logging.info(f'Loaded model weights, but got ValueError in model loading: {str(e)}')
        logging.info(f'Loaded {"and froze " if freeze else ""}{loaded} layers from {model_layers}.')
    m.compile(
        optimizer=opt, loss=[tm.loss for tm in tensor_maps_out],
        metrics={tm.output_name(): tm.metrics for tm in tensor_maps_out},
    )
    m.summary()
    return m


def _make_multimodal_multitask_model(
        encoders: Dict[TensorMap, Encoder],
        bottle_neck: BottleNeck,
        decoders: Dict[TensorMap, Decoder],  # Assumed to be topologically sorted according to parents hierarchy
) -> Model:
    inputs: Dict[TensorMap, Input] = {}
    encoder_outputs: Dict[TensorMap, Tuple[Tensor, List[Tensor]]] = {}  # TensorMap -> embed, encoder_intermediates
    encoder_intermediates = {}
    for tm, encoder in encoders.items():
        x = Input(shape=tm.shape, name=tm.input_name())
        inputs[tm] = x
        y, intermediates = encoder(x)
        encoder_outputs[tm] = y
        encoder_intermediates[tm] = intermediates

    bottle_neck_outputs = bottle_neck(encoder_outputs)

    decoder_outputs = {}
    for tm, decoder in decoders.items():
        decoder_outputs[tm] = decoder(bottle_neck_outputs[tm], encoder_intermediates, decoder_outputs)

    return Model(inputs=list(inputs.values()), outputs=list(decoder_outputs.values()))


def _transfer_layers_by_name(model_layers: str, freeze_model_layers: str, custom_dict: Dict[str, Any], m: Model):
    # load layers for transfer learning
    loaded = 0
    m.load_weights(model_layers, by_name=True)
    try:
        m_other = load_model(model_layers, custom_objects=custom_dict, compile=False)
        for other_layer in m_other.layers:
            try:
                target_layer = m.get_layer(other_layer.name)
                target_layer.set_weights(other_layer.get_weights())
                loaded += 1
                if freeze_model_layers:
                    target_layer.trainable = False
            except (ValueError, KeyError):
                logging.warning(f'Error loading layer {other_layer.name} from model: {model_layers}. Will still try to load other layers.')
    except ValueError as e:
        logging.info(f'Loaded model weights, but got ValueError in model loading: {str(e)}')
    logging.info(f'Loaded {"and froze " if freeze_model_layers else ""}{loaded} layers from {model_layers}.')


def _load_model_encoders_and_decoders(tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], custom_dict: Dict[str, Any],
                                      optimizer, model_file: str):
    encoders = {}
    decoders = {}
    merger = None
    try:
        for tm in tensor_maps_in:
            encoders[tm] = load_model(f"{os.path.dirname(model_file)}/encoder_{tm.name}.h5", custom_objects=custom_dict, compile=False)
        for tm in tensor_maps_out:
            decoders[tm] = load_model(f"{os.path.dirname(model_file)}/decoder_{tm.name}.h5", custom_objects=custom_dict, compile=False)
        merger = load_model(f"{os.path.dirname(model_file)}/merger.h5", custom_objects=custom_dict, compile=False)
    except OSError as e:
        logging.warning(f'Could not load some model modules, error: {e}')
    logging.info(f"Attempting to load model file from: {model_file}")
    m = load_model(model_file, custom_objects=custom_dict, compile=False)
    m.compile(optimizer=optimizer, loss=[tm.loss for tm in tensor_maps_out],
              metrics={tm.output_name(): tm.metrics for tm in tensor_maps_out})
    m.summary()
    logging.info(f"Loaded encoders, decoders and model file from: {model_file}")
    return m, encoders, decoders, merger


def make_paired_autoencoder_model(
        pairs: List[Tuple[TensorMap, TensorMap]],
        pair_loss: str = 'cosine',
        pair_loss_weight: float = 1.0,
        multimodal_merge: str = 'average',
        **kwargs
) -> Model:
    custom_dict = _get_custom_objects(kwargs['tensor_maps_out'])
    opt = get_optimizer(
        kwargs['optimizer'], kwargs['learning_rate'], steps_per_epoch=kwargs['training_steps'],
        learning_rate_schedule=kwargs['learning_rate_schedule'], optimizer_kwargs=kwargs.get('optimizer_kwargs'),
    )
    if 'model_file' in kwargs and kwargs['model_file'] is not None:
        return _load_model_encoders_and_decoders(kwargs['tensor_maps_in'], kwargs['tensor_maps_out'], custom_dict, opt, kwargs['model_file'])

    inputs = {tm: Input(shape=tm.shape, name=tm.input_name()) for tm in kwargs['tensor_maps_in']}
    real_serial_layers = kwargs['model_layers']
    kwargs['model_layers'] = None
    multimodal_activations = []
    encoders = {}
    decoders = {}
    outputs = {}
    losses = []
    for left, right in pairs:
        if left in encoders:
            encode_left = encoders[left]
        else:
            kwargs['tensor_maps_in'] = [left]
            left_model = make_multimodal_multitask_model(**kwargs)
            encode_left = make_hidden_layer_model(left_model, [left], kwargs['hidden_layer'])
        h_left = encode_left(inputs[left])

        if right in encoders:
            encode_right = encoders[right]
        else:
            kwargs['tensor_maps_in'] = [right]
            right_model = make_multimodal_multitask_model(**kwargs)
            encode_right = make_hidden_layer_model(right_model, [right], kwargs['hidden_layer'])
        h_right = encode_right(inputs[right])

        if pair_loss == 'cosine':
            loss_layer = CosineLossLayer(pair_loss_weight)
        elif pair_loss == 'euclid':
            loss_layer = L2LossLayer(pair_loss_weight)

        multimodal_activations.extend(loss_layer([h_left, h_right]))
        encoders[left] = encode_left
        encoders[right] = encode_right

    kwargs['tensor_maps_in'] = list(inputs.keys())
    if multimodal_merge == 'average':
        multimodal_activation = Average()(multimodal_activations)
    elif multimodal_merge == 'concatenate':
        multimodal_activation = Concatenate()(multimodal_activations)
        multimodal_activation = Dense(units=kwargs['dense_layers'][0], use_bias=False)(multimodal_activation)
        multimodal_activation = _activation_layer(kwargs['activation'])(multimodal_activation)
    else:
        raise NotImplementedError(f'No merge architecture for method: {multimodal_merge}')
    latent_inputs = Input(shape=(kwargs['dense_layers'][-1]), name='input_concept_space')

    # build decoder models
    for tm in kwargs['tensor_maps_out']:
        if tm.axes() > 1:
            shape = _calc_start_shape(num_upsamples=len(kwargs['dense_blocks']), output_shape=tm.shape,
                                      upsample_rates=[kwargs['pool_x'], kwargs['pool_y'], kwargs['pool_z']],
                                      channels=kwargs['dense_blocks'][-1])

            restructure = FlatToStructure(output_shape=shape, activation=kwargs['activation'],
                                          normalization=kwargs['dense_normalize'])

            decode = ConvDecoder(
                tensor_map_out=tm,
                filters_per_dense_block=kwargs['dense_blocks'][::-1],
                conv_layer_type=kwargs['conv_type'],
                conv_x=kwargs['conv_x'] if tm.axes() > 2 else kwargs['conv_width'],
                conv_y=kwargs['conv_y'],
                conv_z=kwargs['conv_z'],
                block_size=kwargs['block_size'],
                activation=kwargs['activation'],
                normalization=kwargs['conv_normalize'],
                regularization=kwargs['conv_regularize'],
                regularization_rate=kwargs['conv_regularize_rate'],
                upsample_x=kwargs['pool_x'],
                upsample_y=kwargs['pool_y'],
                upsample_z=kwargs['pool_z'],
                u_connect_parents=[tm_in for tm_in in kwargs['tensor_maps_in'] if tm in kwargs['u_connect'][tm_in]],
            )
            reconstruction = decode(restructure(latent_inputs), {}, {})
        else:
            dense_block = FullyConnected(
                widths=kwargs['dense_layers'],
                activation=kwargs['activation'],
                normalization=kwargs['dense_normalize'],
                regularization=kwargs['dense_regularize'],
                regularization_rate=kwargs['dense_regularize_rate'],
                is_encoder=False,
            )
            decode = DenseDecoder(tensor_map_out=tm, parents=tm.parents, activation=kwargs['activation'])
            reconstruction = decode(dense_block(latent_inputs), {}, {})

        decoder = Model(latent_inputs, reconstruction, name=tm.output_name())
        decoders[tm] = decoder
        outputs[tm.output_name()] = decoder(multimodal_activation)
        losses.append(tm.loss)

    m = Model(inputs=list(inputs.values()), outputs=list(outputs.values()))
    my_metrics = {tm.output_name(): tm.metrics for tm in kwargs['tensor_maps_out']}
    m.compile(optimizer=opt, loss=losses, metrics=my_metrics)
    m.summary()

    if real_serial_layers is not None:
        m.load_weights(real_serial_layers, by_name=True)
        logging.info(f"Loaded model weights from:{real_serial_layers}")

    return m, encoders, decoders


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Predicting ~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def embed_model_predict(model, tensor_maps_in, embed_layer, test_data, batch_size):
    embed_model = make_hidden_layer_model(model, tensor_maps_in, embed_layer)
    return embed_model.predict(test_data, batch_size=batch_size)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Model Builders ~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


def _get_tensor_maps_for_characters(
    tensor_maps_in: List[TensorMap], base_model: Model, language_layer: str, language_prefix: str, embed_name='embed',
    embed_size=64, burn_in=100,
):
    embed_model = make_hidden_layer_model(base_model, tensor_maps_in, embed_name)
    tm_embed = TensorMap(embed_name, shape=(embed_size,), interpretation=Interpretation.EMBEDDING, parents=tensor_maps_in.copy(), model=embed_model)

    if PARTNERS_READ_TEXT in language_layer:
        tm_char = TensorMap(
            f'{language_layer}{LANGUAGE_MODEL_SUFFIX}', Interpretation.LANGUAGE, shape=(len(PARTNERS_CHAR_2_IDX),),
            channel_map=PARTNERS_CHAR_2_IDX, cacheable=False,
        )
        tm_burn_in = TensorMap(
            language_layer, Interpretation.LANGUAGE, shape=(burn_in, len(PARTNERS_CHAR_2_IDX)), path_prefix=language_prefix,
            dependent_map=tm_char, cacheable=False,
        )
        logging.info(f'From language layer: {language_layer} created tensor maps for Partners language data.')
    else:
        tm_char = TensorMap(
            f'{language_layer}{LANGUAGE_MODEL_SUFFIX}', Interpretation.LANGUAGE, shape=(len(ECG_CHAR_2_IDX),), channel_map=ECG_CHAR_2_IDX,
            cacheable=False,
        )
        tm_burn_in = TensorMap(
            language_layer, Interpretation.LANGUAGE, shape=(burn_in, len(ECG_CHAR_2_IDX)), path_prefix=language_prefix,
            dependent_map=tm_char, cacheable=False,
        )
        logging.info(f'From language layer: {language_layer} created tensor maps for UKB ECG language data.')

    return [tm_embed, tm_burn_in], [tm_char]


def get_model_inputs_outputs(
    model_files: List[str],
    tensor_maps_in: List[TensorMap],
    tensor_maps_out: List[TensorMap],
) -> Dict[str, Dict[str, TensorMap]]:
    """Organizes given input and output tensors as nested dictionary.

    Returns:
        dict: The nested dictionary of tensors.
            The inner dictionary is keyed by tensor type ('input' or 'output').
            The outer dictionary is keyed by 'model_file'.

            {
                'model_file_1':
                    {
                        'input': [tensor1, tensor2],
                        'output': [tensor3, tensor4]
                    },
                'model_file_2':
                    {
                        'input': [tensor2, tensor5],
                        'output': [tensor4, tensor6]
                    }
            }

    """

    input_prefix = "input"
    output_prefix = "output"
    got_tensor_maps_for_characters = False
    models_inputs_outputs = dict()

    for model_file in model_files:
        custom = _get_custom_objects(tensor_maps_out)
        logging.info(f'custom keys: {list(custom.keys())}')
        m = load_model(model_file, custom_objects=custom, compile=False)
        model_inputs_outputs = defaultdict(list)
        for input_tensor_map in tensor_maps_in:
            try:
                m.get_layer(input_tensor_map.input_name())
                model_inputs_outputs[input_prefix].append(input_tensor_map)
            except ValueError:
                pass
        for output_tensor_map in tensor_maps_out:
            try:
                m.get_layer(output_tensor_map.output_name())
                model_inputs_outputs[output_prefix].append(output_tensor_map)
            except ValueError:
                pass
        if not got_tensor_maps_for_characters:
            try:
                m.get_layer('input_ecg_rest_text_ecg_text')
                # TODO: Is this broken?
                char_maps_in, char_maps_out = _get_tensor_maps_for_characters(tensor_maps_in, m)
                model_inputs_outputs[input_prefix].extend(char_maps_in)
                tensor_maps_in.extend(char_maps_in)
                model_inputs_outputs[output_prefix].extend(char_maps_out)
                tensor_maps_out.extend(char_maps_out)
                got_tensor_maps_for_characters = True
                logging.info(f"Doing char model dance:{[tm.input_name() for tm in tensor_maps_in]}")
                logging.info(f"Doing char model dance out:{[tm.output_name() for tm in tensor_maps_out]}")
            except ValueError:
                pass
        models_inputs_outputs[model_file] = model_inputs_outputs

    return models_inputs_outputs
