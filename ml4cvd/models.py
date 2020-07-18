# models.py
# This file defines model factories.
# Model factories connect input TensorMaps to output TensorMaps with computational graphs.

# Imports
import os
import time
import logging
import numpy as np
from enum import Enum, auto
from itertools import chain
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
from tensorflow.keras.layers import MaxPooling2D, MaxPooling3D, AveragePooling1D, AveragePooling2D, AveragePooling3D, Layer
from tensorflow.keras.layers import SeparableConv1D, SeparableConv2D, DepthwiseConv2D, Concatenate, Add
import tensorflow_probability as tfp

from ml4cvd.metrics import get_metric_dict
from ml4cvd.plots import plot_metric_history
from ml4cvd.TensorMap import TensorMap, Interpretation
from ml4cvd.optimizers import get_optimizer, NON_KERAS_OPTIMIZERS
from ml4cvd.defines import JOIN_CHAR, IMAGE_EXT, MODEL_EXT, ECG_CHAR_2_IDX, PARTNERS_CHAR_2_IDX, PARTNERS_READ_TEXT

CHANNEL_AXIS = -1  # Set to 1 for Theano backend
LANGUAGE_MODEL_SUFFIX = '_next_character'
tfd = tfp.distributions


class BottleneckType(Enum):
    FlattenRestructure = auto()  # All decoder outputs are flattened to put into embedding
    GlobalAveragePoolStructured = auto()  # Structured (not flat) decoder outputs are global average pooled
    Variational = auto()  # All decoder outputs are flattened then variationally sampled to put into embedding
    NoBottleNeck = auto()  # only works when everything is u_connected


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


class ResidualBlock:
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
        self.conv_layers = [
            conv_layer(filters=num_filters, kernel_size=kernel, padding='same', dilation_rate=2**i if dilate else 1)
            for i, (num_filters, kernel) in enumerate(zip(filters_per_conv, kernels))
        ]
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


class DenseConvolutionalBlock:
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


class FullyConnectedBlock:
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
    return K.mean(x, axis=tuple(range(1, len(x.shape) - 1)))


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
        self.fully_connected = FullyConnectedBlock(
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
        self.fully_connected = FullyConnectedBlock(
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
        self.res_block = ResidualBlock(
            dimension=dimension, filters_per_conv=res_filters, conv_layer_type=conv_layer_type, conv_x=res_x,
            conv_y=res_y, conv_z=res_z, activation=activation, normalization=normalization,
            regularization=regularization, regularization_rate=regularization_rate, dilate=dilate,
        )

        dense_x, dense_y, dense_z = conv_x[num_res:], conv_y[num_res:], conv_z[num_res:]
        self.dense_blocks = [
            DenseConvolutionalBlock(
                dimension=dimension, conv_layer_type=conv_layer_type, filters=filters, conv_x=[x]*block_size, conv_y=[y]*block_size,
                conv_z=[z]*block_size, block_size=block_size, activation=activation, normalization=normalization,
                regularization=regularization, regularization_rate=regularization_rate,
            ) for filters, x, y, z in zip(filters_per_dense_block, dense_x, dense_y, dense_z)
        ]
        self.pools = _pool_layers_from_kind_and_dimension(dimension, pool_type, len(filters_per_dense_block) + 1, pool_x, pool_y, pool_z)

    def __call__(self, x: Tensor) -> Tuple[Tensor, List[Tensor]]:
        intermediates = []
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
            DenseConvolutionalBlock(
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
            x = dense_block(x)
            x = upsample(x)
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
            [VariationalDiagNormal],
        )
    }
    return {**custom_objects, **get_metric_dict(tensor_maps_out)}


def _repeat_dimension(dim: List[int], name: str, num_filters_needed: int) -> List[int]:
        if len(dim) != num_filters_needed:
            logging.warning(
                f'Number of {name} dimensions for convolutional kernel sizes ({len(dim)}) '
                f'do not match number of convolutional layers/blocks ({num_filters_needed}), '
                f'matching values to fit {num_filters_needed} convolutional layers/blocks.',
            )
            repeat = num_filters_needed // len(dim) + 1
            dim = (dim * repeat)[:num_filters_needed]
        return dim


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
        conv_dilate: bool = None,
        u_connect: DefaultDict[TensorMap, Set[TensorMap]] = None,
        pool_x: int = None,
        pool_y: int = None,
        pool_z: int = None,
        pool_type: str = None,
        training_steps: int = None,
        learning_rate_schedule: str = None,
        **kwargs
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
    :param: dense_regularize_rate: Rate of dense_regularize
    :param conv_layers: List of number of filters in each convolutional layer
    :param dense_blocks: List of number of filters in densenet modules for densenet convolutional models
    :param block_size: Number of layers within each Densenet module for densenet convolutional models
    :param conv_type: Type of convolution to use, e.g. separable
    :param conv_normalize: Type of normalization layer for convolutions, e.g. batch norm
    :param conv_regularize: Type of regularization for convolutions (e.g. dropout)
    :param conv_regularize_rate: Rate of conv_regularize
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
    conv_x = _repeat_dimension(conv_x, 'x', num_filters_needed)
    conv_y = _repeat_dimension(conv_y, 'y', num_filters_needed)
    conv_z = _repeat_dimension(conv_z, 'z', num_filters_needed)

    encoders: Dict[TensorMap: Layer] = {}
    for tm in tensor_maps_in:
        if tm.axes() > 1:
            encoders[tm] = ConvEncoder(
                filters_per_dense_block=dense_blocks,
                dimension=tm.axes(),
                res_filters=conv_layers,
                conv_layer_type=conv_type,
                conv_x=conv_x,
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
            encoders[tm] = FullyConnectedBlock(
                widths=[tm.annotation_units],
                activation=activation,
                normalization=dense_normalize,
                regularization=dense_regularize,
                regularization_rate=dense_regularize_rate,
                is_encoder=True,
            )

    pre_decoder_shapes: Dict[TensorMap, Optional[Tuple[int, ...]]] = {}
    for tm in tensor_maps_out:
        if any([tm in out for out in u_connect.values()]) or tm.axes() == 1:
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
            raise ValueError(f'To use {BottleneckType.NoBottleNeck}, all output TensorMaps must be u-connected to.')
        bottleneck = UConnectBottleNeck(u_connect)
    else:
        raise NotImplementedError(f'Unknown BottleneckType {bottleneck_type}.')

    conv_x, conv_y, conv_z = conv_x[num_res:], conv_y[num_res:], conv_z[num_res:]
    decoders: Dict[TensorMap, Layer] = {}
    for tm in tensor_maps_out:
        if tm.axes() > 1:
            decoders[tm] = ConvDecoder(
                tensor_map_out=tm,
                filters_per_dense_block=dense_blocks[::-1],
                conv_layer_type=conv_type,
                conv_x=conv_x,
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Training ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def train_model_from_generators(
    model: Model,
    generate_train: Iterable,
    generate_valid: Iterable,
    training_steps: int,
    validation_steps: int,
    batch_size: int,
    epochs: int,
    patience: int,
    output_folder: str,
    run_id: str,
    inspect_model: bool,
    inspect_show_labels: bool,
    return_history: bool = False,
    plot: bool = True,
) -> Union[Model, Tuple[Model, History]]:
    """Train a model from tensor generators for validation and training data.

    Training data lives on disk, it will be loaded by generator functions.
    Plots the metric history after training. Creates a directory to save weights, if necessary.
    Measures runtime and plots architecture diagram if inspect_model is True.

    :param model: The model to optimize
    :param generate_train: Generator function that yields mini-batches of training data.
    :param generate_valid: Generator function that yields mini-batches of validation data.
    :param training_steps: Number of mini-batches in each so-called epoch
    :param validation_steps: Number of validation mini-batches to examine after each epoch.
    :param batch_size: Number of training examples in each mini-batch
    :param epochs: Maximum number of epochs to run regardless of Early Stopping
    :param patience: Number of epochs to wait before reducing learning rate.
    :param output_folder: Directory where output file will be stored
    :param run_id: User-chosen string identifying this run
    :param inspect_model: If True, measure training and inference runtime of the model and generate architecture plot.
    :param inspect_show_labels: If True, show labels on the architecture plot.
    :param return_history: If true return history from training and don't plot the training history
    :return: The optimized model.
    """
    model_file = os.path.join(output_folder, run_id, run_id + MODEL_EXT)
    if not os.path.exists(os.path.dirname(model_file)):
        os.makedirs(os.path.dirname(model_file))

    if inspect_model:
        image_p = os.path.join(output_folder, run_id, 'architecture_graph_' + run_id + IMAGE_EXT)
        _inspect_model(model, generate_train, generate_valid, batch_size, training_steps, inspect_show_labels, image_p)

    history = model.fit(
        generate_train, steps_per_epoch=training_steps, epochs=epochs, verbose=1,
        validation_steps=validation_steps, validation_data=generate_valid,
        callbacks=_get_callbacks(patience, model_file),
    )
    generate_train.kill_workers()
    generate_valid.kill_workers()

    logging.info('Model weights saved at: %s' % model_file)
    model = load_model(model_file, custom_objects=_get_custom_objects(generate_train.output_maps))
    if plot:
        plot_metric_history(history, training_steps, run_id, os.path.dirname(model_file))
    if return_history:
        return model, history
    return model


def _get_callbacks(
    patience: int, model_file: str,
) -> List[Callback]:
    callbacks = [
        ModelCheckpoint(filepath=model_file, verbose=1, save_best_only=True),
        EarlyStopping(monitor='val_loss', patience=patience * 3, verbose=1),
        ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=patience, verbose=1),
    ]
    return callbacks


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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Inspections ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def _inspect_model(
    model: Model,
    generate_train: Iterable,
    generate_valid: Iterable,
    batch_size: int,
    training_steps: int,
    inspect_show_labels: bool,
    image_path: str,
) -> Model:
    """Collect statistics on model inference and training times.

    Arguments:
        model: the model to inspect
        generate_train: training data generator function
        generate_valid: Validation data generator function
        batch_size: size of the mini-batches
        training_steps: number of optimization steps to take
        inspect_show_labels: if True, show layer labels on the architecture diagram
        image_path: file path of the architecture diagram

    Returns:
        The slightly optimized keras model
    """
    if image_path:
        _plot_dot_model_in_color(model_to_dot(model, show_shapes=inspect_show_labels, expand_nested=True), image_path, inspect_show_labels)
    t0 = time.time()
    _ = model.fit(generate_train, steps_per_epoch=training_steps, validation_steps=1, validation_data=generate_valid)
    t1 = time.time()
    n = batch_size * training_steps
    train_speed = (t1 - t0) / n
    logging.info(f'Spent:{(t1 - t0):0.2f} seconds training, Samples trained on:{n} Per sample training speed:{train_speed:0.3f} seconds.')
    t0 = time.time()
    _ = model.predict(generate_valid, steps=training_steps, verbose=1)
    t1 = time.time()
    inference_speed = (t1 - t0) / n
    logging.info(f'Spent:{(t1 - t0):0.2f} seconds predicting, Samples inferred:{n} Per sample inference speed:{inference_speed:0.4f} seconds.')
    return model


def _plot_dot_model_in_color(dot, image_path, inspect_show_labels):
    import pydot
    legend = {}
    for n in dot.get_nodes():
        if n.get_label():
            if 'Conv1' in n.get_label():
                legend['Conv1'] = "cyan"
                n.set_fillcolor("cyan")
            elif 'Conv2' in n.get_label():
                legend['Conv2'] = "deepskyblue1"
                n.set_fillcolor("deepskyblue1")
            elif 'Conv3' in n.get_label():
                legend['Conv3'] = "deepskyblue3"
                n.set_fillcolor("deepskyblue3")
            elif 'UpSampling' in n.get_label():
                legend['UpSampling'] = "darkslategray2"
                n.set_fillcolor("darkslategray2")
            elif 'Transpose' in n.get_label():
                legend['Transpose'] = "deepskyblue2"
                n.set_fillcolor("deepskyblue2")
            elif 'BatchNormalization' in n.get_label():
                legend['BatchNormalization'] = "goldenrod1"
                n.set_fillcolor("goldenrod1")
            elif 'output_' in n.get_label():
                n.set_fillcolor("darkolivegreen2")
                legend['Output'] = "darkolivegreen2"
            elif 'softmax' in n.get_label():
                n.set_fillcolor("chartreuse")
                legend['softmax'] = "chartreuse"
            elif 'MaxPooling' in n.get_label():
                legend['MaxPooling'] = "aquamarine"
                n.set_fillcolor("aquamarine")
            elif 'Dense' in n.get_label():
                legend['Dense'] = "gold"
                n.set_fillcolor("gold")
            elif 'Reshape' in n.get_label():
                legend['Reshape'] = "coral"
                n.set_fillcolor("coral")
            elif 'Input' in n.get_label():
                legend['Input'] = "darkolivegreen1"
                n.set_fillcolor("darkolivegreen1")
            elif 'Activation' in n.get_label():
                legend['Activation'] = "yellow"
                n.set_fillcolor("yellow")
        n.set_style("filled")
        if not inspect_show_labels:
            n.set_label('\n')

    for label in legend:
        legend_node = pydot.Node('legend' + label, label=label, shape="box", fillcolor=legend[label])
        dot.add_node(legend_node)

    logging.info('Saving architecture diagram to:{}'.format(image_path))
    dot.write_png(image_path)


def saliency_map(input_tensor: np.ndarray, model: Model, output_layer_name: str, output_index: int) -> np.ndarray:
    """Compute saliency maps of the given model (presumably already trained) on a batch of inputs with respect to the desired output layer and index.

    For example, with a trinary classification layer called quality and classes good, medium, and bad output layer name would be "quality_output"
    and output_index would be 0 to get gradients w.r.t. good, 1 to get gradients w.r.t. medium, and 2 for gradients w.r.t. bad.

    :param input_tensor: A batch of input tensors
    :param model: A trained model expecting those inputs
    :param output_layer_name: The name of the output layer that the derivative will be taken with respect to
    :param output_index: The index within the output layer that the derivative will be taken with respect to

    :return: Array of the gradients same shape as input_tensor
    """
    get_gradients = _gradients_from_output(model, output_layer_name, output_index)
    activation, gradients = get_gradients([input_tensor])
    return gradients


def _gradients_from_output(model, output_layer, output_index):
    K.set_learning_phase(1)
    input_tensor = model.input
    x = model.get_layer(output_layer).output[:, output_index]
    grads = K.gradients(x, input_tensor)[0]
    grads /= (K.sqrt(K.mean(K.square(grads))) + 1e-6)  # normalization trick: we normalize the gradient
    iterate = K.function([input_tensor], [x, grads])
    return iterate


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
    else:
        tm_char = TensorMap(
            f'{language_layer}{LANGUAGE_MODEL_SUFFIX}', Interpretation.LANGUAGE, shape=(len(ECG_CHAR_2_IDX),), channel_map=ECG_CHAR_2_IDX,
            cacheable=False,
        )
        tm_burn_in = TensorMap(
            language_layer, Interpretation.LANGUAGE, shape=(burn_in, len(ECG_CHAR_2_IDX)), path_prefix=language_prefix,
            dependent_map=tm_char, cacheable=False,
        )

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
