# models.py

# Imports
import os
import h5py
import time
import logging
import operator
import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable, Callable

# Keras imports
from keras import layers
import keras.backend as K
from keras.optimizers import Adam
from keras.models import Model, load_model
from keras.utils.vis_utils import model_to_dot
from keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau
from keras.layers import SpatialDropout1D, SpatialDropout2D, SpatialDropout3D, add, concatenate
from keras.layers import Input, Dense, Dropout, BatchNormalization, Activation, Flatten, LSTM, RepeatVector
from keras.layers.convolutional import Conv1D, Conv2D, Conv3D, UpSampling1D, UpSampling2D, UpSampling3D, MaxPooling1D
from keras.layers.convolutional import MaxPooling2D, MaxPooling3D, AveragePooling1D, AveragePooling2D, AveragePooling3D

from ml4cvd.TensorMap import TensorMap
from ml4cvd.metrics import get_metric_dict
from ml4cvd.plots import plot_metric_history
from ml4cvd.defines import JOIN_CHAR, IMAGE_EXT, TENSOR_EXT, ECG_CHAR_2_IDX

CHANNEL_AXIS = -1  # Set to 1 for Theano backend


def make_shallow_model(tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap],
                       learning_rate: float, model_file: str = None, model_layers: str = None) -> Model:
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
    if len(input_tensors) > 1:
        logging.warning('multi input tensors not fully supported')
    for it in input_tensors:
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


def make_waveform_model_unet(tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], learning_rate: float,
                             model_file: str = None, model_layers: str = None) -> Model:
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


def make_character_model_plus(tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], learning_rate: float,
                              base_model: Model, model_layers: str=None) -> Tuple[Model, Model]:
    """Make a ECG captioning model from an ECG embedding model

    The base_model must have an embedding layer, but besides that can have any number of other predicition TensorMaps.
    Input and output tensor maps are set from the command line.
    Model summary printed to output

    :param tensor_maps_in: List of input TensorMaps, only 1 input TensorMap is currently supported,
                            otherwise there are layer name collisions.
    :param tensor_maps_out: List of output TensorMaps
    :param learning_rate: Size of learning steps in SGD optimization
    :param base_model: The model the computes the ECG embedding
    :param model_layers: Optional HD5 model file whose weights will be loaded into this model when layer names match.
    :return: a tuple of the compiled keras model and the character emitting sub-model
    """
    char_maps_in, char_maps_out = _get_tensor_maps_for_characters(tensor_maps_in, base_model)
    tensor_maps_in.extend(char_maps_in)
    tensor_maps_out.extend(char_maps_out)
    char_model = make_character_model(tensor_maps_in, tensor_maps_out, learning_rate)
    losses = []
    my_metrics = {}
    loss_weights = []
    output_layers = []
    for tm in tensor_maps_out:
        losses.append(tm.loss)
        loss_weights.append(tm.loss_weight)
        my_metrics[tm.output_name()] = tm.metrics
        if tm.name == 'ecg_rest_next_char':
            output_layers.append(char_model.get_layer(tm.output_name()))
        else:
            output_layers.append(base_model.get_layer(tm.output_name()))

    m = Model(inputs=base_model.inputs+char_model.inputs, outputs=base_model.outputs+char_model.outputs)
    m.summary()
    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0)
    m.compile(optimizer=opt, loss=losses, loss_weights=loss_weights, metrics=my_metrics)

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        logging.info('Loaded model weights from:{}'.format(model_layers))

    return m, char_model


def make_character_model(tensor_maps_in: List[TensorMap], tensor_maps_out: List[TensorMap], learning_rate: float,
                         model_file: str = None, model_layers: str = None) -> Model:
    """Make a ECG captioning model

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

    input_layers = []
    for it in tensor_maps_in:
        if it.is_hidden_layer():
            embed_in = Input(shape=it.shape, name=it.input_name())
            input_layers.append(embed_in)
        elif it.is_ecg_text():
            burn_in = Input(shape=it.shape, name=it.input_name())
            input_layers.append(burn_in)
            repeater = RepeatVector(it.shape[0])
        else:
            logging.warning(f"character model cant handle {it.name} from group:{it.group}")

    logging.info(f"inputs: {[il.name for il in input_layers]}")
    wave_embeds = repeater(embed_in)
    lstm_in = concatenate([burn_in, wave_embeds], name='concat_embed_and_text')
    lstm_out = LSTM(128)(lstm_in)

    output_layers = []
    for ot in tensor_maps_out:
        if ot.name == 'ecg_rest_next_char':
            output_layers.append(Dense(ot.shape[-1], activation=ot.activation, name=ot.output_name())(lstm_out))

    m = Model(inputs=input_layers, outputs=output_layers)
    m.summary()
    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0)
    m.compile(optimizer=opt, loss='categorical_crossentropy')

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        logging.info('Loaded model weights from:{}'.format(model_layers))

    return m


def make_hidden_layer_model_from_file(parent_file: str, tensor_maps_in: List[TensorMap], output_layer_name: str, tensor_maps_out: List[TensorMap]):
    parent_model = load_model(parent_file, custom_objects=get_metric_dict(tensor_maps_out))
    return make_hidden_layer_model(parent_model, tensor_maps_in, output_layer_name)


def make_hidden_layer_model(parent_model: Model, tensor_maps_in: List[TensorMap], output_layer_name: str):
    parent_inputs = [parent_model.get_layer(tm.input_name()).input for tm in tensor_maps_in]
    dummy_input = {tm.input_name(): np.zeros((1,) + parent_model.get_layer(tm.input_name()).input_shape[1:]) for tm in tensor_maps_in}
    intermediate_layer_model = Model(inputs=parent_inputs, outputs=parent_model.get_layer(output_layer_name).output)
    # If we do not predict here then the graph is disconnected, I do not know why?!
    intermediate_layer_model.predict(dummy_input)
    return intermediate_layer_model


def make_multimodal_to_multilabel_model(model_file: str,
                                        model_layers: str,
                                        model_freeze: str,
                                        tensor_maps_in: List[TensorMap],
                                        tensor_maps_out: List[TensorMap],
                                        activation: str,
                                        dense_layers: List[int],
                                        dropout: float,
                                        mlp_concat: bool,
                                        conv_layers: List[int],
                                        max_pools: List[int],
                                        res_layers: List[int],
                                        dense_blocks: List[int],
                                        block_size: List[int],
                                        conv_bn: bool,
                                        conv_x: int,
                                        conv_y: int,
                                        conv_z: int,
                                        conv_dropout: float,
                                        conv_width: int,
                                        u_connect: bool,
                                        pool_x: int,
                                        pool_y: int,
                                        pool_z: int,
                                        padding: str,
                                        learning_rate: float,
                                        ) -> Model:
    """Make multi-task, multi-modal feed forward neural network for all kinds of prediction

	This model factory can be used to make networks for classification, regression, and segmentation
	The tasks attempted are given by the output TensorMaps.
	The modalities and the first layers in the architecture are determined by the input TensorMaps.

	Hyperparameters are exposed to the command line.
	Model summary printed to output

    :param model_file: HD5 model file to load and return.
    :param model_layers: HD5 model file whose weights will be loaded into this model when layer names match.
    :param model_freeze: HD5 model file whose weights will be loaded and frozen into this model when layer names match.
    :param tensor_maps_in: List of input TensorMaps
    :param tensor_maps_out: List of output TensorMaps
    :param activation: Activation function as a string (e.g. 'relu', 'linear, or 'softmax)
    :param dense_layers: List of number of filters in each dense layer.
    :param dropout: Dropout rate in dense layers
    :param mlp_concat: If True, conncatenate unstructured inputs to each deeper dense layer
    :param conv_layers: List of number of filters in each convolutional layer
    :param max_pools: List of maxpool sizes in X and Y dimensions after convolutional layers
    :param res_layers: List of convolutional layers with residual connections
    :param dense_blocks: List of number of filters in densenet modules for densenet convolutional models
    :param block_size: Number of layers within each Densenet module for densenet convolutional models
    :param conv_bn: if True, Batch normalize convolutional layers
    :param conv_x: Size of X dimension for 2D and 3D convolutional kernels
    :param conv_y: Size of Y dimension for 2D and 3D convolutional kernels
    :param conv_z: Size of Z dimension for 3D convolutional kernels
    :param conv_dropout: Dropout rate in convolutional layers
    :param conv_width: Size of convolutional kernel for 1D models.
    :param u_connect: Include U connections between early and late convolutional layers.
    :param pool_x: Pooling in the X dimension for Convolutional models.
    :param pool_y: Pooling in the Y dimension for Convolutional models.
    :param pool_z: Pooling in the Z dimension for 3D Convolutional models.
    :param padding: Padding string can be 'valid' or 'same'. UNets and residual nets require 'same'.
    :param learning_rate:
    :return: a compiled keras model
	"""
    if model_file is not None:
        logging.info("Attempting to load model file from: {}".format(model_file))
        m = load_model(model_file, custom_objects=get_metric_dict(tensor_maps_out))
        m.summary()
        logging.info("Loaded model file from: {}".format(model_file))
        return m

    input_tensors = [Input(shape=tm.shape, name=tm.input_name()) for tm in tensor_maps_in]
    input_multimodal = []
    upsamplers = []
    channel_axis = -1

    for j, tm in enumerate(tensor_maps_in):
        if len(tm.shape) == 4:
            last_convolution3d = _conv_block3d(input_tensors[j], upsamplers, conv_layers, max_pools, res_layers,
                                               activation, conv_bn, (conv_x, conv_y, conv_z), conv_dropout, padding)
            last_convolution3d = _dense_block3d(last_convolution3d, upsamplers, dense_blocks, block_size, activation,
                                                conv_bn, (conv_x, conv_y, conv_z), (pool_x, pool_y, pool_z), conv_dropout, padding)
            input_multimodal.append(Flatten()(last_convolution3d))
        elif len(tm.shape) == 3:
            last_convolution2d = _conv_block2d(input_tensors[j], upsamplers, conv_layers, max_pools, res_layers,
                                               activation, conv_bn, (conv_x, conv_y), conv_dropout, padding)
            last_convolution2d = _dense_block2d(last_convolution2d, upsamplers, dense_blocks, block_size, activation,
                                                conv_bn, (conv_x, conv_y), (pool_x, pool_y), conv_dropout, padding)
            input_multimodal.append(Flatten()(last_convolution2d))
        elif len(tm.shape) == 2:
            last_convolution1d = _conv_block1d(input_tensors[j], upsamplers,  conv_layers, max_pools, res_layers,
                                               activation, conv_bn, conv_width, conv_dropout, padding)
            last_convolution1d = _dense_block1d(last_convolution1d, upsamplers,  dense_blocks, block_size, activation,
                                                conv_bn, conv_width, conv_dropout, pool_x, padding)
            input_multimodal.append(Flatten()(last_convolution1d))
        else:
            mlp_input = input_tensors[j]
            mlp = Dense(units=tm.annotation_units, activation=activation)(mlp_input)
            input_multimodal.append(mlp)

    if len(input_multimodal) > 1:
        multimodal_activation = concatenate(input_multimodal, axis=channel_axis)
    elif len(input_multimodal) == 1:
        multimodal_activation = input_multimodal[0]
    else:
        raise ValueError('No input activations.')

    for i, hidden_units in enumerate(dense_layers):
        if conv_bn:
            multimodal_activation = BatchNormalization()(multimodal_activation)
        if i == len(dense_layers)-1:
            multimodal_activation = Dense(units=hidden_units, activation=activation, name='embed')(multimodal_activation)
        else:
            multimodal_activation = Dense(units=hidden_units, activation=activation)(multimodal_activation)
        if dropout > 0:
            multimodal_activation = Dropout(dropout)(multimodal_activation)
        if mlp_concat:
            multimodal_activation = concatenate([multimodal_activation, mlp_input], axis=channel_axis)

    losses = []
    my_metrics = {}
    loss_weights = []
    output_predictions = {}
    output_tensor_maps_to_process = tensor_maps_out.copy()
    while len(output_tensor_maps_to_process) > 0:
        tm = output_tensor_maps_to_process.pop(0)

        if not tm.parents is None and any(not p in output_predictions for p in tm.parents):
            output_tensor_maps_to_process.append(tm)
            continue

        losses.append(tm.loss)
        loss_weights.append(tm.loss_weight)
        my_metrics[tm.output_name()] = tm.metrics

        if len(tm.shape) == 4:
            for x, up_conv, upsampler in reversed(upsamplers):
                if u_connect:
                    last_convolution3d = concatenate([up_conv(upsampler(last_convolution3d)), x])
                else:
                    last_convolution3d = upsampler(last_convolution3d)
            conv_label = Conv3D(tm.shape[channel_axis], (1, 1, 1), activation="linear")(last_convolution3d)
            output_predictions[tm.output_name()] = Activation(tm.activation, name=tm.output_name())(conv_label)
        elif len(tm.shape) == 3:
            for x, up_conv, upsampler in reversed(upsamplers):
                if u_connect:
                    last_convolution2d = concatenate([up_conv(upsampler(last_convolution2d)), x])
                else:
                    last_convolution2d = upsampler(last_convolution2d)
            conv_label = Conv2D(tm.shape[channel_axis], (1, 1), activation="linear")(last_convolution2d)
            output_predictions[tm.output_name()] = Activation(tm.activation, name=tm.output_name())(conv_label)
        elif len(tm.shape) == 2:
            for x, up_conv, upsampler in reversed(upsamplers):
                if u_connect:
                    last_convolution1d = concatenate([up_conv(upsampler(last_convolution1d)), x])
                else:
                    last_convolution1d = upsampler(last_convolution1d)
            conv_label = Conv1D(tm.shape[channel_axis], 1, activation="linear")(last_convolution1d)
            output_predictions[tm.output_name()] = Activation(tm.activation, name=tm.output_name())(conv_label)
            sx = _conv_block1d(conv_label, [],  conv_layers, max_pools, res_layers, activation, conv_bn, conv_width, conv_dropout, padding)
            flat_activation = Flatten()(sx)
            for hidden_units in dense_layers:
                flat_activation = Dense(units=hidden_units, activation=activation)(flat_activation)
                if dropout > 0:
                    flat_activation = Dropout(dropout)(flat_activation)
            multimodal_activation = concatenate([multimodal_activation, flat_activation])
        elif tm.parents is not None:
            if False and len(K.int_shape(output_predictions[tm.parents[0]])) > 1:  # TODO: fix False
                output_predictions[tm.output_name()] = Dense(units=tm.shape[0], activation=tm.activation, name=tm.output_name())(multimodal_activation)
            else:
                parented_activation = concatenate([multimodal_activation] + [output_predictions[p] for p in tm.parents])
                parented_activation = Dense(units=tm.annotation_units, activation=activation)(parented_activation)
                parented_activation = concatenate([parented_activation] + [output_predictions[p] for p in tm.parents])
                output_predictions[tm.output_name()] = Dense(units=tm.shape[0], activation=tm.activation, name=tm.output_name())(parented_activation)
        elif tm.is_categorical_any():
            output_predictions[tm.output_name()] = Dense(units=tm.shape[0], activation='softmax', name=tm.output_name())(multimodal_activation)
        else:
            output_predictions[tm.output_name()] = Dense(units=1, activation=tm.activation, name=tm.output_name())(multimodal_activation)

    m = Model(inputs=input_tensors, outputs=list(output_predictions.values()))
    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0)
    m.summary()

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        logging.info('Loaded model weights from:{}'.format(model_layers))

    if model_freeze is not None:
        frozen = 0
        m.load_weights(model_freeze, by_name=True)
        m_freeze = load_model(model_freeze, custom_objects=get_metric_dict(tensor_maps_out))
        frozen_layers = [layer.name for layer in m_freeze.layers]
        for l in m.layers:
            if l.name in frozen_layers:
                l.trainable = False
                frozen += 1
        logging.info('Loaded and froze:{} layers from:{}'.format(frozen, model_freeze))

    m.compile(optimizer=opt, loss=losses, loss_weights=loss_weights, metrics=my_metrics)
    return m


def make_translation_model(model_file: str,
                                        model_layers: str,
                                        model_freeze: str,
                                        tensor_maps_in: List[TensorMap],
                                        tensor_maps_out: List[TensorMap],
                                        activation: str,
                                        dense_layers: List[int],
                                        dropout: float,
                                        mlp_concat: bool,
                                        conv_layers: List[int],
                                        max_pools: List[int],
                                        res_layers: List[int],
                                        dense_blocks: List[int],
                                        block_size: List[int],
                                        conv_bn: bool,
                                        conv_x: int,
                                        conv_y: int,
                                        conv_z: int,
                                        conv_dropout: float,
                                        conv_width: int,
                                        u_connect: bool,
                                        pool_x: int,
                                        pool_y: int,
                                        pool_z: int,
                                        padding: str,
                                        learning_rate: float,
                                        ) -> Model:
    """Make multi-task, multi-modal feed forward neural network for all kinds of prediction

	This model factory can be used to make networks for classification, regression, and segmentation
	The tasks attempted are given by the output TensorMaps.
	The modalities and the first layers in the architecture are determined by the input TensorMaps.

	Hyperparameters are exposed to the command line.
	Model summary printed to output

    :param model_file: HD5 model file to load and return.
    :param model_layers: HD5 model file whose weights will be loaded into this model when layer names match.
    :param model_freeze: HD5 model file whose weights will be loaded and frozen into this model when layer names match.
    :param tensor_maps_in: List of input TensorMaps
    :param tensor_maps_out: List of output TensorMaps
    :param activation: Activation function as a string (e.g. 'relu', 'linear, or 'softmax)
    :param dense_layers: List of number of filters in each dense layer.
    :param dropout: Dropout rate in dense layers
    :param mlp_concat: If True, conncatenate unstructured inputs to each deeper dense layer
    :param conv_layers: List of number of filters in each convolutional layer
    :param max_pools: List of maxpool sizes in X and Y dimensions after convolutional layers
    :param res_layers: List of convolutional layers with residual connections
    :param dense_blocks: List of number of filters in densenet modules for densenet convolutional models
    :param block_size: Number of layers within each Densenet module for densenet convolutional models
    :param conv_bn: if True, Batch normalize convolutional layers
    :param conv_x: Size of X dimension for 2D and 3D convolutional kernels
    :param conv_y: Size of Y dimension for 2D and 3D convolutional kernels
    :param conv_z: Size of Z dimension for 3D convolutional kernels
    :param conv_dropout: Dropout rate in convolutional layers
    :param conv_width: Size of convolutional kernel for 1D models.
    :param u_connect: Include U connections between early and late convolutional layers.
    :param pool_x: Pooling in the X dimension for Convolutional models.
    :param pool_y: Pooling in the Y dimension for Convolutional models.
    :param pool_z: Pooling in the Z dimension for 3D Convolutional models.
    :param padding: Padding string can be 'valid' or 'same'. UNets and residual nets require 'same'.
    :param learning_rate:
    :return: a compiled keras model
	"""
    if model_file is not None:
        logging.info("Attempting to load model file from: {}".format(model_file))
        m = load_model(model_file, custom_objects=get_metric_dict(tensor_maps_out))
        m.summary()
        logging.info("Loaded model file from: {}".format(model_file))
        return m

    input_tensors = [Input(shape=tm.shape, name=tm.input_name()) for tm in tensor_maps_in]
    input_multimodal = []
    channel_axis = -1
    layers = {}

    for j, tm in enumerate(tensor_maps_in):
        if len(tm.shape) == 4:
            last_convolution3d = _conv_block3d(input_tensors[j], [], conv_layers, max_pools, res_layers, activation, conv_bn, (conv_x, conv_y, conv_z), conv_dropout, padding)
            input_multimodal.append(Flatten()(last_convolution3d))
        elif len(tm.shape) == 3:
            last_convolution2d = _conv_block2d_new(input_tensors[j], layers, conv_layers, max_pools, res_layers, activation, conv_bn, (conv_x, conv_y), conv_dropout, padding)
            last_convolution2d = _dense_block2d_new(last_convolution2d, layers, dense_blocks, block_size, activation, conv_bn, (conv_x, conv_y), (pool_x, pool_y), conv_dropout, padding)
        elif len(tm.shape) == 2:
            last_convolution1d = _conv_block1d(input_tensors[j], [],  conv_layers, max_pools, res_layers, activation, conv_bn, conv_width, conv_dropout, padding)
            input_multimodal.append(Flatten()(last_convolution1d))
        else:
            mlp_input = input_tensors[j]
            mlp = Dense(units=tm.annotation_units, activation=activation)(mlp_input)
            input_multimodal.append(mlp)
            for i, hidden_units in enumerate(dense_layers):
                if conv_bn:
                    multimodal_activation = BatchNormalization()(mlp)
                if i == len(dense_layers)-1:
                    multimodal_activation = Dense(units=hidden_units, activation=activation, name='embed')(multimodal_activation)
                else:
                    multimodal_activation = Dense(units=hidden_units, activation=activation)(multimodal_activation)
                if dropout > 0:
                    multimodal_activation = Dropout(dropout)(multimodal_activation)
                if mlp_concat:
                    multimodal_activation = concatenate([multimodal_activation, mlp_input], axis=channel_axis)

    losses = []
    my_metrics = {}
    loss_weights = []
    output_predictions = {}
    output_tensor_maps_to_process = tensor_maps_out.copy()
    print(layers)
    while len(output_tensor_maps_to_process) > 0:
        tm = output_tensor_maps_to_process.pop(0)

        if not tm.parents is None and any(not p in output_predictions for p in tm.parents):
            output_tensor_maps_to_process.append(tm)
            continue

        losses.append(tm.loss)
        loss_weights.append(tm.loss_weight)
        my_metrics[tm.output_name()] = tm.metrics

        if len(tm.shape) == 4:
            conv_label = Conv3D(tm.shape[channel_axis], (1, 1, 1), activation="linear")(last_convolution3d)
            output_predictions[tm.output_name()] = Activation(tm.activation, name=tm.output_name())(conv_label)
        elif len(tm.shape) == 3:
            print('got rev la to be layers:', _get_layer_kind_sorted(layers, 'Pooling2D'))
            print('got rev la to be layers:', _get_layer_kind_sorted(layers, 'Conv2D'))
            all_filters = conv_layers + dense_blocks
            for i, name in enumerate(reversed(_get_layer_kind_sorted(layers, 'Pooling2D'))):
                print("!!!!!!!!!!!!!!!!!!! for named layer", name)
                print(f"Conv2D{JOIN_CHAR}{_get_layer_index_offset_str(name, -1)}")
                early_conv = layers[f"Conv2D{JOIN_CHAR}{_get_layer_index_offset_str(name, -1)}"]
                if u_connect:
                    last_convolution2d = UpSampling2D((pool_x, pool_y))(last_convolution2d)
                    last_convolution2d = Conv2D(filters=all_filters[-(1+i)], kernel_size=(conv_x, conv_y), activation=activation, padding=padding)(last_convolution2d)
                    last_convolution2d = concatenate([last_convolution2d, early_conv])
                else:
                    last_convolution2d = UpSampling2D((pool_x, pool_y))(last_convolution2d)
            conv_label = Conv2D(tm.shape[channel_axis], (1, 1), activation="linear")(last_convolution2d)
            output_predictions[tm.output_name()] = Activation(tm.activation, name=tm.output_name())(conv_label)
        elif len(tm.shape) == 2:
            conv_label = Conv1D(tm.shape[channel_axis], 1, activation="linear")(last_convolution1d)
            output_predictions[tm.output_name()] = Activation(tm.activation, name=tm.output_name())(conv_label)
            sx = _conv_block1d(conv_label, [],  conv_layers, max_pools, res_layers, activation, conv_bn, conv_width, conv_dropout, padding)
            flat_activation = Flatten()(sx)
            for hidden_units in dense_layers:
                flat_activation = Dense(units=hidden_units, activation=activation)(flat_activation)
                if dropout > 0:
                    flat_activation = Dropout(dropout)(flat_activation)
            multimodal_activation = concatenate([multimodal_activation, flat_activation])
        elif tm.parents is not None:
            if len(K.int_shape(output_predictions[tm.parents[0]])) > 1:
                output_predictions[tm.output_name()] = Dense(units=tm.shape[0], activation=tm.activation, name=tm.output_name())(multimodal_activation)
            else:
                parented_activation = concatenate([multimodal_activation] + [output_predictions[p] for p in tm.parents])
                parented_activation = Dense(units=tm.annotation_units, activation=activation)(parented_activation)
                parented_activation = concatenate([parented_activation] + [output_predictions[p] for p in tm.parents])
                output_predictions[tm.output_name()] = Dense(units=tm.shape[0], activation=tm.activation, name=tm.output_name())(parented_activation)
        elif tm.is_categorical_any():
            output_predictions[tm.output_name()] = Dense(units=tm.shape[0], activation='softmax', name=tm.output_name())(multimodal_activation)
        else:
            output_predictions[tm.output_name()] = Dense(units=1, activation=tm.activation, name=tm.output_name())(multimodal_activation)

    m = Model(inputs=input_tensors, outputs=list(output_predictions.values()))
    opt = Adam(lr=learning_rate, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.0)
    m.summary()

    if model_layers is not None:
        m.load_weights(model_layers, by_name=True)
        logging.info('Loaded model weights from:{}'.format(model_layers))

    if model_freeze is not None:
        frozen = 0
        m.load_weights(model_freeze, by_name=True)
        m_freeze = load_model(model_freeze, custom_objects=get_metric_dict(tensor_maps_out))
        frozen_layers = [layer.name for layer in m_freeze.layers]
        for l in m.layers:
            if l.name in frozen_layers:
                l.trainable = False
                frozen += 1
        logging.info('Loaded and froze:{} layers from:{}'.format(frozen, model_freeze))

    m.compile(optimizer=opt, loss=losses, loss_weights=loss_weights, metrics=my_metrics)
    return m



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Training ~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def train_model_from_generators(model: Model,
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
                                inspect_show_labels: bool) -> Model:
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
    :return: The optimized model.
	"""
    model_file = os.path.join(output_folder, run_id, run_id + TENSOR_EXT)
    if not os.path.exists(os.path.dirname(model_file)):
        os.makedirs(os.path.dirname(model_file))

    if inspect_model:
        image_p = os.path.join(output_folder, run_id, 'architecture_graph_' + run_id + IMAGE_EXT)
        _inspect_model(model, generate_train, generate_valid, batch_size, training_steps, inspect_show_labels, image_p)

    history = model.fit_generator(generate_train, steps_per_epoch=training_steps, epochs=epochs, verbose=1,
                                  validation_steps=validation_steps, validation_data=generate_valid,
                                  callbacks=_get_callbacks(patience, model_file))

    plot_metric_history(history, run_id, os.path.dirname(model_file))
    logging.info('Model weights saved at: %s' % model_file)

    return model


def _get_callbacks(patience: int, model_file: str) -> List[Callable]:
    callbacks = [
        ModelCheckpoint(filepath=model_file, verbose=1, save_best_only=True),
        EarlyStopping(monitor='val_loss', patience=patience * 3, verbose=1),
        ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=patience, verbose=1)
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
def _conv_block3d(x: K.placeholder,
                  upsamplers: List[Tuple[K.placeholder, Conv3D, K.placeholder]],
                  conv_layers: List[int],
                  max_pools: List[int],
                  res_layers: List[int],
                  activation: str,
                  conv_bn: bool,
                  kernel: Tuple[int, int, int],
                  conv_dropout: float,
                  padding: str):
    max_pool_diff = len(conv_layers) - len(max_pools)
    residual_diff = len(conv_layers) - len(res_layers)
    for i, c in enumerate(conv_layers):
        residual2d = x
        if conv_bn and i > 0:
            x = BatchNormalization(axis=CHANNEL_AXIS)(x)
            x = Activation(activation)(x)
            x = Conv3D(filters=c, kernel_size=kernel, activation='linear', padding=padding)(x)
        else:
            x = Conv3D(filters=c, kernel_size=kernel, activation=activation, padding=padding)(x)
        if conv_dropout > 0:
            x = SpatialDropout3D(conv_dropout)(x)
        if i >= max_pool_diff:
            pool_size = (
                max_pools[i - max_pool_diff], max_pools[i - max_pool_diff], max_pools[i - max_pool_diff])
            up_conv = Conv3D(filters=c, kernel_size=kernel, activation=activation, padding=padding)
            upsamplers.append((residual2d, up_conv, UpSampling3D(pool_size)))
            x = MaxPooling3D(pool_size=pool_size)(x)
            if i >= residual_diff:
                residual2d = MaxPooling3D(pool_size=pool_size)(residual2d)
        if i >= residual_diff:
            if K.int_shape(x)[CHANNEL_AXIS] == K.int_shape(residual2d)[CHANNEL_AXIS]:
                x = add([x, residual2d])
            else:
                residual2d = Conv3D(filters=K.int_shape(x)[CHANNEL_AXIS], kernel_size=(1, 1, 1))(residual2d)
                x = add([x, residual2d])
    return x


def _dense_block3d(x: K.placeholder,
                   upsamplers: List[Tuple[K.placeholder, Conv3D, K.placeholder]],
                   dense_blocks: List[int],
                   block_size: int,
                   activation: str,
                   conv_bn: bool,
                   kernel: Tuple[int, int, int],
                   pool_size: Tuple[int, int, int],
                   conv_dropout: float,
                   padding: str,
                   pool_mode: str='average'):
    for db_filters in dense_blocks:
        for i in range(block_size):
            residual3d = x
            if conv_bn and i > 0:
                x = BatchNormalization(axis=CHANNEL_AXIS)(x)
                x = Activation(activation)(x)
                x = Conv3D(filters=db_filters, kernel_size=kernel, activation='linear', padding=padding)(x)
            else:
                x = Conv3D(filters=db_filters, kernel_size=kernel, activation=activation, padding=padding)(x)
            if conv_dropout > 0:
                x = SpatialDropout3D(conv_dropout)(x)

            if i == 0:
                up_conv = Conv3D(filters=db_filters, kernel_size=kernel, activation=activation, padding=padding)
                upsamplers.append((residual3d, up_conv, UpSampling3D(pool_size)))
                if pool_mode == 'average':
                    x = AveragePooling3D(pool_size, strides=pool_size)(x)
                elif pool_mode == 'max':
                    x = MaxPooling3D(pool_size, strides=pool_size)(x)
                dense_connections = [x]
            else:
                dense_connections += [x]
                x = concatenate(dense_connections, axis=CHANNEL_AXIS)
    return x


def _conv_block2d(x: K.placeholder,
                  upsamplers: List[Tuple[K.placeholder, Conv2D, K.placeholder]],
                  conv_layers: List[int],
                  max_pools: List[int],
                  res_layers: List[int],
                  activation: str,
                  conv_bn: bool,
                  kernel: Tuple[int, int],
                  conv_dropout: float,
                  padding: str):
    max_pool_diff = len(conv_layers) - len(max_pools)
    residual_diff = len(conv_layers) - len(res_layers)
    for i, c in enumerate(conv_layers):
        residual2d = x
        if conv_bn and i > 0:
            x = BatchNormalization(axis=CHANNEL_AXIS)(x)
            x = Activation(activation)(x)
            x = Conv2D(filters=c, kernel_size=kernel, activation='linear', padding=padding)(x)
        else:
            x = Conv2D(filters=c, kernel_size=kernel, activation=activation, padding=padding)(x)
        if conv_dropout > 0:
            x = SpatialDropout2D(conv_dropout)(x)
        if i >= max_pool_diff:
            pool_size = (max_pools[i - max_pool_diff], max_pools[i - max_pool_diff])
            up_conv = Conv2D(filters=c, kernel_size=kernel, activation=activation, padding=padding)
            upsamplers.append((residual2d, up_conv, UpSampling2D(pool_size)))
            x = MaxPooling2D(pool_size=pool_size)(x)
            if i >= residual_diff:
                residual2d = MaxPooling2D(pool_size=pool_size)(residual2d)
        if i >= residual_diff:
            if K.int_shape(x)[CHANNEL_AXIS] == K.int_shape(residual2d)[CHANNEL_AXIS]:
                x = add([x, residual2d])
            else:
                residual2d = Conv2D(filters=K.int_shape(x)[CHANNEL_AXIS], kernel_size=(1, 1))(residual2d)
                x = add([x, residual2d])
    return x


def _conv_block2d_new(x: K.placeholder,
                      layers: Dict[str, K.placeholder],
                      conv_layers: List[int],
                      max_pools: List[int],
                      res_layers: List[int],
                      activation: str,
                      conv_bn: bool,
                      kernel: Tuple[int, int],
                      conv_dropout: float,
                      padding: str):
    max_pool_diff = len(conv_layers) - len(max_pools)
    residual_diff = len(conv_layers) - len(res_layers)
    for i, c in enumerate(conv_layers):
        residual2d = x
        if conv_bn and i > 0:
            x = layers[f"BatchNormalization_{str(len(layers))}"] = BatchNormalization(axis=CHANNEL_AXIS)(x)
            x = layers[f"Activation_{str(len(layers))}"] = Activation(activation)(x)
            x = layers[f"Conv2D_{str(len(layers))}"] = Conv2D(filters=c, kernel_size=kernel, activation='linear', padding=padding)(x)
        else:
            x = layers[f"Conv2D_{str(len(layers))}"] = Conv2D(filters=c, kernel_size=kernel, activation=activation, padding=padding)(x)
        if conv_dropout > 0:
            x = layers[f"SpatialDropout2D_{str(len(layers))}"] = SpatialDropout2D(conv_dropout)(x)
        if i >= max_pool_diff:
            pool_size = (max_pools[i - max_pool_diff], max_pools[i - max_pool_diff])
            x = layers[f"MaxPooling2D_{str(len(layers))}"] = MaxPooling2D(pool_size=pool_size)(x)
            if i >= residual_diff:
                residual2d = layers[f"MaxPooling2D_{str(len(layers))}"] = MaxPooling2D(pool_size=pool_size)(residual2d)
        if i >= residual_diff:
            if K.int_shape(x)[CHANNEL_AXIS] == K.int_shape(residual2d)[CHANNEL_AXIS]:
                x = layers[f"add_{str(len(layers))}"] = add([x, residual2d])
            else:
                residual2d = layers[f"Conv2D_{str(len(layers))}"] = Conv2D(filters=K.int_shape(x)[CHANNEL_AXIS], kernel_size=(1, 1))(residual2d)
                x = layers[f"add_{str(len(layers))}"] = add([x, residual2d])

    return _get_last_layer(layers)


def _get_last_layer(named_layers):
    max_index = -1
    max_layer = ''
    for k in named_layers:
        cur_index = int(k.split('_')[-1])
        if max_index < cur_index:
            max_index = cur_index
            max_layer = k
    return named_layers[max_layer]


def _get_last_layer_by_kind(named_layers, kind):
    max_index = -1
    for k in named_layers:
        if kind in k:
            max_index = max(max_index, int(k.split('_')[-1]))
    return named_layers[kind + JOIN_CHAR + str(max_index)]


def _get_layer_index_offset_str(named_layer, offset):
    return str(int(named_layer.split('_')[-1]) + offset)


def _get_layer_kind_sorted(named_layers, kind):
    return [k for k in sorted(list(named_layers.keys()), key=lambda x: int(x.split('_')[-1])) if kind in k]


def _dense_block2d(x: K.placeholder,
                   upsamplers: List[Tuple[K.placeholder, Conv2D, K.placeholder]],
                   dense_blocks: List[int],
                   block_size: int,
                   activation: str,
                   conv_bn: bool,
                   kernel: Tuple[int, int],
                   pool_size: Tuple[int, int],
                   conv_dropout: float,
                   padding: str,
                   pool_mode: str='average'):
    for db_filters in dense_blocks:
        for i in range(block_size):
            residual2d = x
            if conv_bn and i > 0:
                x = BatchNormalization(axis=CHANNEL_AXIS)(x)
                x = Activation(activation)(x)
                x = Conv2D(filters=db_filters, kernel_size=kernel, activation='linear', padding=padding)(x)
            else:
                x = Conv2D(filters=db_filters, kernel_size=kernel, activation=activation, padding=padding)(x)
            if conv_dropout > 0:
                x = SpatialDropout2D(conv_dropout)(x)

            if i == 0:
                up_conv = Conv2D(filters=db_filters, kernel_size=kernel, activation=activation, padding=padding)
                upsamplers.append((residual2d, up_conv, UpSampling2D(pool_size)))
                if pool_mode == 'average':
                    x = AveragePooling2D(pool_size, strides=pool_size)(x)
                elif pool_mode == 'max':
                    x = MaxPooling2D(pool_size, strides=pool_size)(x)
                dense_connections = [x]
            else:
                dense_connections += [x]
                x = concatenate(dense_connections, axis=CHANNEL_AXIS)
    return x


def _dense_block2d_new(x: K.placeholder,
                       layers: Dict[str, K.placeholder],
                       dense_blocks: List[int],
                       block_size: int,
                       activation: str,
                       conv_bn: bool,
                       kernel: Tuple[int, int],
                       pool_size: Tuple[int, int],
                       conv_dropout: float,
                       padding: str,
                       pool_mode='average'):
    for db_filters in dense_blocks:
        for i in range(block_size):
            if conv_bn and i > 0:
                x = layers[f"BatchNormalization_{str(len(layers))}"] = BatchNormalization(axis=CHANNEL_AXIS)(x)
                x = layers[f"Activation_{str(len(layers))}"] = Activation(activation)(x)
                x = layers[f"Conv2D_{str(len(layers))}"] = Conv2D(filters=db_filters, kernel_size=kernel, activation='linear', padding=padding)(x)
            else:
                x = layers[f"Conv2D_{str(len(layers))}"] = Conv2D(filters=db_filters, kernel_size=kernel, activation=activation, padding=padding)(x)
            if conv_dropout > 0:
                x = layers[f"SpatialDropout2D_{str(len(layers))}"] = SpatialDropout2D(conv_dropout)(x)

            if i == 0:
                if pool_mode == 'average':
                    x = layers[f"AveragePooling2D{JOIN_CHAR}{str(len(layers))}"] = AveragePooling2D(pool_size, strides=pool_size)(x)
                elif pool_mode == 'max':
                    x = layers[f"MaxPooling2D{JOIN_CHAR}{str(len(layers))}"] = MaxPooling2D(pool_size, strides=pool_size)(x)
                dense_connections = [x]
            else:
                dense_connections += [x]
                x = layers[f"concatenate{JOIN_CHAR}{str(len(layers))}"] = concatenate(dense_connections, axis=CHANNEL_AXIS)
    return _get_last_layer(layers)

def _conv_block1d(x: K.placeholder,
                  upsamplers: List[Tuple[K.placeholder, Conv1D, K.placeholder]],
                  conv_layers: List[int],
                  max_pools: List[int],
                  res_layers: List[int],
                  activation: str,
                  conv_bn: bool,
                  conv_width: int,
                  conv_dropout: float,
                  padding: str):
    residual_diff = len(conv_layers) - len(res_layers)
    max_pool_diff = len(conv_layers) - len(max_pools)
    for i, c in enumerate(conv_layers):
        residual1d = x
        if conv_bn and i > 0:
            x = BatchNormalization(axis=CHANNEL_AXIS)(x)
            x = Activation(activation)(x)
            x = Conv1D(filters=c, kernel_size=conv_width, activation='linear', padding=padding)(x)
        else:
            x = Conv1D(filters=c, kernel_size=conv_width, activation=activation, padding=padding)(x)
        if conv_dropout > 0:
            x = SpatialDropout1D(conv_dropout)(x)
        if i >= max_pool_diff:
            pool_size = max_pools[i - max_pool_diff]
            up_conv = Conv1D(filters=c, kernel_size=conv_width, activation=activation, padding=padding)
            upsamplers.append((residual1d, up_conv, UpSampling1D(pool_size)))
            x = MaxPooling1D(pool_size=pool_size)(x)
            if i >= residual_diff:
                residual1d = MaxPooling1D(pool_size=pool_size)(residual1d)
        if i >= residual_diff:
            if K.int_shape(x)[CHANNEL_AXIS] == K.int_shape(residual1d)[CHANNEL_AXIS]:
                x = add([x, residual1d])
            else:
                residual1d = Conv1D(filters=K.int_shape(x)[CHANNEL_AXIS], kernel_size=1)(residual1d)
                x = add([x, residual1d])
    return x


def _dense_block1d(x: K.placeholder,
                   upsamplers: List[Tuple[K.placeholder, Conv1D, K.placeholder]],
                   dense_blocks: List[int],
                   block_size: int,
                   activation: str,
                   conv_bn: bool,
                   conv_width: int,
                   conv_dropout: float,
                   pool_x: int,
                   padding: str,
                   pool_mode: str='max'):
    for db_filters in dense_blocks:
        for i in range(block_size):
            residual1d = x
            if conv_bn and i > 0:
                x = BatchNormalization(axis=CHANNEL_AXIS)(x)
                x = Activation(activation)(x)
                x = Conv1D(filters=db_filters, kernel_size=conv_width, activation='linear', padding=padding)(x)
            else:
                x = Conv1D(filters=db_filters, kernel_size=conv_width, activation=activation, padding=padding)(x)
            if conv_dropout > 0:
                x = SpatialDropout1D(conv_dropout)(x)

            if i == 0:
                up_conv = Conv1D(filters=db_filters, kernel_size=conv_width, activation=activation, padding=padding)
                upsamplers.append((residual1d, up_conv, UpSampling1D(pool_x)))
                if pool_mode == 'average':
                    x = AveragePooling1D(pool_x, strides=pool_x)(x)
                elif pool_mode == 'max':
                    x = MaxPooling1D(pool_x, strides=pool_x)(x)
                dense_connections = [x]
            else:
                dense_connections += [x]
                x = concatenate(dense_connections, axis=CHANNEL_AXIS)
    return x


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Inspections ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def _inspect_model(model: Model,
                  generate_train: Iterable,
                  generate_valid: Iterable,
                  batch_size: int,
                  training_steps: int,
                  inspect_show_labels: bool,
                  image_path: str) -> Model:
    """Collect statistics on model inference and training times.

	Arguments
	    model: the model to inspect
		generate_train: training data generator function
		generate_valid: Validation data generator function
		batch_size: size of the mini-batches
		training_steps: number of optimization steps to take
		inspect_show_labels: if True, show layer labels on the architecture diagram
		image_path: file path of the architecture diagram

	Returns
		The slightly optimized keras model
	"""
    if image_path:
        _plot_dot_model_in_color(model_to_dot(model, show_shapes=inspect_show_labels), image_path, inspect_show_labels)

    t0 = time.time()
    _ = model.fit_generator(generate_train, steps_per_epoch=training_steps, validation_steps=1,
                            validation_data=generate_valid)
    t1 = time.time()
    train_speed = (t1 - t0) / (batch_size * training_steps)
    logging.info('Spent:{} seconds training, batch_size:{} steps:{} Per example training speed:{}'.format(
        (t1 - t0), batch_size, training_steps, train_speed))

    t0 = time.time()
    _ = model.predict_generator(generate_valid, steps=training_steps, verbose=1)
    t1 = time.time()
    inference_speed = (t1 - t0) / (batch_size * training_steps)
    logging.info('Spent:{} seconds predicting. Per tensor inference speed:{}'.format((t1 - t0), inference_speed))

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

    for l in legend:
        legend_node = pydot.Node('legend'+l, label=l, shape="box", fillcolor=legend[l])
        dot.add_node(legend_node)

    logging.info('Saving architecture diagram to:{}'.format(image_path))
    dot.write_png(image_path)


def _get_tensor_maps_for_characters(tensor_maps_in: List[TensorMap], base_model: Model):
    embed_model = make_hidden_layer_model(base_model, tensor_maps_in, 'embed')
    tm_embed = TensorMap('embed', shape=(64,), group='hidden_layer', required_inputs=tensor_maps_in.copy(), model=embed_model)
    tm_char = TensorMap('ecg_rest_next_char', shape=(len(ECG_CHAR_2_IDX),), channel_map=ECG_CHAR_2_IDX, activation='softmax', loss='categorical_crossentropy', loss_weight=1.0)
    tm_burn_in = TensorMap('ecg_rest_text', shape=(100, len(ECG_CHAR_2_IDX)), group='ecg_text', channel_map={'context': 0, 'alphabet': 1}, dependent_map=tm_char)
    return [tm_embed, tm_burn_in], [tm_char]


def get_model_inputs_outputs(model_files: List[str],
                             tensor_maps_in: List[TensorMap],
                             tensor_maps_out: List[TensorMap]) -> Dict[str, Dict[str, TensorMap]]:
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
        with h5py.File(model_file, 'r') as hd5:
            model_inputs_outputs = defaultdict(list)
            for input_tensor_map in tensor_maps_in:
                if input_tensor_map.input_name() in hd5["model_weights"]:
                    model_inputs_outputs[input_prefix].append(input_tensor_map)
            for output_tensor_map in tensor_maps_out:
                if output_tensor_map.output_name() in hd5["model_weights"]:
                    model_inputs_outputs[output_prefix].append(output_tensor_map)
            if not got_tensor_maps_for_characters and 'input_ecg_rest_text_ecg_text' in hd5["model_weights"]:
                m = load_model(model_file, custom_objects=get_metric_dict(tensor_maps_out))
                char_maps_in, char_maps_out = _get_tensor_maps_for_characters(tensor_maps_in, m)
                model_inputs_outputs[input_prefix].extend(char_maps_in)
                tensor_maps_in.extend(char_maps_in)
                model_inputs_outputs[output_prefix].extend(char_maps_out)
                tensor_maps_out.extend(char_maps_out)
                got_tensor_maps_for_characters = True
                logging.info(f"Doing char model dance:{[tm.input_name() for tm in tensor_maps_in]}")
                logging.info(f"Doing char model dance out:{[tm.output_name() for tm in tensor_maps_out]}")

        models_inputs_outputs[model_file] = model_inputs_outputs

    return models_inputs_outputs
