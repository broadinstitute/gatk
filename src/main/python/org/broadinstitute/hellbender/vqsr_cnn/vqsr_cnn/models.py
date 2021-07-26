import os
import json

# Keras Imports
import tensorflow as tf
import tensorflow.keras.backend as K
from tensorflow.keras import layers
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
#
from . import plots
from . import defines
from . import arguments
from . import tensor_maps


def start_session_get_args_and_model(intra_ops, inter_ops, semantics_json, weights_hd5=None, tensor_type=None):
    #K.clear_session()
    #K.get_session().close()
    #cfg = K.tf.ConfigProto(intra_op_parallelism_threads=intra_ops, inter_op_parallelism_threads=inter_ops)
    #cfg.gpu_options.allow_growth = True
    #K.set_session(K.tf.Session(config=cfg))
    return args_and_model_from_semantics(semantics_json, weights_hd5, tensor_type)


def args_and_model_from_semantics(semantics_json, weights_hd5=None, tensor_type=None):
    args = arguments.parse_args()

    if semantics_json is not None and os.path.exists(semantics_json):
        model = set_args_and_get_model_from_semantics(args, semantics_json, weights_hd5)
    else:
        model = load_model(weights_hd5, custom_objects=get_metric_dict(args.labels))
        args.tensor_name = tensor_type

    return args, model


def set_args_and_get_model_from_semantics(args, semantics_json, weights_hd5=None):
    """Recreate a model from a json file specifying model semantics.

    Update the args namespace from the semantics file values.
    Assert that the serialized tensor map and the recreated one are the same.

    Arguments:
        args.tensor_name: String which indicates tensor map to use or None
        args.window_size: sites included in the tensor map
        args.read_limit: Maximum reads included in the tensor map
        args.annotations: List of annotations or None
        semantics_json: Semantics json file (created with serialize_model_semantics())

    Returns:
        The Keras model
    """
    with open(semantics_json, 'r') as infile:
        semantics = json.load(infile)

    if 'model_version' in semantics:
        assert(args.model_version == semantics['model_version'])

    if 'input_tensor_map' in semantics:
        args.tensor_name = semantics['input_tensor_map_name']
        args.window_size = semantics['window_size']
        args.read_limit = semantics['read_limit']
        tm = tensor_maps.get_tensor_channel_map_from_args(args)
        assert(len(tm) == len(semantics['input_tensor_map']))
        for key in tm:
            assert(tm[key] == semantics['input_tensor_map'][key])

    if 'input_annotations' in semantics:
        args.annotations = semantics['input_annotations']
        args.annotation_set = semantics['input_annotation_set']

    args.input_symbols = semantics['input_symbols']
    args.labels = semantics['output_labels']

    if 'channels_last' in semantics:
        args.channels_last = semantics['channels_last']
        if args.channels_last:
            K.set_image_data_format('channels_last')
        else:
            K.set_image_data_format('channels_first')

    if weights_hd5 is None:
        weights_hd5 = os.path.join(os.path.dirname(semantics_json), semantics['architecture'])

    print('Updated arguments:', args, '\nWeight file from:', weights_hd5)
    model = load_model(weights_hd5, custom_objects=get_metric_dict(args.labels))
    model.summary()
    return model


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Models ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def build_default_1d_annotation_model(args):
    return build_reference_annotation_1d_model_from_args(args,
                                                         conv_width=7,
                                                         conv_layers=[256, 216, 128, 64, 32],
                                                         conv_dropout=0.1,
                                                         conv_batch_normalize=True,
                                                         spatial_dropout=True,
                                                         max_pools=[],
                                                         padding='same',
                                                         annotation_units=64,
                                                         annotation_shortcut=True,
                                                         fc_layers=[64, 64],
                                                         fc_dropout=0.2,
                                                         annotation_batch_normalize=True,
                                                         fc_batch_normalize=False)


def build_1d_annotation_model_from_args(args):
    return build_reference_annotation_1d_model_from_args(args,
                                                         conv_width=args.conv_width,
                                                         conv_layers=args.conv_layers,
                                                         conv_dropout=args.conv_dropout,
                                                         conv_batch_normalize=args.conv_batch_normalize,
                                                         spatial_dropout=args.spatial_dropout,
                                                         max_pools=args.max_pools,
                                                         padding=args.padding,
                                                         annotation_units=args.annotation_units,
                                                         annotation_shortcut=args.annotation_shortcut,
                                                         fc_layers=args.fc_layers,
                                                         fc_dropout=args.fc_dropout,
                                                         fc_batch_normalize=args.fc_batch_normalize)


def build_2d_annotation_model_from_args(args):
    return read_tensor_2d_annotation_model_from_args(args,
                                                     conv_width = args.conv_width,
                                                     conv_height = args.conv_height,
                                                     conv_layers = args.conv_layers,
                                                     conv_dropout = args.conv_dropout,
                                                     conv_batch_normalize = args.conv_batch_normalize,
                                                     spatial_dropout = args.spatial_dropout,
                                                     max_pools = args.max_pools,
                                                     padding = args.padding,
                                                     annotation_units = args.annotation_units,
                                                     annotation_shortcut = args.annotation_shortcut,
                                                     fc_layers = args.fc_layers,
                                                     fc_dropout = args.fc_dropout,
                                                     fc_batch_normalize = args.fc_batch_normalize)


def build_default_2d_annotation_model(args):
    return read_tensor_2d_annotation_model_from_args(args,
                                                     conv_width = 25,
                                                     conv_height = 25,
                                                     conv_layers = [64, 48, 32, 24],
                                                     conv_dropout = 0.1,
                                                     conv_batch_normalize = False,
                                                     spatial_dropout = True,
                                                     max_pools = [(3,1),(3,1)],
                                                     padding='valid',
                                                     annotation_units = 64,
                                                     annotation_shortcut = False,
                                                     fc_layers = [24],
                                                     fc_dropout = 0.3,
                                                     fc_batch_normalize = False)


def read_tensor_2d_annotation_model_from_args(args,
                                              conv_width = 6,
                                              conv_height = 6,
                                              conv_layers = [128, 128, 128, 128],
                                              conv_dropout = 0.0,
                                              conv_batch_normalize = False,
                                              spatial_dropout = True,
                                              residual_layers = [],
                                              max_pools = [(3,1), (3,3)],
                                              padding='valid',
                                              annotation_units = 16,
                                              annotation_shortcut = False,
                                              annotation_batch_normalize = True,
                                              fc_layers = [64],
                                              fc_dropout = 0.0,
                                              fc_batch_normalize = False,
                                              kernel_initializer='glorot_normal',
                                              kernel_single_channel=True,
                                              fc_initializer='glorot_normal'):
    '''Builds Read Tensor 2d CNN model with variant annotations mixed in for classifying variants.

    Arguments specify widths and depths of each layer.
    2d Convolutions followed by dense connection mixed with annotation values.
    Dynamically sets input channels based on args via defines.total_input_channels_from_args(args)
    Uses the functional API. Supports theano or tensorflow channel ordering.
    Prints out model summary.

    Arguments
        args.window_size: Length in base-pairs of sequence centered at the variant to use as input.
        args.labels: The output labels (e.g. SNP, NOT_SNP, INDEL, NOT_INDEL)
        args.weights_hd5: An existing model file to load weights from
        args.channels_last: Theano->False or Tensorflow->True channel ordering flag
        conv_layers: list of number of convolutional filters in each layer
        batch_normalization: Boolean whether to apply batch normalization or not
    Returns
        The keras model
    '''
    in_channels = tensor_maps.total_input_channels_from_args(args)

    if K.image_data_format() == 'channels_last':
        in_shape = (args.read_limit, args.window_size, in_channels)
        concat_axis = -1
    else:
        in_shape = (in_channels, args.read_limit, args.window_size)
        concat_axis = 1

    x = read_tensor_in = Input(shape=in_shape, name=args.tensor_name)

    max_pool_diff = max(0, len(conv_layers)-len(max_pools))

    # Add convolutional layers
    for i,f in enumerate(conv_layers):
        if kernel_single_channel and i%2 == 0:
            cur_kernel = (conv_width, 1)
        elif kernel_single_channel:
            cur_kernel = (1, conv_height)
        else:
            cur_kernel = (conv_width, conv_height)

        if conv_batch_normalize:
            x = Conv2D(int(f), cur_kernel, activation='linear', padding=padding, kernel_initializer=kernel_initializer)(x)
            x = BatchNormalization(axis=concat_axis)(x)
            x = Activation('relu')(x)
        else:
            x = Conv2D(int(f), cur_kernel, activation='relu', padding=padding, kernel_initializer=kernel_initializer)(x)

        if conv_dropout > 0 and spatial_dropout:
            x = SpatialDropout2D(conv_dropout)(x)
        elif conv_dropout > 0:
            x = Dropout(conv_dropout)(x)

        if i >= max_pool_diff:
            x = MaxPooling2D(max_pools[i-max_pool_diff])(x)

    for i,r in enumerate(residual_layers):
        if kernel_single_channel and i%2 == 0:
            cur_kernel = (conv_width, 1)
        elif kernel_single_channel:
            cur_kernel = (1, conv_height)
        else:
            cur_kernel = (conv_width, conv_height)

        y = Conv2D(r.filters[0], (1, 1), strides=r.strides)(x)
        y = BatchNormalization(axis=concat_axis)(y)
        y = Activation('relu')(y)

        y = Conv2D(r.filters[1], cur_kernel, padding='same')(y)
        y = BatchNormalization(axis=concat_axis)(y)
        y = Activation('relu')(y)

        y = Conv2D(r.filters[2], (1, 1))(y)
        y = BatchNormalization(axis=concat_axis)(y)

        if r.identity:
            x = layers.add([y, x])
        else:
            shortcut = Conv2D(r.filters[2], (1, 1), strides=r.strides)(x)
            shortcut = BatchNormalization(axis=concat_axis)(shortcut)
            x = layers.add([y, shortcut])

        x = Activation('relu')(x)

    x = Flatten()(x)

    # Mix the variant annotations in
    annotations = annotations_in = Input(shape=(len(args.annotations),), name=args.annotation_set)
    if annotation_batch_normalize:
        annotations_in = BatchNormalization(axis=-1)(annotations)

    annotations_mlp = Dense(units=annotation_units, kernel_initializer=fc_initializer, activation='relu')(annotations_in)
    x = layers.concatenate([x, annotations_mlp], axis=concat_axis)

    # Fully connected layers
    for fc_units in fc_layers:

        if fc_batch_normalize:
            x = Dense(units=fc_units, kernel_initializer=fc_initializer, activation='linear')(x)
            x = BatchNormalization(axis=1)(x)
            x = Activation('relu')(x)
        else:
            x = Dense(units=fc_units, kernel_initializer=fc_initializer, activation='relu')(x)

        if fc_dropout > 0:
            x = Dropout(fc_dropout)(x)

    if annotation_shortcut:
        x = layers.concatenate([x, annotations_in], axis=concat_axis)

    # Softmax output
    prob_output = Dense(units=len(args.labels), kernel_initializer=fc_initializer, activation='softmax')(x)

    # Map inputs to outputs
    model = Model(inputs=[read_tensor_in, annotations], outputs=[prob_output])

    adamo = Adam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.)
    model.compile(loss='categorical_crossentropy', optimizer=adamo, metrics=get_metrics(args.labels))
    model.summary()

    if os.path.exists(args.weights_hd5):
        model.load_weights(args.weights_hd5, by_name=True)
        print('Loaded model weights from:', args.weights_hd5)

    return model


def build_reference_annotation_1d_model_from_args(args,
                                                  conv_width = 6,
                                                  conv_layers = [128, 128, 128, 128],
                                                  conv_dropout = 0.0,
                                                  conv_batch_normalize = False,
                                                  spatial_dropout = True,
                                                  max_pools = [],
                                                  padding='valid',
                                                  activation = 'relu',
                                                  annotation_units = 16,
                                                  annotation_shortcut = False,
                                                  annotation_batch_normalize = True,
                                                  fc_layers = [64],
                                                  fc_dropout = 0.0,
                                                  fc_batch_normalize = False,
                                                  fc_initializer = 'glorot_normal',
                                                  kernel_initializer = 'glorot_normal',
                                                  alpha_dropout = False
                                                  ):
    '''Build Reference 1d CNN model for classifying variants.

    Architecture specified by parameters.
    Dynamically sets input channels based on args via defines.total_input_channels_from_args(args)
    Uses the functional API.
    Prints out model summary.

    Arguments
        args.annotations: The variant annotations, perhaps from a HaplotypeCaller VCF.
        args.labels: The output labels (e.g. SNP, NOT_SNP, INDEL, NOT_INDEL)

    Returns
        The keras model
    '''
    in_channels = tensor_maps.total_input_channels_from_args(args)
    concat_axis = -1
    x = reference = Input(shape=(args.window_size, in_channels), name=args.tensor_name)

    max_pool_diff = len(conv_layers)-len(max_pools)
    for i,c in enumerate(conv_layers):

        if conv_batch_normalize:
            x = Conv1D(filters=c, kernel_size=conv_width, activation='linear', padding=padding, kernel_initializer=kernel_initializer)(x)
            x = BatchNormalization(axis=concat_axis)(x)
            x = Activation(activation)(x)
        else:
            x = Conv1D(filters=c, kernel_size=conv_width, activation=activation, padding=padding, kernel_initializer=kernel_initializer)(x)

        if conv_dropout > 0 and alpha_dropout:
            x = AlphaDropout(conv_dropout)(x)
        elif conv_dropout > 0 and spatial_dropout:
            x = SpatialDropout1D(conv_dropout)(x)
        elif conv_dropout > 0:
            x = Dropout(conv_dropout)(x)

        if i >= max_pool_diff:
            x = MaxPooling1D(max_pools[i-max_pool_diff])(x)

    f = Flatten()(x)

    annotations = annotations_in = Input(shape=(len(args.annotations),), name=args.annotation_set)
    if annotation_batch_normalize:
        annotations_in = BatchNormalization(axis=concat_axis)(annotations_in)
    annotation_mlp = Dense(units=annotation_units, kernel_initializer=fc_initializer, activation=activation)(annotations_in)

    x = layers.concatenate([f, annotation_mlp], axis=1)
    for fc in fc_layers:
        if fc_batch_normalize:
            x = Dense(units=fc, activation='linear', kernel_initializer=fc_initializer)(x)
            x = BatchNormalization(axis=1)(x)
            x = Activation(activation)(x)
        else:
            x = Dense(units=fc, activation=activation, kernel_initializer=fc_initializer)(x)

        if fc_dropout > 0 and alpha_dropout:
            x = AlphaDropout(fc_dropout)(x)
        elif fc_dropout > 0:
            x = Dropout(fc_dropout)(x)

    if annotation_shortcut:
        x = layers.concatenate([x, annotations_in], axis=1)

    prob_output = Dense(units=len(args.labels), activation='softmax', name='softmax_predictions')(x)

    model = Model(inputs=[reference, annotations], outputs=[prob_output])

    adam = Adam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.)
    model.compile(optimizer=adam, loss='categorical_crossentropy', metrics=get_metrics(args.labels))
    model.summary()

    if os.path.exists(args.weights_hd5):
        model.load_weights(args.weights_hd5, by_name=True)
        print('Loaded model weights from:', args.weights_hd5)

    return model


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Optimizing ~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def train_model_from_generators(args, model, generate_train, generate_valid, save_weight_hd5):
    '''Train an image model for classifying variants.

    Training data lives on disk, it will be loaded by generator functions.
    Plots the metric history after training. Creates a directory to save weights at if necessary.

    Arguments
        args.batch_size: size of the mini-batches
        args.patience: Maximum number of epochs to run without validation loss improvement
        args.epochs: Maximum number of epochs to run regardless of Early Stopping
        args.training_steps: Number of mini-batches in each so-called epoch
        args.validation_steps: Number of validation mini-batches to examine after each epoch.
        model: the model to optimize
        generate_train: training data generator function
        valid_tuple: Validation data data generator function
        save_weight_hd5: path to save the model weights at

    Returns
        The now optimized keras model
    '''
    if not os.path.exists(os.path.dirname(save_weight_hd5)):
        os.makedirs(os.path.dirname(save_weight_hd5))
    serialize_model_semantics(args, save_weight_hd5)

    history = model.fit_generator(generate_train,
                                  steps_per_epoch=args.training_steps, epochs=args.epochs, verbose=1,
                                  validation_steps=args.validation_steps, validation_data=generate_valid,
                                  callbacks=get_callbacks(args, save_weight_hd5))
    print('Training complete, model weights saved at: %s' % save_weight_hd5)
    if args.image_dir:
        plots.plot_metric_history(history, plots.weight_path_to_title(save_weight_hd5), prefix=args.image_dir)




    return model


def get_callbacks(args, save_weight_hd5):
    callbacks = []

    callbacks.append(ModelCheckpoint(filepath=save_weight_hd5, verbose=1, save_best_only=True))
    callbacks.append(EarlyStopping(monitor='val_loss', patience=args.patience*4, verbose=1))
    callbacks.append(ReduceLROnPlateau(monitor='val_loss', patience=args.patience, verbose=1))

    if args.tensor_board:
        callbacks.append(TensorBoard())

    return callbacks


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Metrics ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def precision(y_true, y_pred):
    '''Calculates the precision, a metric for multi-label classification of
    how many selected items are relevant.
    '''
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def recall(y_true, y_pred):
    '''Calculates the recall, a metric for multi-label classification of
    how many relevant items are selected.
    '''
    true_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall


def per_class_recall(labels):
    recall_fxns = []

    for label_key in labels:
        label_idx = labels[label_key]
        fxn = 'def '+ label_key + '_recall(y_true, y_pred):\n'
        fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
        fxn += '\tpossible_positives = K.sum(K.round(K.clip(y_true, 0, 1)), axis=0)\n'
        fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'

        exec(fxn)
        recall_fxn = eval(label_key + '_recall')
        recall_fxns.append(recall_fxn)

    return recall_fxns


def per_class_precision(labels):
    precision_fxns = []

    for label_key in labels:
        label_idx = labels[label_key]
        fxn = 'def '+ label_key + '_precision(y_true, y_pred):\n'
        fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
        fxn += '\tpredicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0)\n'
        fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'

        exec(fxn)
        precision_fxn = eval(label_key + '_precision')
        precision_fxns.append(precision_fxn)
    return precision_fxns


def get_metric_dict(labels=defines.SNP_INDEL_LABELS):
    metrics = {'precision':precision, 'recall':recall}
    precision_fxns = per_class_precision(labels)
    recall_fxns = per_class_recall(labels)
    for i,label_key in enumerate(labels.keys()):
        metrics[label_key+'_precision'] = precision_fxns[i]
        metrics[label_key+'_recall'] = recall_fxns[i]
    return metrics


def per_class_recall_3d(labels):
    recall_fxns = []

    for label_key in labels:
        label_idx = labels[label_key]
        fxn = 'def '+ label_key + '_recall(y_true, y_pred):\n'
        fxn += '\ttrue_positives = K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0)\n'
        fxn += '\tpossible_positives = K.sum(K.sum(K.round(K.clip(y_true, 0, 1)), axis=0), axis=0)\n'
        fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'

        exec(fxn)
        recall_fxn = eval(label_key + '_recall')
        recall_fxns.append(recall_fxn)

    return recall_fxns


def per_class_precision_3d(labels):
    precision_fxns = []

    for label_key in labels:
        label_idx = labels[label_key]
        fxn = 'def '+ label_key + '_precision(y_true, y_pred):\n'
        fxn += '\ttrue_positives = K.sum(K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0), axis=0)\n'
        fxn += '\tpredicted_positives = K.sum(K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0), axis=0)\n'
        fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'

        exec(fxn)
        precision_fxn = eval(label_key + '_precision')
        precision_fxns.append(precision_fxn)

    return precision_fxns


def get_metrics(classes=None, dim=2):
    if classes and dim == 2:
        return per_class_precision(classes) + per_class_recall(classes)
    elif classes and dim == 3:
        return per_class_precision_3d(classes) + per_class_recall_3d(classes)
    else:
        return [precision, recall]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Serialization ~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def serialize_model_semantics(args, architecture_hd5):
    """Save a json file specifying model semantics, I/O contract.

    Arguments
        args.tensor_name: String which indicates tensor map to use (from defines.py) or None
        args.window_size: sites included in the tensor map
        args.read_limit: Maximum reads included in the tensor map
        args.annotations: List of annotations or None
        args.id: the id of the run will be the name of the semantics file
        architecture_hd5: Keras model and weights hd5 file (created with save_model())
    """
    semantics = {
        'id': args.id,
        'output_labels': args.labels,
        'architecture': os.path.basename(architecture_hd5),
        'input_symbols': args.input_symbols,
        'model_version': args.model_version,
        'gatk_version': args.gatk_version,
    }

    if args.tensor_name:
        semantics['input_tensor_map_name'] = args.tensor_name
        semantics['input_tensor_map'] = tensor_maps.get_tensor_channel_map_from_args(args)
        semantics['window_size'] = args.window_size
        semantics['read_limit'] = args.read_limit

    if args.annotation_set and args.annotation_set != '_':
        semantics['input_annotations'] = args.annotations
        semantics['input_annotation_set'] = args.annotation_set

    if args.data_dir:
        semantics['data_dir'] = args.data_dir

    semantics['channels_last'] = args.channels_last

    json_file_name = args.output_dir + args.id + '.json'
    with open(json_file_name, 'w') as outfile:
        json.dump(semantics, outfile)

    print('Saved model semantics at:', json_file_name)