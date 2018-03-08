import os
import json

# Keras Imports
from keras import layers
from keras import metrics
import keras.backend as K
from keras.optimizers import Adam
from keras.models import Model, load_model
from keras.layers.convolutional import Conv1D, Conv2D,  MaxPooling1D, MaxPooling2D
from keras.callbacks import ModelCheckpoint, EarlyStopping, TensorBoard, ReduceLROnPlateau
from keras.layers import Input, Dense, Dropout, BatchNormalization, SpatialDropout2D, Activation, Flatten

from . import plots
from . import defines
from . import arguments
from . import tensor_maps


def args_and_model_from_semantics(semantics_json):
    args = arguments.parse_args()
    model = set_args_and_get_model_from_semantics(args, semantics_json)
    return args, model


def set_args_and_get_model_from_semantics(args, semantics_json):
    '''Recreate a model from a json file specifying model semantics.

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
    '''
    with open(semantics_json, 'r') as infile:
        semantics = json.load(infile)

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

    weight_path_hd5 = os.path.join(os.path.dirname(semantics_json), semantics['architecture'])
    print('Loading keras weight file from:', weight_path_hd5)
    model = load_model(weight_path_hd5, custom_objects=get_metric_dict(args.labels))
    model.summary()
    return model


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Models ~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def build_reference_annotation_model(args):
    '''Build Reference 1d CNN model for classifying variants with skip connected annotations.

    Convolutions followed by dense connection, concatenated with annotations.
    Dynamically sets input channels based on args via tensor_maps.get_tensor_channel_map_from_args(args)
    Uses the functional API.
    Prints out model summary.

    Arguments
        args.tensor_name: The name of the tensor mapping which data goes to which channels
        args.annotation_set: The variant annotation set, perhaps from a HaplotypeCaller VCF.
        args.labels: The output labels (e.g. SNP, NOT_SNP, INDEL, NOT_INDEL)

    Returns
        The keras model
    '''
    if K.image_data_format() == 'channels_last':
        channel_axis = -1
    else:
        channel_axis = 1

    channel_map = tensor_maps.get_tensor_channel_map_from_args(args)
    reference = Input(shape=(args.window_size, len(channel_map)), name=args.tensor_name)
    conv_width = 12
    conv_dropout = 0.1
    fc_dropout = 0.2
    x = Conv1D(filters=256, kernel_size=conv_width, activation="relu", kernel_initializer='he_normal')(reference)
    x = Conv1D(filters=256, kernel_size=conv_width, activation="relu", kernel_initializer='he_normal')(x)
    x = Dropout(conv_dropout)(x)
    x = Conv1D(filters=128, kernel_size=conv_width, activation="relu", kernel_initializer='he_normal')(x)
    x = Dropout(conv_dropout)(x)
    x = Flatten()(x)

    annotations = Input(shape=(len(args.annotations),), name=args.annotation_set)
    annos_normed = BatchNormalization(axis=channel_axis)(annotations)
    annos_normed_x = Dense(units=40, kernel_initializer='normal', activation='relu')(annos_normed)

    x = layers.concatenate([x, annos_normed_x], axis=channel_axis)
    x = Dense(units=40, kernel_initializer='normal', activation='relu')(x)
    x = Dropout(fc_dropout)(x)
    x = layers.concatenate([x, annos_normed], axis=channel_axis)

    prob_output = Dense(units=len(args.labels), kernel_initializer='glorot_normal', activation='softmax')(x)

    model = Model(inputs=[reference, annotations], outputs=[prob_output])

    adamo = Adam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.)

    model.compile(optimizer=adamo, loss='categorical_crossentropy', metrics=get_metrics(args.labels))
    model.summary()

    if os.path.exists(args.weights_hd5):
        model.load_weights(args.weights_hd5, by_name=True)
        print('Loaded model weights from:', args.weights_hd5)

    return model



def build_read_tensor_2d_and_annotations_model(args):
    '''Build Read Tensor 2d CNN model with variant annotations mixed in for classifying variants.

    2d Convolutions followed by dense connection mixed with annotation values.
    Dynamically sets input channels based on args via defines.total_input_channels_from_args(args)
    Uses the functional API. Supports theano or tensorflow channel ordering via K.image_data_format().
    Prints out model summary.

    Arguments
        args.window_size: Length in base-pairs of sequence centered at the variant to use as input.
        args.labels: The output labels (e.g. SNP, NOT_SNP, INDEL, NOT_INDEL)

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

    read_tensor = Input(shape=in_shape, name=args.tensor_name)

    read_conv_width = 16
    conv_dropout = 0.2
    fc_dropout = 0.3
    x = Conv2D(216, (read_conv_width, 1), padding='valid', activation="relu")(read_tensor)
    x = Conv2D(160, (1, read_conv_width), padding='valid', activation="relu")(x)
    x = Conv2D(128, (read_conv_width, 1), padding='valid', activation="relu")(x)
    x = MaxPooling2D((2,1))(x)
    x = Conv2D(96, (1, read_conv_width), padding='valid', activation="relu")(x)
    x = MaxPooling2D((2,1))(x)
    x = Dropout(conv_dropout)(x)
    x = Conv2D(64, (read_conv_width, 1), padding='valid', activation="relu")(x)
    x = MaxPooling2D((2,1))(x)
    x = Dropout(conv_dropout)(x)

    x = Flatten()(x)

    # Mix the variant annotations in
    annotations = Input(shape=(len(args.annotations),), name=args.annotation_set)
    annotations_bn = BatchNormalization(axis=1)(annotations)
    alt_input_mlp = Dense(units=16, kernel_initializer='glorot_normal', activation='relu')(annotations_bn)
    x = layers.concatenate([x, alt_input_mlp], axis=concat_axis)

    x = Dense(units=32, kernel_initializer='glorot_normal', activation='relu')(x)
    x = layers.concatenate([x, annotations_bn], axis=concat_axis)
    x = Dropout(fc_dropout)(x)

    prob_output = Dense(units=len(args.labels), kernel_initializer='glorot_normal', activation='softmax')(x)

    model = Model(inputs=[read_tensor, annotations], outputs=[prob_output])

    adamo = Adam(lr=0.0001, beta_1=0.9, beta_2=0.999, epsilon=1e-08, clipnorm=1.)
    model.compile(loss='categorical_crossentropy', optimizer=adamo, metrics=get_metrics(args.labels))

    model.summary()

    if os.path.exists(args.weights_hd5):
        model.load_weights(args.weights_hd5, by_name=True)
        print('Loaded model weights from:', args.weights_hd5)

    return model


def build_tiny_2d_annotation_model(args):
    return read_tensor_2d_annotation_model_from_args(args,
                                                     conv_width = 11,
                                                     conv_height = 5,
                                                     conv_layers = [32, 32],
                                                     conv_dropout = 0.0,
                                                     spatial_dropout = False,
                                                     max_pools = [(2,1),(8,1)],
                                                     padding='valid',
                                                     annotation_units = 10,
                                                     annotation_shortcut = False,
                                                     fc_layers = [16],
                                                     fc_dropout = 0.0,
                                                     batch_normalization = False)

def read_tensor_2d_annotation_model_from_args(args,
                                              conv_width = 6,
                                              conv_height = 6,
                                              conv_layers = [128, 128, 128, 128],
                                              conv_dropout = 0.0,
                                              spatial_dropout = True,
                                              max_pools = [(3,1), (3,3)],
                                              padding='valid',
                                              annotation_units = 16,
                                              annotation_shortcut = False,
                                              fc_layers = [64],
                                              fc_dropout = 0.0,
                                              kernel_initializer='he_normal',
                                              fc_initializer='glorot_normal',
                                              batch_normalization = False):
    '''Builds Read Tensor 2d CNN model with variant annotations mixed in for classifying variants.

    Arguments specify widths and depths of each layer.
    2d Convolutions followed by dense connection mixed with annotation values.
    Dynamically sets input channels based on args via defines.total_input_channels_from_args(args)
    Uses the functional API. Supports theano or tensorflow channel ordering.
    Prints out model summary.

    Arguments:
        args.window_size: Length in base-pairs of sequence centered at the variant to use as input.
        args.labels: The output labels (e.g. SNP, NOT_SNP, INDEL, NOT_INDEL)
        args.weights_hd5: An existing model file to load weights from
        args.channels_last: Theano->False or Tensorflow->True channel ordering flag
        conv_layers: list of number of convolutional filters in each layer
        batch_normalization: Boolean whether to apply batch normalization or not

    Returns:
        The keras model
    '''
    in_channels = tensor_maps.total_input_channels_from_args(args)
    if args.channels_last:
        in_shape = (args.read_limit, args.window_size, in_channels)
        concat_axis = -1
    else:
        in_shape = (in_channels, args.read_limit, args.window_size)
        concat_axis = 1

    x = read_tensor_in = Input(shape=in_shape, name=args.tensor_name)

    # Add convolutional layers
    max_pool_diff = len(conv_layers)-len(max_pools)
    for i,f in enumerate(conv_layers):
        if i%2 == 0:
            cur_kernel = (conv_width, 1)
        else:
            cur_kernel = (1, conv_height)

        if batch_normalization:
            x = Conv2D(f, cur_kernel, activation='linear', padding=padding, kernel_initializer=kernel_initializer)(x)
            x = BatchNormalization(axis=concat_axis)(x)
            x = Activation('relu')(x)
        else:
            x = Conv2D(f, cur_kernel, activation='relu', padding=padding, kernel_initializer=kernel_initializer)(x)

        if conv_dropout > 0 and spatial_dropout:
            x = SpatialDropout2D(conv_dropout)(x)
        elif conv_dropout > 0:
            x = Dropout(conv_dropout)(x)

        if i >= max_pool_diff:
            x = MaxPooling2D(max_pools[i-max_pool_diff])(x)

    x = Flatten()(x)

    # Mix the variant annotations in
    annotations = Input(shape=(len(args.annotations),), name=args.annotation_set)
    annotations_bn = BatchNormalization(axis=-1)(annotations)
    alt_input_mlp = Dense(units=annotation_units, kernel_initializer=fc_initializer, activation='relu')(annotations_bn)
    x = layers.concatenate([x, alt_input_mlp], axis=concat_axis)

    # Fully connected layers
    for fc_units in fc_layers:
        x = Dense(units=fc_units, kernel_initializer=fc_initializer, activation='relu')(x)
        if fc_dropout > 0:
            x = Dropout(fc_dropout)(x)

    if annotation_shortcut:
        x = layers.concatenate([x, annotations_bn], axis=concat_axis)

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

    history = model.fit_generator(generate_train,
                                  steps_per_epoch=args.training_steps, epochs=args.epochs, verbose=1,
                                  validation_steps=args.validation_steps, validation_data=generate_valid,
                                  callbacks=get_callbacks(args, save_weight_hd5))
    if args.image_dir:
        plots.plot_metric_history(history, plots.weight_path_to_title(save_weight_hd5), prefix=args.image_dir)

    serialize_model_semantics(args, save_weight_hd5)
    print('Model weights saved at: %s' % save_weight_hd5)

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
        return [metrics.categorical_accuracy] + per_class_precision(classes) + per_class_recall(classes)
    elif classes and dim == 3:
        return [metrics.categorical_accuracy] + per_class_precision_3d(classes) + per_class_recall_3d(classes)
    else:
        return [metrics.categorical_accuracy, precision, recall]




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Serialization ~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def serialize_model_semantics(args, architecture_hd5):
    '''Save a json file specifying model semantics, I/O contract.

    Arguments
        args.tensor_name: String which indicates tensor map to use (from defines.py) or None
        args.window_size: sites included in the tensor map
        args.read_limit: Maximum reads included in the tensor map
        args.annotations: List of annotations or None
        args.id: the id of the run will be the name of the semantics file
        architecture_hd5: Keras model and weights hd5 file (created with save_model())
    '''
    semantics = {
        'id' : args.id,
        'output_labels' : args.labels,
        'architecture' : os.path.basename(architecture_hd5),
        'input_symbols' : args.input_symbols,
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