# Training models built by the model factory

import os
import logging
from typing import List, Tuple, Iterable, Union

from tensorflow.keras.callbacks import History
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping, ReduceLROnPlateau, Callback

from ml4h.plots import plot_metric_history
from ml4h.defines import IMAGE_EXT, MODEL_EXT
from ml4h.models.inspect import _inspect_model
from ml4h.models.model_factory import _get_custom_objects


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
    save_last_model: bool = False,
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
    :param plot: If true, plots the metrics for train and validation set at the end of each epoch
    :param save_last_model: If true saves the model weights from last epoch otherwise saves model with best validation loss

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
        callbacks=_get_callbacks(patience, model_file, save_last_model),
    )
    generate_train.kill_workers()
    generate_valid.kill_workers()

    logging.info('Model weights saved at: %s' % model_file)
    custom_dict = _get_custom_objects(generate_train.output_maps)
    model = load_model(model_file, custom_objects=custom_dict, compile=False)
    model.compile(optimizer='adam', loss=custom_dict['loss'])
    if plot:
        plot_metric_history(history, training_steps, run_id, os.path.dirname(model_file))
    if return_history:
        return model, history
    return model


def _get_callbacks(
    patience: int, model_file: str, save_last_model: bool
) -> List[Callback]:
    callbacks = [
        ModelCheckpoint(filepath=model_file, verbose=1, save_best_only=not save_last_model),
        EarlyStopping(monitor='val_loss', patience=patience * 3, verbose=1),
        ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=patience, verbose=1),
    ]
    return callbacks

