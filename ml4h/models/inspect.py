
# Imports
import time
import logging
import numpy as np
from typing import Iterable

import tensorflow.keras.backend as K
from tensorflow.keras.models import Model
from tensorflow.keras.utils import model_to_dot


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
    inference_steps = max(1, training_steps//8)
    _ = model.predict(generate_valid, steps=inference_steps, verbose=1)
    t1 = time.time()
    inference_speed = (t1 - t0) / (batch_size * inference_steps)

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
    grads /= 80
    iterate = K.function([input_tensor], [x, grads])
    return iterate
