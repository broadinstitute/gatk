import os
import numpy as np
import pandas as pd
from typing import Optional
from tensorflow.keras import optimizers
from tensorflow.keras import backend as K
from tensorflow.keras.models import Model
from tensorflow_addons.optimizers import RectifiedAdam, TriangularCyclicalLearningRate, Triangular2CyclicalLearningRate

from ml4h.plots import plot_find_learning_rate
from ml4h.tensor_generators import TensorGenerator


def get_optimizer(name: str, learning_rate: float, steps_per_epoch: int = None, learning_rate_schedule: str = None, optimizer_kwargs=None):
    if not optimizer_kwargs:
        optimizer_kwargs = {}
    name = str.lower(name)
    rate_or_schedule = _get_learning_rate_schedule(learning_rate, learning_rate_schedule, steps_per_epoch)
    try:
        opt = optimizers.get(name)
        opt.__init__(rate_or_schedule, **optimizer_kwargs)
        return opt
    except ValueError:
        pass
    if name in NON_KERAS_OPTIMIZERS:
        return NON_KERAS_OPTIMIZERS[name](rate_or_schedule, **optimizer_kwargs)
    raise ValueError(f'Unknown optimizer {name}.')


def _get_learning_rate_schedule(learning_rate: float, learning_rate_schedule: str = None, steps_per_epoch: int = None):
    if learning_rate_schedule is None:
        return learning_rate
    if learning_rate_schedule == 'triangular':
        return TriangularCyclicalLearningRate(
            initial_learning_rate=learning_rate / 5, maximal_learning_rate=learning_rate,
            step_size=steps_per_epoch * 5,
        )
    if learning_rate_schedule == 'triangular2':
        return Triangular2CyclicalLearningRate(
            initial_learning_rate=learning_rate / 5, maximal_learning_rate=learning_rate,
            step_size=steps_per_epoch * 5,
        )
    else:
        raise ValueError(f'Learning rate schedule "{learning_rate_schedule}" unknown.')


NON_KERAS_OPTIMIZERS = {
    'radam': RectifiedAdam,
}


def find_learning_rate(model: Model, generate_train: TensorGenerator, steps: int, output_folder: str = None) -> Optional[float]:
    """
    Finds the learning rate for the model that will decrease the loss most quickly.
    Recommended to set batch size as large as possible.
    If output_folder provided, then plot and csv will be generated.
    Based on https://sgugger.github.io/how-do-you-find-a-good-learning-rate.html
    """
    beta = .98  # parameter for smoothing of loss. Larger is more smooth.
    lrs = np.geomspace(1e-7, 1e2, steps)
    avg_loss = 0
    best_loss = np.inf
    optimizer = model.optimizer
    losses, smoothed_losses = [], []
    for i, lr in enumerate(lrs):
        K.set_value(optimizer.learning_rate, lr)
        history = model.fit(generate_train, verbose=0, steps_per_epoch=1, epochs=1)
        loss = history.history['loss'][0]
        losses.append(loss)
        avg_loss = beta * avg_loss + (1 - beta) * loss
        smoothed_loss = avg_loss / (1 - beta**(i + 1))
        smoothed_losses.append(smoothed_loss)
        best_loss = min(best_loss, smoothed_loss)
        if smoothed_loss > 2 * best_loss:
            break
    best_lr = _choose_best_learning_rate(np.array(smoothed_losses), lrs)
    if output_folder:
        pd.DataFrame({'loss': losses, 'learning_rate': lrs[:len(losses)], 'smoothed_loss': smoothed_losses, 'picked_lr': best_lr}).to_csv(os.path.join(output_folder, 'find_learning_rate.csv'), index=False)
        plot_find_learning_rate(
            learning_rates=lrs[:len(losses)], losses=losses, smoothed_losses=smoothed_losses, picked_learning_rate=best_lr,
            figure_path=output_folder,
        )
    return best_lr


def _choose_best_learning_rate(smoothed_loss: np.ndarray, lr: np.ndarray) -> float:
    burn_in = len(lr) // 10
    best_delta = -np.inf
    best_lr = None
    mean_width = 4
    for i in range(burn_in, len(smoothed_loss)):
        if smoothed_loss[i - 1] > smoothed_loss[0]:
            continue  # loss delta from previously high loss should be ignored
        delta = smoothed_loss[i - mean_width - 1: i - 1].mean() - smoothed_loss[i - mean_width: i].mean()  # positive is good, means loss decreased a lot
        if delta > best_delta:
            best_delta = delta
            best_lr = lr[i]
    return float(best_lr) if best_lr else None
