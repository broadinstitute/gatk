import tensorflow as tf
from tensorflow.keras import optimizers
import os
os.environ['TF_KERAS'] = '1'  # TODO: gross but somehow necessary for RAdam
from keras_radam import RAdam


def get_optimizer(name: str, lr: float, optimizer_kwargs=None):
    if not optimizer_kwargs:
        optimizer_kwargs = {}
    name = str.lower(name)
    try:
        opt = optimizers.get(name)
        opt.__init__(lr=lr, **optimizer_kwargs)
        return opt
    except ValueError:
        pass
    if name in NON_KERAS_OPTIMIZERS:
        return NON_KERAS_OPTIMIZERS[name](lr, **optimizer_kwargs)
    raise ValueError(f'Unknown optimizer {name}')


class RAdam(RAdam):
    """
    This is for backwards compatability with old models when learning_rate was lr in keras optimizers
    TODO: replace models in tests.py so this is unnecessary
    """

    def __init__(self, lr=None, learning_rate=None, **kwargs):
        super().__init__(learning_rate=lr or learning_rate)


NON_KERAS_OPTIMIZERS = {
    'radam': RAdam,
}
