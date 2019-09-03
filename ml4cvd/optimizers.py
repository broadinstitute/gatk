from keras import optimizers
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


NON_KERAS_OPTIMIZERS = {
    'radam': RAdam,
}
