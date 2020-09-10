# MNIST Hand written digit tensor maps
import os
import sys
import gzip
import pickle
from typing import Dict

import h5py
import numpy as np

from ml4h.TensorMap import TensorMap, Interpretation


def mnist_image_from_hd5(tm: TensorMap, hd5: h5py.File, dependents: Dict = {}) -> np.ndarray:
    return np.array(hd5['mnist_image'])


mnist_image = TensorMap('mnist_image', shape=(28, 28, 1), tensor_from_file=mnist_image_from_hd5)


def mnist_label_from_hd5(tm: TensorMap, hd5: h5py.File, dependents: Dict = {}) -> np.ndarray:
    one_hot = np.zeros(tm.shape, dtype=np.float32)
    one_hot[int(hd5['mnist_label'][0])] = 1.0
    return one_hot


mnist_label = TensorMap(
    'mnist_label', Interpretation.CATEGORICAL, tensor_from_file=mnist_label_from_hd5,
    channel_map={f'digit_{i}': i for i in range(10)},
)


def mnist_label_as_time_to_event(tm: TensorMap, hd5: h5py.File, dependents: Dict = {}) -> np.ndarray:
    tensor = np.zeros(tm.shape, dtype=np.float32)
    label = float(hd5['mnist_label'][0])
    tensor[0] = 1.0 if np.random.rand() > (label / 10) else 0.0
    tensor[1] = np.random.randint(1, 3650)
    return tensor


mnist_time_to_event = TensorMap(
    'mnist_time_to_event', Interpretation.TIME_TO_EVENT,
    tensor_from_file=mnist_label_as_time_to_event,
)


def mnist_label_as_survival_curve(tm: TensorMap, hd5: h5py.File, dependents: Dict = {}) -> np.ndarray:
    label = float(hd5['mnist_label'][0])
    has_disease = 1.0 if np.random.rand() > (label / 10) else 0.0
    days_follow_up = np.random.randint(1, 3650)

    intervals = int(tm.shape[0] / 2)
    days_per_interval = tm.days_window / intervals
    survival_then_censor = np.zeros(tm.shape, dtype=np.float32)
    for i, day_delta in enumerate(np.arange(0, tm.days_window, days_per_interval)):
        survival_then_censor[i] = float(day_delta < days_follow_up)
        if day_delta <= days_follow_up < day_delta + days_per_interval:
            survival_then_censor[intervals + i] = has_disease
    return survival_then_censor


mnist_survival_curve = TensorMap(
    'mnist_survival_curve', Interpretation.SURVIVAL_CURVE, shape=(50,),
    tensor_from_file=mnist_label_as_survival_curve,
)


def mnist_as_hd5(hd5_folder):
    train, _, _ = load_data('mnist.pkl.gz')
    mnist_images = train[0].reshape((-1, 28, 28, 1))
    if not os.path.exists(hd5_folder):
        os.makedirs(hd5_folder)
    for i, mnist_example in enumerate(mnist_images):
        with h5py.File(os.path.join(hd5_folder, f'{i}.hd5'), 'w') as hd5:
            hd5.create_dataset('mnist_image', data=mnist_example)
            hd5.create_dataset('mnist_label', data=[train[1][i]])
        if (i+1) % 5000 == 0:
            print(f'Wrote {i+1} MNIST images and labels as HD5 files')


def load_data(dataset):
    ''' Loads the dataset
    :param dataset: the path to the dataset (here MNIST)'''
    data_dir, data_file = os.path.split(dataset)
    if data_dir == "" and not os.path.isfile(dataset):
        # Check if dataset is in the data directory.
        new_path = os.path.join("data", dataset)
        if os.path.isfile(new_path) or data_file == 'mnist.pkl.gz':
            dataset = new_path

    if (not os.path.isfile(dataset)) and data_file == 'mnist.pkl.gz':
        from urllib.request import urlretrieve
        origin = ('http://www.iro.umontreal.ca/~lisa/deep/data/mnist/mnist.pkl.gz')
        print('Downloading data from %s' % origin)
        if not os.path.exists(os.path.dirname(dataset)):
            os.makedirs(os.path.dirname(dataset))
        urlretrieve(origin, dataset)

    print('loading data...')
    f = gzip.open(dataset, 'rb')
    if sys.version_info[0] == 3:
        u = pickle._Unpickler(f)
        u.encoding = 'latin1'
        train_set, valid_set, test_set = u.load()
    else:
        train_set, valid_set, test_set = pickle.load(f)
    f.close()

    return train_set, valid_set, test_set
