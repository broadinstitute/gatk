# tensor_generators.py
#
# On-the-fly data generation of tensors for training or prediction.
#
# October 2018
# Sam Friedman 
# sam@broadinstitute.org

# Python 2/3 friendly
from __future__ import print_function

# Imports
import os
import csv
import h5py
import time
import logging
import traceback
import numpy as np
from collections import Counter
from multiprocessing import Process, Queue
from itertools import chain
from typing import List, Dict, Tuple, Set, Optional, Iterator, Callable, Any

from ml4cvd.defines import TENSOR_EXT
from ml4cvd.TensorMap import TensorMap

np.set_printoptions(threshold=np.inf)


TENSOR_GENERATOR_TIMEOUT = 64
TENSOR_GENERATOR_MAX_Q_SIZE = 32

Path = str
PathIterator = Iterator[Path]
Batch = Dict[Path, np.ndarray]
BatchFunction = Callable[[Batch, Batch, bool, List[Path], 'kwargs'], Any]


class _ShufflePaths(Iterator):

    def __init__(self, paths: List[Path]):
        self.paths = paths
        self.paths.sort()
        self.idx = 0

    def __next__(self):
        path = self.paths[self.idx]
        self.idx += 1
        if self.idx >= len(self.paths):
            self.idx = 0
            np.random.shuffle(self.paths)
        return path


class _WeightedPaths(Iterator):

    def __init__(self, paths: List[PathIterator], weights: List[float]):
        self.paths = paths
        self.weights = weights
        if len(paths) != len(weights):
            raise ValueError('Weights must be the same length as paths.')

    def __next__(self) -> str:
        return np.random.choice(np.random.choice(self.paths, self.weights))


class TensorGenerator:
    def __init__(self, batch_size, input_maps, output_maps, paths, num_workers, cache_size,
                 weights=None, keep_paths=False, mixup=0.0, name='worker', siamese=False, augment=False):
        """
        :param paths: If weights is provided, paths should be a list of path lists the same length as weights
        """
        self.augment = augment
        self.run_on_main_thread = num_workers == 0
        self.q = None
        self._started = False
        self.workers = []
        self.worker_instances = []
        self.batch_size, self.input_maps, self.output_maps, self.num_workers, self.cache_size, self.weights, self.name, self.keep_paths = \
            batch_size, input_maps, output_maps, num_workers, cache_size, weights, name, keep_paths
        if num_workers == 0:
            num_workers = 1  # The one worker is the main thread
        if weights is None:
            worker_paths = np.array_split(paths, num_workers)
            self.true_epoch_lens = list(map(len, worker_paths))
            self.path_iters = [_ShufflePaths(p) for p in worker_paths]
        else:
            # split each path list into paths for each worker.
            # E.g. for two workers: [[p1, p2], [p3, p4, p5]] -> [[[p1], [p2]], [[p3, p4], [p5]]
            split_paths = [np.array_split(a, num_workers) for a in paths]
            # Next, each list of paths gets given to each worker. E.g. [[[p1], [p3, p4]], [[p2], [p5]]]
            worker_paths = np.swapaxes(split_paths, 0, 1)
            self.true_epoch_lens = [max(map(len, p)) for p in worker_paths]
            self.path_iters = [_WeightedPaths(p, weights) for p in worker_paths]

        self.batch_function_kwargs = {}
        if mixup > 0:
            self.batch_function = _mixup_batch
            self.batch_size *= 2
            self.batch_function_kwargs = {'alpha': mixup}
        elif siamese:
            self.batch_function = _make_batch_siamese
        else:
            self.batch_function = _identity_batch

    def _init_workers(self):
        self.q = Queue(min(self.batch_size, TENSOR_GENERATOR_MAX_Q_SIZE))
        self._started = True
        for i, (path_iter, iter_len) in enumerate(zip(self.path_iters, self.true_epoch_lens)):
            name = f'{self.name}_{i}'
            worker_instance = _MultiModalMultiTaskWorker(
                self.q,
                self.input_maps, self.output_maps,
                path_iter, iter_len,
                self.batch_function, self.batch_size, self.keep_paths, self.batch_function_kwargs,
                self.cache_size,
                name,
                self.augment,
            )
            self.worker_instances.append(worker_instance)
            if not self.run_on_main_thread:
                process = Process(target=worker_instance.multiprocessing_worker, name=name,
                                  args=())
                process.start()
                self.workers.append(process)
        logging.info(f"Started {i} {self.name.replace('_', ' ')}s with cache size {self.cache_size/1e9}GB.")

    def __next__(self) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Optional[List[str]]]:
        if not self._started:
            self._init_workers()
        logging.debug(f'Currently there are {self.q.qsize()} queued batches.')
        if self.run_on_main_thread:
            return next(self.worker_instances[0])
        else:
            return self.q.get(TENSOR_GENERATOR_TIMEOUT)

    def kill_workers(self):
        if self._started and not self.run_on_main_thread:
            for worker in self.workers:
                worker.terminate()
            logging.info(f'Stopped {len(self.workers)} workers.')
        self.workers = []

    def __iter__(self):  # This is so python type annotations recognize TensorGenerator as an iterator
        return self

    def __del__(self):
        self.kill_workers()


class TensorMapArrayCache:
    """
    Caches numpy arrays created by tensor maps up to a maximum number of bytes
    """

    def __init__(self, max_size, input_tms: List[TensorMap], output_tms: List[TensorMap], max_rows: Optional[int] = np.inf):
        input_tms = [tm for tm in input_tms if tm.cacheable]
        output_tms = [tm for tm in output_tms if tm.cacheable]
        self.max_size = max_size
        self.data = {}
        self.row_size = sum(np.zeros(tm.shape, dtype=np.float32).nbytes for tm in set(input_tms + output_tms))
        self.nrows = min(int(max_size / self.row_size), max_rows) if self.row_size else 0
        self.autoencode_names: Dict[str, str] = {}
        for tm in input_tms:
            self.data[tm.input_name()] = np.zeros((self.nrows,) + tm.shape, dtype=np.float32)
        for tm in output_tms:
            if tm in input_tms:  # Useful for autoencoders
                self.autoencode_names[tm.output_name()] = tm.input_name()
            else:
                self.data[tm.output_name()] = np.zeros((self.nrows,) + tm.shape, dtype=np.float32)
        self.files_seen = Counter()  # name -> max position filled in cache
        self.key_to_index = {}  # file_path, name -> position in self.data
        self.hits = 0
        self.failed_paths: Set[str] = set()

    def _fix_key(self, key: Tuple[str, str]) -> Tuple[str, str]:
        file_path, name = key
        return file_path, self.autoencode_names.get(name, name)

    def __setitem__(self, key: Tuple[str, str], value) -> bool:
        """
        :param key: should be a tuple file_path, name
        """
        file_path, name = self._fix_key(key)
        if key in self.key_to_index:  # replace existing value
            self.data[name][self.key_to_index[key]] = value
            return True
        if self.files_seen[name] >= self.nrows:  # cache already full
            return False
        self.key_to_index[key] = self.files_seen[name]
        self.data[name][self.key_to_index[key]] = value
        self.files_seen[name] += 1
        return True

    def __getitem__(self, key: Tuple[str, str]):
        """
        :param key: should be a tuple file_path, name
        """
        file_path, name = self._fix_key(key)
        val = self.data[name][self.key_to_index[file_path, name]]
        self.hits += 1
        return val

    def __contains__(self, key: Tuple[str, str]):
        return self._fix_key(key) in self.key_to_index

    def __len__(self):
        return sum(self.files_seen.values())

    def average_fill(self):
        return np.mean(list(self.files_seen.values()) or [0]) / self.nrows if self.nrows else 0

    def __str__(self):
        hits = f"The cache has had {self.hits} hits."
        fullness = ' - '.join(f"{name} has {count} / {self.nrows} tensors" for name, count in self.files_seen.items())
        return f'{hits} {fullness}.'


class _MultiModalMultiTaskWorker:

    def __init__(self,
                 q: Queue,
                 input_maps: List[TensorMap], output_maps: List[TensorMap],
                 path_iter: PathIterator, true_epoch_len: int,
                 batch_function: BatchFunction, batch_size: int, return_paths: bool, batch_func_kwargs: Dict,
                 cache_size: float,
                 name: str,
                 augment: bool,
                 ):
        self.q = q
        self.input_maps = input_maps
        self.output_maps = output_maps
        self.path_iter = path_iter
        self.true_epoch_len = true_epoch_len
        self.batch_function = batch_function
        self.batch_size = batch_size
        self.return_paths = return_paths
        self.batch_func_kwargs = batch_func_kwargs
        self.cache_size = cache_size
        self.name = name
        self.augment = augment

        self.stats = Counter()
        self.epoch_stats = Counter()
        self.start = time.time()
        self.paths_in_batch = []
        self.in_batch = {tm.input_name(): np.zeros((batch_size,) + tm.shape) for tm in input_maps}
        self.out_batch = {tm.output_name(): np.zeros((batch_size,) + tm.shape) for tm in output_maps}

        self.cache = TensorMapArrayCache(cache_size, input_maps, output_maps, true_epoch_len)
        logging.info(f'{name} initialized cache of size {self.cache.row_size * self.cache.nrows / 1e9:.3f} GB.')

        self.dependents = {}
        self.idx = 0

    def _handle_tm(self, tm: TensorMap, is_input: bool, path: Path) -> h5py.File:
        name = tm.input_name() if is_input else tm.output_name()
        batch = self.in_batch if is_input else self.out_batch
        idx = self.stats['batch_index']
        if tm in self.dependents:
            batch[name][idx] = self.dependents[tm]
            if tm.cacheable:
                self.cache[path, name] = self.dependents[tm]
            return self.hd5
        if (path, name) in self.cache:
            batch[name][idx] = self.cache[path, name]
            return self.hd5
        if self.hd5 is None:  # Don't open hd5 if everything is in the self.cache
            self.hd5 = h5py.File(path, 'r')
        tensor = tm.postprocess_tensor(tm.tensor_from_file(tm, self.hd5, self.dependents), augment=self.augment)
        batch[name][idx] = tensor
        if tm.cacheable:
            self.cache[path, name] = tensor
        return self.hd5

    def _handle_tensor_path(self, path: Path) -> None:
        hd5 = None
        if path in self.cache.failed_paths:
            self.epoch_stats['skipped_paths'] += 1
            return
        try:
            self.dependents = {}
            self.hd5 = None
            for tm in self.input_maps:
                hd5 = self._handle_tm(tm, True, path)
            for tm in self.output_maps:
                hd5 = self._handle_tm(tm, False, path)
            self.paths_in_batch.append(path)
            self.stats['Tensors presented'] += 1
            self.stats['batch_index'] += 1
        except (IndexError, KeyError, ValueError, OSError, RuntimeError) as e:
            error_name = type(e).__name__
            self.stats[f"{error_name} while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
            self.epoch_stats[f"{error_name}: {e}"] += 1
            self.cache.failed_paths.add(path)
            _log_first_error(self.stats, path)
        finally:
            if hd5 is not None:
                hd5.close()

    def _on_epoch_end(self):
        self.stats['epochs'] += 1
        for k in self.stats:
            logging.debug(f"{k}: {self.stats[k]}")
        error_info = '\n\t\t'.join([f'[{error}] - {count}'
                                    for error, count in sorted(self.epoch_stats.items(), key=lambda x: x[1], reverse=True)])
        info_string = '\n\t'.join([
            f"The following errors occurred:\n\t\t{error_info}",
            f"Generator looped & shuffled over {self.true_epoch_len} paths.",
            f"{int(self.stats['Tensors presented']/self.stats['epochs'])} tensors were presented.",
            f"{self.epoch_stats['skipped_paths']} paths were skipped because they previously failed.",
            str(self.cache),
            f"{(time.time() - self.start):.2f} seconds elapsed.",
        ])
        logging.info(f"Worker {self.name} - In true epoch {self.stats['epochs']}:\n\t{info_string}")
        if self.stats['Tensors presented'] == 0:
            raise ValueError(f"Completed an epoch but did not find any tensors to yield")
        if 'test' in self.name:
            logging.warning(f'Test worker {self.name} completed a full epoch. Test results may be double counting samples.')
        self.start = time.time()
        self.epoch_stats = Counter()

    def multiprocessing_worker(self):
        for i, path in enumerate(self.path_iter):
            self._handle_tensor_path(path)
            if self.stats['batch_index'] == self.batch_size:

                out = self.batch_function(self.in_batch, self.out_batch, self.return_paths, self.paths_in_batch, **self.batch_func_kwargs)
                self.q.put(out)
                self.paths_in_batch = []
                self.stats['batch_index'] = 0
                self.in_batch = {tm.input_name(): np.zeros((self.batch_size,) + tm.shape) for tm in self.input_maps}
                self.out_batch = {tm.output_name(): np.zeros((self.batch_size,) + tm.shape) for tm in self.output_maps}
            if i > 0 and i % self.true_epoch_len == 0:
                self._on_epoch_end()

    def __next__(self):
        while self.stats['batch_index'] < self.batch_size:
            path = next(self.path_iter)
            self._handle_tensor_path(path)
            if self.idx > 0 and self.idx % self.true_epoch_len == 0:
                self._on_epoch_end()
            self.idx += 1
        self.stats['batch_index'] = 0
        out = self.batch_function(self.in_batch, self.out_batch, self.return_paths, self.paths_in_batch, **self.batch_func_kwargs)
        self.paths_in_batch = []
        return out


def big_batch_from_minibatch_generator(generator: TensorGenerator, minibatches: int):
    """Collect minibatches into bigger batches

    Returns a dicts of numpy arrays like the same kind as generator but with more examples.

    Arguments:
        generator: TensorGenerator of minibatches
        minibatches: number of times to call generator and collect a minibatch

    Returns:
        A tuple of dicts mapping tensor names to big batches of numpy arrays mapping.
    """
    first_batch = next(generator)
    saved_tensors = {}
    batch_size = None
    for key, batch_array in chain(first_batch[0].items(), first_batch[1].items()):
        shape = (batch_array.shape[0] * minibatches,) + batch_array.shape[1:]
        saved_tensors[key] = np.zeros(shape)
        batch_size = batch_array.shape[0]
        saved_tensors[key][:batch_size] = batch_array

    keep_paths = generator.keep_paths
    if keep_paths:
        paths = first_batch[2]

    input_tensors, output_tensors = list(first_batch[0]), list(first_batch[1])
    for i in range(1, minibatches):
        logging.debug(f'big_batch_from_minibatch {100 * i / minibatches:.2f}% done.')
        next_batch = next(generator)
        s, t = i * batch_size, (i + 1) * batch_size
        for key in input_tensors:
            saved_tensors[key][s:t] = next_batch[0][key]
        for key in output_tensors:
            saved_tensors[key][s:t] = next_batch[1][key]
        if keep_paths:
            paths.extend(next_batch[2])

    for key, array in saved_tensors.items():
        logging.info(f"Made a big batch of tensors with key:{key} and shape:{array.shape}.")
    inputs = {key: saved_tensors[key] for key in input_tensors}
    outputs = {key: saved_tensors[key] for key in output_tensors}
    if keep_paths:
        return inputs, outputs, paths
    else:
        return inputs, outputs


def get_test_train_valid_paths(tensors, valid_ratio, test_ratio, test_modulo, test_csv):
    """Return 3 disjoint lists of tensor paths.

    The paths are split in training, validation and testing lists
    apportioned according to valid_ratio and test_ratio

    Arguments:
        tensors: directory containing tensors
        valid_ratio: rate of tensors in validation list
        test_ratio: rate of tensors in testing list
        test_modulo: if greater than 1, all sample ids modulo this number will be used for testing regardless of test_ratio and valid_ratio

    Returns:
        A tuple of 3 lists of hd5 tensor file paths
    """
    test_paths = []
    train_paths = []
    valid_paths = []

    assert valid_ratio > 0 and test_ratio > 0 and valid_ratio+test_ratio < 1.0

    if test_csv is not None:
        lol = list(csv.reader(open(test_csv, 'r')))
        test_dict = {l[0]: True for l in lol}
        logging.info(f'Using external test set with {len(test_dict)} examples from file:{test_csv}')
        test_ratio = 0.0
        test_modulo = 0

    for root, dirs, files in os.walk(tensors):
        for name in files:
            if os.path.splitext(name)[-1].lower() != TENSOR_EXT:
                continue

            if test_csv is not None and os.path.splitext(name)[0] in test_dict:
                test_paths.append(os.path.join(root, name))
                continue

            dice = np.random.rand()
            if dice < test_ratio or (test_modulo > 1 and int(os.path.splitext(name)[0]) % test_modulo == 0):
                test_paths.append(os.path.join(root, name))
            elif dice < (valid_ratio+test_ratio):
                valid_paths.append(os.path.join(root, name))
            else:
                train_paths.append(os.path.join(root, name))

    logging.info(f"Found {len(train_paths)} train, {len(valid_paths)} validation, and {len(test_paths)} testing tensors at: {tensors}")
    if len(train_paths) == 0 and len(valid_paths) == 0 and len(test_paths) == 0:
        raise ValueError(f"Not enough tensors at {tensors}\n")
    return train_paths, valid_paths, test_paths


def get_test_train_valid_paths_split_by_csvs(tensors, balance_csvs, valid_ratio, test_ratio, test_modulo, test_csv):
    stats = Counter()
    sample2group = {}
    for i, b_csv in enumerate(balance_csvs):
        lol = list(csv.reader(open(b_csv, 'r'), delimiter=','))
        logging.info(f"Class Balance CSV Header: {list(enumerate(lol[0]))}")

        for row in lol[1:]:
            sample_id = row[0]
            sample2group[sample_id] = i+1  # group 0 means background class
            stats['group_'+str(i+1)] += 1
    logging.info(f"Balancing with CSVs of Sample IDs stats: {stats}")

    test_paths = [[] for _ in range(len(balance_csvs)+1)]
    train_paths = [[] for _ in range(len(balance_csvs)+1)]
    valid_paths = [[] for _ in range(len(balance_csvs)+1)]
    for root, dirs, files in os.walk(tensors):
        for name in files:
            splits = os.path.splitext(name)
            if splits[-1].lower() != TENSOR_EXT:
                continue

            group = 0
            sample_id = os.path.basename(splits[0])
            if sample_id in sample2group:
                group = sample2group[sample_id]

            dice = np.random.rand()
            if dice < test_ratio or (test_modulo > 1 and int(os.path.splitext(name)[0]) % test_modulo == 0):
                test_paths[group].append(os.path.join(root, name))
            elif dice < (valid_ratio+test_ratio):
                valid_paths[group].append(os.path.join(root, name))
            else:
                train_paths[group].append(os.path.join(root, name))

    for i in range(len(train_paths)):
        if len(train_paths[i]) == 0 or len(valid_paths[i]) == 0 or len(test_paths[i]) == 0:
            my_error = f"Not enough tensors at {tensors}\nGot {len(train_paths[i])} train {len(valid_paths[i])} valid and {len(test_paths[i])} test."
            raise ValueError(my_error)
        if i == 0:
            logging.info(f"Found {len(train_paths[i])} train {len(valid_paths[i])} valid and {len(test_paths[i])} test tensors outside the CSVs.")
        else:
            logging.info(f"CSV:{balance_csvs[i-1]}\nhas: {len(train_paths[i])} train, {len(valid_paths[i])} valid, {len(test_paths[i])} test tensors.")
    return train_paths, valid_paths, test_paths


def test_train_valid_tensor_generators(tensor_maps_in: List[TensorMap],
                                       tensor_maps_out: List[TensorMap],
                                       tensors: str,
                                       batch_size: int,
                                       valid_ratio: float,
                                       test_ratio: float,
                                       test_modulo: int,
                                       num_workers: int,
                                       cache_size: float,
                                       balance_csvs: List[str],
                                       keep_paths: bool = False,
                                       keep_paths_test: bool = True,
                                       mixup_alpha: float = -1.0,
                                       test_csv: str = None,
                                       siamese: bool = False,
                                       **kwargs) -> Tuple[TensorGenerator, TensorGenerator, TensorGenerator]:
    """ Get 3 tensor generator functions for training, validation and testing data.

    :param tensor_maps_in: list of TensorMaps that are input names to a model
    :param tensor_maps_out: list of TensorMaps that are output from a model
    :param tensors: directory containing tensors
    :param batch_size: number of examples in each mini-batch
    :param valid_ratio: rate of tensors to use for validation
    :param test_ratio: rate of tensors to use for testing
    :param test_modulo: if greater than 1, all sample ids modulo this number will be used for testing regardless of test_ratio and valid_ratio
    :param num_workers: number of processes spun off for training and testing. Validation uses half as many workers
    :param cache_size: size in bytes of maximum cache for EACH worker
    :param balance_csvs: if not empty, generator will provide batches balanced amongst the Sample ID in these CSVs.
    :param keep_paths: also return the list of tensor files loaded for training and validation tensors
    :param keep_paths_test:  also return the list of tensor files loaded for testing tensors
    :param mixup_alpha: If positive, mixup batches and use this value as shape parameter alpha
    :param test_csv: CSV file of sample ids to use for testing if set will ignore test_ration and test_modulo
    :param siamese: if True generate input for a siamese model i.e. a left and right input tensors for every input TensorMap
    :return: A tuple of three generators. Each yields a Tuple of dictionaries of input and output numpy arrays for training, validation and testing.
    """
    generate_train, generate_valid, generate_test = None, None, None
    if len(balance_csvs) > 0:
        train_paths, valid_paths, test_paths = get_test_train_valid_paths_split_by_csvs(tensors, balance_csvs, valid_ratio, test_ratio, test_modulo, test_csv)
        weights = [1.0/(len(balance_csvs)+1) for _ in range(len(balance_csvs)+1)]
    else:
        train_paths, valid_paths, test_paths = get_test_train_valid_paths(tensors, valid_ratio, test_ratio, test_modulo, test_csv)
        weights = None
    generate_train = TensorGenerator(batch_size, tensor_maps_in, tensor_maps_out, train_paths, num_workers, cache_size, weights, keep_paths, mixup_alpha, name='train_worker', siamese=siamese, augment=True)
    generate_valid = TensorGenerator(batch_size, tensor_maps_in, tensor_maps_out, valid_paths, num_workers // 2, cache_size, weights, keep_paths, name='validation_worker', siamese=siamese, augment=False)
    generate_test = TensorGenerator(batch_size, tensor_maps_in, tensor_maps_out, test_paths, num_workers, 0, weights, keep_paths or keep_paths_test, name='test_worker', siamese=siamese, augment=False)
    return generate_train, generate_valid, generate_test


def _log_first_error(stats: Counter, tensor_path: str):
    for k in stats:
        if 'Error' in k and stats[k] == 1:
            stats[k] += 1  # Increment so we only see these messages once
            logging.debug(f"At tensor path: {tensor_path}")
            logging.debug(f"Got first error: {k}")


def _identity_batch(in_batch: Batch, out_batch: Batch, return_paths: bool, paths: List[Path]):
    return (in_batch, out_batch, paths) if return_paths else (in_batch, out_batch)


def _mixup_batch(in_batch: Batch, out_batch: Batch, return_paths: bool, paths: List[Path], alpha: float = 1.0, permute_first: bool = False):
    full_batch = in_batch.values().__iter__().__next__().shape[0]
    half_batch = full_batch // 2

    if permute_first:
        permuted = np.random.permutation(full_batch)
        for k in in_batch:
            in_batch[k] = in_batch[k][permuted, ...]
        for k in out_batch:
            out_batch[k] = out_batch[k][permuted, ...]

    mixed_ins = {k: np.zeros((half_batch,) + in_batch[k].shape[1:]) for k in in_batch}
    mixed_outs = {k: np.zeros((half_batch,) + out_batch[k].shape[1:]) for k in out_batch}
    for i in range(half_batch):
        weight0 = np.random.beta(alpha, alpha)
        weight1 = 1 - weight0
        for k in in_batch:
            mixed_ins[k][i] = (in_batch[k][i, ...] * weight0) + (in_batch[k][half_batch + i, ...] * weight1)
        for k in out_batch:
            mixed_outs[k][i] = (out_batch[k][i, ...] * weight0) + (out_batch[k][half_batch + i, ...] * weight1)

    return _identity_batch(mixed_ins, mixed_outs, return_paths, paths[:half_batch])


def _make_batch_siamese(in_batch: Batch, out_batch: Batch, return_paths: bool, paths: List[Path]):
    full_batch = in_batch.values().__iter__().__next__().shape[0]
    half_batch = full_batch // 2

    siamese_in = {k+'_left': np.zeros((half_batch,) + in_batch[k].shape[1:]) for k in in_batch}
    siamese_in.update({k+'_right': np.zeros((half_batch,) + in_batch[k].shape[1:]) for k in in_batch})
    siamese_out = {'output_siamese': np.zeros((half_batch, 1))}

    for i in range(half_batch):
        for k in in_batch:
            siamese_in[k+'_left'][i] = in_batch[k][i, ...]
            siamese_in[k+'_right'][i] = in_batch[k][half_batch + i, ...]
        random_task_key = np.random.choice(list(out_batch.keys()))
        siamese_out['output_siamese'][i] = 0 if np.array_equal(out_batch[random_task_key][i], out_batch[random_task_key][i+half_batch]) else 1

    return _identity_batch(siamese_in, siamese_out, return_paths, paths)
