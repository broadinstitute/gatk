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
import math
import h5py
import time
import logging
import traceback
import numpy as np
import pandas as pd
from collections import Counter
from multiprocessing import Process, Queue
from itertools import chain
from typing import List, Dict, Tuple, Set, Optional, Iterator, Callable, Any, Union

from ml4h.defines import TENSOR_EXT
from ml4h.TensorMap import TensorMap

np.set_printoptions(threshold=np.inf)


DEFAULT_VALID_RATIO = 0.2
DEFAULT_TEST_RATIO = 0.1

TENSOR_GENERATOR_TIMEOUT = 64
TENSOR_GENERATOR_MAX_Q_SIZE = 32

# TensorGenerator batch indices
BATCH_INPUT_INDEX, BATCH_OUTPUT_INDEX, BATCH_SAMPLE_WEIGHTS_INDEX, BATCH_PATHS_INDEX = 0, 1, 2, 3

Path = str
PathIterator = Iterator[Path]
Batch = Dict[Path, np.ndarray]
BatchFunction = Callable[[Batch, Batch, bool, List[Path], 'kwargs'], Any]


class _ShufflePaths(Iterator):

    def __init__(self, paths: List[Path]):
        self.paths = paths
        np.random.shuffle(self.paths)
        self.idx = 0

    def __next__(self):
        if self.idx >= len(self.paths):
            self.idx = 0
            np.random.shuffle(self.paths)
        path = self.paths[self.idx]
        self.idx += 1
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
    def __init__(
        self, batch_size: int, input_maps: List[TensorMap], output_maps: List[TensorMap],
        paths: Union[List[str], List[List[str]]], num_workers: int, cache_size: float, weights: List[float] = None,
        keep_paths: bool = False, mixup: float = 0.0, name: str = 'worker', siamese: bool = False,
        augment: bool = False, sample_weight: TensorMap = None,
    ):
        """
        :param paths: If weights is provided, paths should be a list of path lists the same length as weights
        """
        self.augment = augment
        self.run_on_main_thread = num_workers == 0
        self.q = None
        self.stats_q = None
        self._started = False
        self.workers = []
        self.worker_instances = []
        if num_workers == 0:
            num_workers = 1  # The one worker is the main thread
        self.batch_size, self.input_maps, self.output_maps, self.num_workers, self.cache_size, self.weights, self.name, self.keep_paths = \
            batch_size, input_maps, output_maps, num_workers, cache_size, weights, name, keep_paths
        self.true_epochs = 0
        self.stats_string = ""
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
        elif sample_weight:
            self.input_maps = input_maps[:] + [sample_weight]
            self.batch_function = _weighted_batch
            self.batch_function_kwargs = {'sample_weight': sample_weight}
        else:
            self.batch_function = _identity_batch

    def _init_workers(self):
        self.q = Queue(min(self.batch_size, TENSOR_GENERATOR_MAX_Q_SIZE))
        self.stats_q = Queue(len(self.worker_instances))
        self._started = True
        for i, (path_iter, iter_len) in enumerate(zip(self.path_iters, self.true_epoch_lens)):
            name = f'{self.name}_{i}'
            worker_instance = _MultiModalMultiTaskWorker(
                self.q,
                self.stats_q,
                self.num_workers,
                self.input_maps, self.output_maps,
                path_iter, iter_len,
                self.batch_function, self.batch_size, self.keep_paths, self.batch_function_kwargs,
                self.cache_size,
                name,
                self.augment,
            )
            self.worker_instances.append(worker_instance)
            if not self.run_on_main_thread:
                process = Process(
                    target=worker_instance.multiprocessing_worker, name=name,
                    args=(),
                )
                process.start()
                self.workers.append(process)
        logging.info(f"Started {i + 1} {self.name.replace('_', ' ')}s with cache size {self.cache_size/1e9}GB.")

    def set_worker_paths(self, paths: List[Path]):
        """In the single worker case, set the worker's paths."""
        if not self._started:
            self._init_workers()
        if not self.run_on_main_thread:
            raise ValueError('Cannot sort paths of multiprocessing workers. num_workers must be 0.')
        self.worker_instances[0].path_iter.paths = paths

    def __next__(self) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Optional[List[str]]]:
        if not self._started:
            self._init_workers()
        if self.stats_q.qsize() == self.num_workers:
            self.aggregate_and_print_stats()
        if self.run_on_main_thread:
            return next(self.worker_instances[0])
        else:
            return self.q.get(TENSOR_GENERATOR_TIMEOUT)

    def aggregate_and_print_stats(self):
        stats = Counter()
        self.true_epochs += 1
        cur_worker = 0
        while self.stats_q.qsize() != 0:
            cur_worker += 1
            worker_stats = self.stats_q.get().copy()
            for k in worker_stats:
                if stats[k] == 0 and cur_worker == 1 and ('_max' in k or '_min' in k):
                    stats[k] = worker_stats[k]
                elif '_max' in k:
                    stats[k] = max(stats[k], worker_stats[k])
                elif '_min' in k:
                    stats[k] = min(stats[k], worker_stats[k])
                else:
                    stats[k] += worker_stats[k]

        all_errors = [
            f'[{error}] - {count:.0f}'
            for error, count in sorted(stats.items(), key=lambda x: x[1], reverse=True) if 'Error' in error
        ]
        if len(all_errors) > 0:
            error_info = f'The following errors were raised:\n\t\t' + '\n\t\t'.join(all_errors)
        else:
            error_info = 'No errors raised.'

        eps = 1e-7
        for tm in self.input_maps + self.output_maps:
            if self.true_epochs != 1:
                break
            if tm.is_categorical() and tm.axes() == 1:
                n = stats[f'{tm.name}_n'] + eps
                self.stats_string = f'{self.stats_string}\nCategorical TensorMap: {tm.name} has {n:.0f} total examples.'
                for channel, index in tm.channel_map.items():
                    examples = stats[f'{tm.name}_index_{index:.0f}']
                    self.stats_string = f'{self.stats_string}\n\tLabel {channel} {examples} examples, {100 * (examples / n):0.2f}% of total.'
            elif tm.is_continuous() and tm.axes() == 1:
                sum_squared = stats[f'{tm.name}_sum_squared']
                n = stats[f'{tm.name}_n'] + eps
                n_sum = stats[f'{tm.name}_sum']
                mean = n_sum / n
                std = np.sqrt((sum_squared/n)-(mean*mean))
                self.stats_string = f'{self.stats_string}\nContinuous TensorMap: {tm.name} has {n:.0f} total examples.\n\tMean: {mean:0.2f}, '
                self.stats_string = f"{self.stats_string}Standard Deviation: {std:0.2f}, Max: {stats[f'{tm.name}_max']:0.2f}, Min: {stats[f'{tm.name}_min']:0.2f}"
            elif tm.is_time_to_event():
                sum_squared = stats[f'{tm.name}_sum_squared']
                n = stats[f'{tm.name}_n'] + eps
                n_sum = stats[f'{tm.name}_sum']
                mean = n_sum / n
                std = np.sqrt((sum_squared/n)-(mean*mean))
                self.stats_string = f"{self.stats_string}\nTime to event TensorMap: {tm.name} Total events: {stats[f'{tm.name}_events']}, "
                self.stats_string = f"{self.stats_string}\n\tMean Follow Up: {mean:0.2f}, Standard Deviation: {std:0.2f}, "
                self.stats_string = f"{self.stats_string}\n\tMax Follow Up: {stats[f'{tm.name}_max']:0.2f}, Min Follow Up: {stats[f'{tm.name}_min']:0.2f}"

        info_string = '\n\t'.join([
            f"Generator looped & shuffled over {sum(self.true_epoch_lens)} paths. Epoch: {self.true_epochs:.0f}",
            f"{stats['Tensors presented']:0.0f} tensors were presented.",
            f"{stats['skipped_paths']} paths were skipped because they previously failed.",
            f"{error_info}",
            f"{self.stats_string}",
        ])
        logging.info(f"\n!!!!>~~~~~~~~~~~~ {self.name} completed true epoch {self.true_epochs} ~~~~~~~~~~~~<!!!!\nAggregated information string:\n\t{info_string}")

    def kill_workers(self):
        if self._started and not self.run_on_main_thread:
            for worker in self.workers:
                worker.terminate()
            logging.info(f'Stopped {len(self.workers)} workers. {self.stats_string}')
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
        self.row_size = sum(np.zeros(tm.static_shape(), dtype=np.float32).nbytes for tm in set(input_tms + output_tms))
        self.nrows = min(int(max_size / self.row_size), max_rows) if self.row_size else 0
        self.autoencode_names: Dict[str, str] = {}
        for tm in input_tms:
            self.data[tm.input_name()] = np.zeros((self.nrows,) + tm.static_shape(), dtype=np.float32)
        for tm in output_tms:
            if tm in input_tms:  # Useful for autoencoders
                self.autoencode_names[tm.output_name()] = tm.input_name()
            else:
                self.data[tm.output_name()] = np.zeros((self.nrows,) + tm.static_shape(), dtype=np.float32)
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

    def __init__(
        self,
        q: Queue,
        stats_q: Queue,
        num_workers: int,
        input_maps: List[TensorMap], output_maps: List[TensorMap],
        path_iter: PathIterator, true_epoch_len: int,
        batch_function: BatchFunction, batch_size: int, return_paths: bool, batch_func_kwargs: Dict,
        cache_size: float,
        name: str,
        augment: bool,
    ):
        self.q = q
        self.stats_q = stats_q
        self.num_workers = num_workers
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

        self.in_batch = {tm.input_name(): np.zeros((batch_size,) + tm.static_shape()) for tm in input_maps}
        self.out_batch = {tm.output_name(): np.zeros((batch_size,) + tm.static_shape()) for tm in output_maps}

        self.cache = TensorMapArrayCache(cache_size, input_maps, output_maps, true_epoch_len)
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
            self._collect_stats(tm, self.dependents[tm])
            return self.hd5
        if (path, name) in self.cache:
            batch[name][idx] = self.cache[path, name]
            return self.hd5
        if self.hd5 is None:  # Don't open hd5 if everything is in the self.cache
            self.hd5 = h5py.File(path, 'r')
        tensor = tm.postprocess_tensor(tm.tensor_from_file(tm, self.hd5, self.dependents), augment=self.augment, hd5=self.hd5)
        slices = tuple(slice(min(tm.static_shape()[i], tensor.shape[i])) for i in range(len(tensor.shape)))
        batch[name][(idx,)+slices] = tensor[slices]
        if tm.cacheable:
            self.cache[path, name] = batch[name][idx]
        self._collect_stats(tm, tensor)
        return self.hd5

    def _collect_stats(self, tm, tensor):
        if tm.is_time_to_event():
            self.epoch_stats[f'{tm.name}_events'] += tensor[0]
            self._collect_continuous_stats(tm, tensor[1])
        if tm.is_categorical() and tm.axes() == 1:
            self.epoch_stats[f'{tm.name}_index_{np.argmax(tensor):.0f}'] += 1
        if tm.is_continuous() and tm.axes() == 1:
            self._collect_continuous_stats(tm, tm.rescale(tensor)[0])
        self.epoch_stats[f'{tm.name}_n'] += 1

    def _collect_continuous_stats(self, tm, rescaled):
        if 0.0 == self.epoch_stats[f'{tm.name}_max'] == self.epoch_stats[f'{tm.name}_min']:
            self.epoch_stats[f'{tm.name}_max'] = rescaled
            self.epoch_stats[f'{tm.name}_min'] = rescaled
        self.epoch_stats[f'{tm.name}_max'] = max(rescaled, self.epoch_stats[f'{tm.name}_max'])
        self.epoch_stats[f'{tm.name}_min'] = min(rescaled, self.epoch_stats[f'{tm.name}_min'])
        self.epoch_stats[f'{tm.name}_sum'] += rescaled
        self.epoch_stats[f'{tm.name}_sum_squared'] += rescaled * rescaled

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
            self.epoch_stats['Tensors presented'] += 1
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
        self.epoch_stats['epochs'] = self.stats['epochs']
        while self.stats_q.qsize() == self.num_workers:
            continue
        self.stats_q.put(self.epoch_stats)
        if self.stats['Tensors presented'] == 0:
            logging.error(f"Completed an epoch but did not find any tensors to yield")
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
                self.in_batch = {tm.input_name(): np.zeros((self.batch_size,) + tm.static_shape()) for tm in self.input_maps}
                self.out_batch = {tm.output_name(): np.zeros((self.batch_size,) + tm.static_shape()) for tm in self.output_maps}
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
    for key, batch_array in chain(first_batch[BATCH_INPUT_INDEX].items(), first_batch[BATCH_OUTPUT_INDEX].items()):
        shape = (batch_array.shape[0] * minibatches,) + batch_array.shape[1:]
        saved_tensors[key] = np.zeros(shape)
        batch_size = batch_array.shape[0]
        saved_tensors[key][:batch_size] = batch_array

    keep_paths = generator.keep_paths
    if keep_paths:
        paths = first_batch[BATCH_PATHS_INDEX]

    input_tensors, output_tensors = list(first_batch[BATCH_INPUT_INDEX]), list(first_batch[BATCH_OUTPUT_INDEX])
    for i in range(1, minibatches):
        logging.debug(f'big_batch_from_minibatch {100 * i / minibatches:.2f}% done.')
        next_batch = next(generator)
        s, t = i * batch_size, (i + 1) * batch_size
        for key in input_tensors:
            saved_tensors[key][s:t] = next_batch[BATCH_INPUT_INDEX][key]
        for key in output_tensors:
            saved_tensors[key][s:t] = next_batch[BATCH_OUTPUT_INDEX][key]
        if keep_paths:
            paths.extend(next_batch[BATCH_PATHS_INDEX])

    for key, array in saved_tensors.items():
        logging.info(f"Made a big batch of tensors with key:{key} and shape:{array.shape}.")
    inputs = {key: saved_tensors[key] for key in input_tensors}
    outputs = {key: saved_tensors[key] for key in output_tensors}
    if keep_paths:
        return inputs, outputs, paths
    else:
        return inputs, outputs


def _get_train_valid_test_discard_ratios(
        valid_ratio: float,
        test_ratio: float,
        train_csv: str,
        valid_csv: str,
        test_csv: str,
) -> Tuple[int, int, int, int]:

    if valid_csv is not None:
        valid_ratio = 0
    if test_csv is not None:
        test_ratio = 0
    if train_csv is not None:
        train_ratio = 0
        discard_ratio = 1.0 - valid_ratio - test_ratio
    else:
        train_ratio = 1.0 - valid_ratio - test_ratio
        discard_ratio = 0

    if not math.isclose(train_ratio + valid_ratio + test_ratio + discard_ratio, 1.0):
        raise ValueError(f'ratios do not sum to 1, train/valid/test/discard = {train_ratio}/{valid_ratio}/{test_ratio}/{discard_ratio}')
    logging.debug(f'train/valid/test/discard ratios: {train_ratio}/{valid_ratio}/{test_ratio}/{discard_ratio}')

    return train_ratio, valid_ratio, test_ratio, discard_ratio


def _sample_csv_to_set(sample_csv: Optional[str] = None) -> Union[None, Set[str]]:
    if sample_csv is None:
        return None

    # Read CSV to dataframe and assume no header
    df = pd.read_csv(sample_csv, header=None)

    # If first row and column is castable to int, there is no header
    try:
        int(df.iloc[0].values[0])
    # If fails, must be header; overwrite column name with first row and remove first row
    except ValueError:
        df.columns = df.iloc[0]
        df = df[1:]

    # Declare set of possible MRN column names
    possible_mrn_col_names = {"sampleid", "medrecn", "mrn", "patient_id"}

    # Find intersection between CSV columns and possible MRN column names
    matches = set(df.columns).intersection(possible_mrn_col_names)

    # If no matches, assume the first column is MRN
    if not matches:
        mrn_col_name = df.columns[0]
    else:
         # Get first string from set of matches to use as column name
        mrn_col_name = next(iter(matches))

    if len(matches) > 1:
        logging.warning(
            f"{sample_csv} has more than one potential column for MRNs. Inferring most likely column name, but recommend explicitly setting MRN column name.",
        )

    # Isolate this column from the dataframe, and cast to strings
    sample_ids = df[mrn_col_name].apply(str)

    return set(sample_ids)



def get_train_valid_test_paths(
        tensors: str,
        sample_csv: str,
        valid_ratio: float,
        test_ratio: float,
        train_csv: str,
        valid_csv: str,
        test_csv: str,
) -> Tuple[List[str], List[str], List[str]]:
    """
    Return 3 disjoint lists of tensor paths.

    The paths are split in training, validation, and testing lists.
    If no arguments are given, paths are split into train/valid/test in the ratio 0.7/0.2/0.1.
    Otherwise, at least 2 arguments are required to specify train/valid/test sets.

    :param tensors: path to directory containing tensors
    :param sample_csv: path to csv containing sample ids, only consider sample ids for splitting
                       into train/valid/test sets if they appear in sample_csv
    :param valid_ratio: rate of tensors in validation list, mutually exclusive with valid_csv
    :param test_ratio: rate of tensors in testing list, mutually exclusive with test_csv
    :param train_csv: path to csv containing sample ids to reserve for training list
    :param valid_csv: path to csv containing sample ids to reserve for validation list, mutually exclusive with valid_ratio
    :param test_csv: path to csv containing sample ids to reserve for testing list, mutually exclusive with test_ratio

    :return: tuple of 3 lists of hd5 tensor file paths
    """
    train_paths = []
    valid_paths = []
    test_paths = []
    discard_paths = []

    train_ratio, valid_ratio, test_ratio, discard_ratio = _get_train_valid_test_discard_ratios(
        valid_ratio=valid_ratio,
        test_ratio=test_ratio,
        train_csv=train_csv,
        valid_csv=valid_csv,
        test_csv=test_csv,
    )

    choices = {
        'train': (train_paths, train_ratio),
        'valid': (valid_paths, valid_ratio),
        'test': (test_paths, test_ratio),
        'discard': (discard_paths, discard_ratio),
    }

    # parse csv's to disjoint sets, None if csv was None
    sample_set = _sample_csv_to_set(sample_csv)

    train_set = _sample_csv_to_set(train_csv)
    valid_set = _sample_csv_to_set(valid_csv)
    test_set = _sample_csv_to_set(test_csv)

    if train_set is not None and valid_set is not None and not train_set.isdisjoint(valid_set):
        raise ValueError('train and validation samples overlap')
    if train_set is not None and test_set is not None and not train_set.isdisjoint(test_set):
        raise ValueError('train and test samples overlap')
    if valid_set is not None and test_set is not None and not valid_set.isdisjoint(test_set):
        raise ValueError('validation and test samples overlap')

    # find tensors and split them among train/valid/test
    for root, dirs, files in os.walk(tensors):
        for name in files:
            path = os.path.join(root, name)
            split = os.path.splitext(name)
            sample_id = split[0]

            if split[-1].lower() != TENSOR_EXT:
                continue
            elif sample_set is not None and sample_id not in sample_set:
                continue
            elif train_set is not None and sample_id in train_set:
                train_paths.append(path)
            elif valid_set is not None and sample_id in valid_set:
                valid_paths.append(path)
            elif test_set is not None and sample_id in test_set:
                test_paths.append(path)
            else:
                choice = np.random.choice([k for k in choices], p=[choices[k][1] for k in choices])
                choices[choice][0].append(path)

    logging.info(f'Found {len(train_paths)} train, {len(valid_paths)} validation, and {len(test_paths)} testing tensors at: {tensors}')
    logging.debug(f'Discarded {len(discard_paths)} tensors due to given ratios')
    if len(train_paths) == 0 or len(valid_paths) == 0 or len(test_paths) == 0:
        raise ValueError(
            f'Not enough tensors at {tensors}\n'
            f'Found {len(train_paths)} training, {len(valid_paths)} validation, and {len(test_paths)} testing tensors\n'
            f'Discarded {len(discard_paths)} tensors',
        )

    return train_paths, valid_paths, test_paths


def get_train_valid_test_paths_split_by_csvs(
        tensors: str,
        balance_csvs: List[str],
        sample_csv: str,
        valid_ratio: float,
        test_ratio: float,
        train_csv: str,
        valid_csv: str,
        test_csv: str,
) -> Tuple[List[List[str]], List[List[str]], List[List[str]]]:
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

    train_paths = [[] for _ in range(len(balance_csvs)+1)]
    valid_paths = [[] for _ in range(len(balance_csvs)+1)]
    test_paths = [[] for _ in range(len(balance_csvs)+1)]

    _train, _valid, _test = get_train_valid_test_paths(
        tensors=tensors,
        sample_csv=sample_csv,
        valid_ratio=valid_ratio,
        test_ratio=test_ratio,
        train_csv=train_csv,
        valid_csv=valid_csv,
        test_csv=test_csv,
    )

    for paths, split_list in [(_train, train_paths), (_valid, valid_paths), (_test, test_paths)]:
        for path in paths:
            split = os.path.splitext(os.path.basename(path))
            sample_id = split[0]

            group = 0
            if sample_id in sample2group:
                group = sample2group[sample_id]
            split_list[group].append(path)

    for i in range(len(train_paths)):
        if len(train_paths[i]) == 0 or len(valid_paths[i]) == 0 or len(test_paths[i]) == 0:
            my_error = f"Not enough tensors at {tensors}\nGot {len(train_paths[i])} train {len(valid_paths[i])} valid and {len(test_paths[i])} test."
            raise ValueError(my_error)
        if i == 0:
            logging.info(f"Found {len(train_paths[i])} train {len(valid_paths[i])} valid and {len(test_paths[i])} test tensors outside the CSVs.")
        else:
            logging.info(f"CSV:{balance_csvs[i-1]}\nhas: {len(train_paths[i])} train, {len(valid_paths[i])} valid, {len(test_paths[i])} test tensors.")
    return train_paths, valid_paths, test_paths


def test_train_valid_tensor_generators(
    tensor_maps_in: List[TensorMap],
    tensor_maps_out: List[TensorMap],
    tensor_maps_protected: List[TensorMap],
    tensors: str,
    batch_size: int,
    num_workers: int,
    training_steps: int,
    validation_steps: int,
    cache_size: float,
    balance_csvs: List[str],
    keep_paths: bool = False,
    keep_paths_test: bool = True,
    mixup_alpha: float = -1.0,
    sample_csv: str = None,
    valid_ratio: float = None,
    test_ratio: float = None,
    train_csv: str = None,
    valid_csv: str = None,
    test_csv: str = None,
    siamese: bool = False,
    sample_weight: TensorMap = None,
    **kwargs
) -> Tuple[TensorGenerator, TensorGenerator, TensorGenerator]:
    """ Get 3 tensor generator functions for training, validation and testing data.
    :param tensor_maps_in: list of TensorMaps that are input names to a model
    :param tensor_maps_out: list of TensorMaps that are output from a model
    :param tensor_maps_protected: list of TensorMaps that are sensitive to bias from a model
                                    only added to the test set
    :param tensors: directory containing tensors
    :param batch_size: number of examples in each mini-batch
    :param num_workers: number of processes spun off for training and testing. Validation uses half as many workers
    :param cache_size: size in bytes of maximum cache for EACH worker
    :param balance_csvs: if not empty, generator will provide batches balanced amongst the Sample ID in these CSVs.
    :param keep_paths: also return the list of tensor files loaded for training and validation tensors
    :param keep_paths_test:  also return the list of tensor files loaded for testing tensors
    :param mixup_alpha: If positive, mixup batches and use this value as shape parameter alpha
    :param sample_csv: CSV file of sample ids, sample ids are considered for train/valid/test only if it is in sample_csv
    :param valid_ratio: rate of tensors to use for validation, mutually exclusive with valid_csv
    :param test_ratio: rate of tensors to use for testing, mutually exclusive with test_csv
    :param train_csv: CSV file of sample ids to use for training
    :param valid_csv: CSV file of sample ids to use for validation, mutually exclusive with valid_ratio
    :param test_csv: CSV file of sample ids to use for testing, mutually exclusive with test_ratio
    :param siamese: if True generate input for a siamese model i.e. a left and right input tensors for every input TensorMap
    :param sample_weight: TensorMap that outputs a sample weight for the other tensors
    :return: A tuple of three generators. Each yields a Tuple of dictionaries of input and output numpy arrays for training, validation and testing.
    """
    generate_train, generate_valid, generate_test = None, None, None
    if len(balance_csvs) > 0:
        train_paths, valid_paths, test_paths = get_train_valid_test_paths_split_by_csvs(
            tensors=tensors,
            balance_csvs=balance_csvs,
            sample_csv=sample_csv,
            valid_ratio=valid_ratio,
            test_ratio=test_ratio,
            train_csv=train_csv,
            valid_csv=valid_csv,
            test_csv=test_csv,
        )
        weights = [1.0/(len(balance_csvs)+1) for _ in range(len(balance_csvs)+1)]
    else:
        train_paths, valid_paths, test_paths = get_train_valid_test_paths(
            tensors=tensors,
            sample_csv=sample_csv,
            valid_ratio=valid_ratio,
            test_ratio=test_ratio,
            train_csv=train_csv,
            valid_csv=valid_csv,
            test_csv=test_csv,
        )
        weights = None

    num_train_workers = int(training_steps / (training_steps + validation_steps) * num_workers) or (1 if num_workers else 0)
    num_valid_workers = int(validation_steps / (training_steps + validation_steps) * num_workers) or (1 if num_workers else 0)
    generate_train = TensorGenerator(
        batch_size, tensor_maps_in, tensor_maps_out, train_paths, num_train_workers, cache_size, weights,
        keep_paths, mixup_alpha, name='train_worker', siamese=siamese, augment=True, sample_weight=sample_weight,
    )
    generate_valid = TensorGenerator(
        batch_size, tensor_maps_in, tensor_maps_out, valid_paths, num_valid_workers, cache_size, weights,
        keep_paths, name='validation_worker', siamese=siamese, augment=False,
    )
    generate_test = TensorGenerator(
        batch_size, tensor_maps_in, tensor_maps_out+tensor_maps_protected, test_paths, num_workers, 0, weights,
        keep_paths or keep_paths_test, name='test_worker', siamese=siamese, augment=False,
    )
    return generate_train, generate_valid, generate_test


def _log_first_error(stats: Counter, tensor_path: str):
    for k in stats:
        if 'Error' in k and stats[k] == 1:
            stats[k] += 1  # Increment so we only see these messages once
            logging.debug(f"At tensor path: {tensor_path}")
            logging.debug(f"Got first error: {k}")


def _identity_batch(in_batch: Batch, out_batch: Batch, return_paths: bool, paths: List[Path]):
    sample_weights = [None] * len(out_batch)
    return (in_batch, out_batch, sample_weights, paths) if return_paths else (in_batch, out_batch, sample_weights)


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


def _weighted_batch(in_batch: Batch, out_batch: Batch, return_paths: bool, paths: List[Path], sample_weight: TensorMap):
    sample_weights = [in_batch.pop(sample_weight.input_name()).flatten()] * len(out_batch)
    return (in_batch, out_batch, sample_weights, paths) if return_paths else (in_batch, out_batch, sample_weights)
