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
from random import choices, shuffle
from contextlib import contextmanager
from multiprocessing import Process, Queue
from typing import List, Dict, Tuple, Set, Optional


from ml4cvd.defines import TENSOR_EXT
from ml4cvd.TensorMap import TensorMap

np.set_printoptions(threshold=np.inf)


TENSOR_GENERATOR_TIMEOUT = 64
TENSOR_GENERATOR_MAX_Q_SIZE = 32


class TensorGenerator:
    def __init__(self, batch_size, input_maps, output_maps, paths, num_workers, cache_size, weights=None, keep_paths=False, mixup=0.0, name='worker', siamese=False):
        """
        :param paths: If weights is provided, paths should be a list of path lists the same length as weights
        """
        self.run_on_main_thread = num_workers == 0
        self.q = None
        self._started = False
        self.workers = []
        self.batch_size, self.input_maps, self.output_maps, self.num_workers, self.cache_size, self.weights, self.keep_paths, self.mixup, self.name, self.siamese, = \
            batch_size, input_maps, output_maps, num_workers, cache_size, weights, keep_paths, mixup, name, siamese
        if num_workers == 0:
            raise NotImplementedError('Running on main thread does not work yet.')
            num_workers = 1  # The one worker is the main thread
        if weights is None:
            self.worker_path_lists = np.array_split(paths, num_workers)
        else:
            # split each path list into paths for each worker.
            # E.g. for two workers: [[p1, p2], [p3, p4, p5]] -> [[[p1], [p2]], [[p3, p4], [p5]]
            split_paths = [np.array_split(a, num_workers) for a in paths]
            # Next, each list of paths gets given to each worker. E.g. [[[p1], [p3, p4]], [[p2], [p5]]]
            self.worker_path_lists = np.swapaxes(split_paths, 0, 1)

    def _init_workers(self):
        self.q = Queue(TENSOR_GENERATOR_MAX_Q_SIZE)
        self._started = True
        if self.run_on_main_thread:
            logging.info(f"Starting {self.name} on main thread.")
            self.workers.append(multimodal_multitask_worker(
                self.q, self.batch_size, self.input_maps, self.output_maps, self.worker_path_lists[0], self.keep_paths,
                self.cache_size, self.mixup, self.name, self.siamese, self.weights, main_thread=True,
            ))
            return
        for i, worker_paths in enumerate(self.worker_path_lists):
            name = f'{self.name}_{i}'
            logging.info(f"Starting {name}.")
            process = Process(target=multimodal_multitask_worker, name=name,
                              args=(
                                  self.q, self.batch_size, self.input_maps, self.output_maps, worker_paths,
                                  self.keep_paths, self.cache_size, self.mixup, name, self.siamese, self.weights, False,
                              ))
            process.start()
            self.workers.append(process)

    def __next__(self) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Optional[List[str]]]:
        if not self._started:
            self._init_workers()
        logging.debug(f'Currently there are {self.q.qsize()} queued batches.')
        if self.run_on_main_thread:
            return next(self.workers[0])
        else:
            return self.q.get(TENSOR_GENERATOR_TIMEOUT)

    def kill_workers(self):
        if self._started and not self.run_on_main_thread:
            for worker in self.workers:
                logging.info(f'Stopping {worker.name}.')
                worker.terminate()
        self.workers = []

    def __iter__(self):  # This is so python type annotations recognize TensorGenerator as an iterator
        return self

    def __del__(self):
        self.kill_workers()


class TensorMapArrayCache:
    """
    Caches numpy arrays created by tensor maps up to a maximum number of bytes
    """

    def __init__(self, max_size, input_tms: List[TensorMap], output_tms: List[TensorMap], max_rows: int):
        self.max_size = max_size
        self.data = {}
        self.row_size = sum(np.zeros(tm.shape, dtype=np.float32).nbytes for tm in input_tms + output_tms)
        self.nrows = min(int(max_size / self.row_size), max_rows)
        for tm in filter(lambda tm: tm.cacheable, input_tms):
            self.data[tm.input_name()] = np.zeros((self.nrows,) + tm.shape, dtype=np.float32)
        for tm in filter(lambda tm: tm.cacheable, output_tms):
            if tm in input_tms:  # Useful for autoencoders
                self.data[tm.output_name()] = self.data[tm.input_name()]
            self.data[tm.output_name()] = np.zeros((self.nrows,) + tm.shape, dtype=np.float32)
        self.files_seen = Counter()  # name -> max position filled in cache
        self.key_to_index = {}  # file_path, name -> position in self.data
        self.hits = 0
        self.failed_paths: Set[str] = set()

    def __setitem__(self, key: Tuple[str, str], value):
        """
        :param key: should be a tuple file_path, name
        """
        file_path, name = key
        if key in self.key_to_index:
            self.data[name][self.key_to_index[key]] = value
            return True
        if self.files_seen[name] >= self.nrows:
            return False
        self.key_to_index[key] = self.files_seen[name]
        self.data[name][self.key_to_index[key]] = value
        self.files_seen[name] += 1
        return True

    def __getitem__(self, key: Tuple[str, str]):
        """
        :param key: should be a tuple file_path, name
        """
        file_path, name = key
        val = self.data[name][self.key_to_index[file_path, name]]
        self.hits += 1
        return val

    def __contains__(self, key: Tuple[str, str]):
        return key in self.key_to_index

    def __len__(self):
        return sum(self.files_seen.values())

    def average_fill(self):
        return np.mean(list(self.files_seen.values()) or [0]) / self.nrows


def _handle_tm(tm: TensorMap, hd5: h5py.File, cache: TensorMapArrayCache, is_input: bool, tp: str, idx: int, batch: Dict, dependents: Dict) -> h5py.File:
    name = tm.input_name() if is_input else tm.output_name()
    if tm in dependents:
        batch[name][idx] = dependents[tm]
        if tm.cacheable:
            cache[tp, name] = dependents[tm]
        return hd5
    if (tp, name) in cache:
        batch[name][idx] = cache[tp, name]
        return hd5
    if hd5 is None:  # Don't open hd5 if everything is in the cache
        hd5 = h5py.File(tp, 'r')
    tensor = tm.tensor_from_file(tm, hd5, dependents)
    batch[name][idx] = tensor
    if tm.cacheable:
        cache[tp, name] = tensor
    return hd5


def multimodal_multitask_worker(q: Queue, batch_size: int, input_maps, output_maps, generator_paths, keep_paths, cache_size,
                                mixup_alpha, name, siamese, weights=None, main_thread=False):
    """Generalized data generator of input and output tensors for feed-forward networks.

    The `modes` are the different inputs, and the `tasks` are given by the outputs.
    Infinitely loops over all examples yielding batch_size input and output tensor mappings.

    :param q: Multiprocessing queue to push finished batches to
    :param batch_size: number of examples in each minibatch
    :param input_maps: list of TensorMaps that are input names to a model
    :param output_maps: list of TensorMaps that are output from a model
    :param generator_paths: list of hd5 tensors shuffled after every loop or epoch
    :param keep_paths: If true will also yield the paths to the tensors used in this batch
    :param cache_size: cache size for this worker's cache
    :param siamese: whether to make siamese ready batches
    :param mixup_alpha: If positive, mixup batches and use this value as shape parameter alpha
    :param name: Name of the worker for logs
    :param weights: Weight to sample each list of paths from generator paths. Must be same length as generator_paths.
    :param main_thread: If true, the worker yields like a generator instead of pushing to a multiprocessing queue

    What gets enqueued:
        A tuple of dicts for the tensor mapping tensor names to numpy arrays
        {input_1:input_tensor1, input_2:input_tensor2, ...}, {output_name1:output_tensor1, ...}
        if include_paths_in_batch is True the 3rd element of tuple will be the list of paths
    """
    if len(generator_paths) == 0:
        raise ValueError(f'No paths provided for worker {name}.')

    stats = Counter()
    paths_in_batch = []
    if mixup_alpha > 0:
        batch_size *= 2
    in_batch = {tm.input_name(): np.zeros((batch_size,) + tm.shape) for tm in input_maps}
    out_batch = {tm.output_name(): np.zeros((batch_size,) + tm.shape) for tm in output_maps}

    max_rows = len(generator_paths) if weights is None else sum(map(len, generator_paths))
    cache = TensorMapArrayCache(cache_size, input_maps, output_maps, max_rows)
    logging.info(f'{name} initialized cache of size {cache.row_size * cache.nrows / 1e9:.3f} GB.')

    while True:
        simple_stats = Counter()
        start = time.time()
        if weights is not None:
            if len(weights) != len(generator_paths):
                raise ValueError('weights must be the same length as train_paths.')
            epoch_len = max((len(p) for p in generator_paths))
            paths = []
            for path_list, weight in zip(generator_paths, weights):
                paths += choices(path_list, k=int(weight * epoch_len))
            shuffle(paths)
        else:
            paths = generator_paths
        all_cacheable = all((tm.cacheable for tm in input_maps + output_maps))
        for tp in paths:
            if tp in cache.failed_paths and all_cacheable:
                simple_stats['skipped_paths'] += 1
                continue
            try:
                hd5 = None
                dependents = {}
                for tm in input_maps:
                    hd5 = _handle_tm(tm, hd5, cache, True, tp, stats['batch_index'], in_batch, dependents)
                for tm in output_maps:
                    hd5 = _handle_tm(tm, hd5, cache, False, tp, stats['batch_index'], out_batch, dependents)
                paths_in_batch.append(tp)
                stats['batch_index'] += 1
                stats['Tensors presented'] += 1
                if stats['batch_index'] == batch_size:
                    if siamese and keep_paths:
                        output = _make_batch_siamese(in_batch, out_batch) + (paths_in_batch,)
                    elif siamese:
                        output = _make_batch_siamese(in_batch, out_batch)
                    elif mixup_alpha > 0 and keep_paths:
                        output = _mixup_batch(in_batch, out_batch, mixup_alpha) + (paths_in_batch[:batch_size//2],)
                    elif mixup_alpha > 0:
                        output = _mixup_batch(in_batch, out_batch, mixup_alpha)
                    elif keep_paths:
                        output = in_batch, out_batch, paths_in_batch
                    else:
                        output = in_batch, out_batch
                    if main_thread:
                        raise NotImplementedError('Running on main thread does not work yet.')
                    else:
                        q.put(output)
                    stats['batch_index'] = 0
                    paths_in_batch = []
            except IndexError as e:
                stats[f"IndexError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
                simple_stats[str(e)] += 1
                cache.failed_paths.add(tp)
            except KeyError as e:
                stats[f"KeyError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
                simple_stats[str(e)] += 1
                cache.failed_paths.add(tp)
            except ValueError as e:
                stats[f"ValueError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
                simple_stats[str(e)] += 1
                cache.failed_paths.add(tp)
            except OSError as e:
                stats[f"OSError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
                simple_stats[str(e)] += 1
                cache.failed_paths.add(tp)
            except RuntimeError as e:
                stats[f"RuntimeError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
                simple_stats[str(e)] += 1
                cache.failed_paths.add(tp)
            finally:
                if hd5 is not None:
                    hd5.close()
                _log_first_error(stats, tp)
        stats['epochs'] += 1
        np.random.shuffle(generator_paths)
        for k in stats:
            logging.debug(f"{k}: {stats[k]}")
        error_info = '\n\t\t'.join([f'[{error}] - {count}'
                                    for error, count in sorted(simple_stats.items(), key=lambda x: x[1], reverse=True)])
        info_string = '\n\t'.join([
            f"The following errors occurred:\n\t\t{error_info}",
            f"Generator looped & shuffled over {len(generator_paths)} tensors.",
            f"{int(stats['Tensors presented']/stats['epochs'])} tensors were presented.",
            f"The cache holds {len(cache)} out of a possible {len(generator_paths) * (len(input_maps) + len(output_maps))} tensors and is {100 * cache.average_fill():.0f}% full.",
            f"So far there have been {cache.hits} cache hits.",
            f"{simple_stats['skipped_paths']} paths were skipped because they previously failed.",
            f"{(time.time() - start):.2f} seconds elapsed.",
        ])
        logging.info(f"Worker {name} - In true epoch {stats['epochs']}:\n\t{info_string}")
        if stats['Tensors presented'] == 0:
            raise ValueError(f"Completed an epoch but did not find any tensors to yield")


def big_batch_from_minibatch_generator(tensor_maps_in, tensor_maps_out, generator: TensorGenerator, minibatches, keep_paths=True, siamese=False):
    """Collect minibatches into bigger batches

    Returns a dicts of numpy arrays like the same kind as generator but with more examples.

    Arguments:
        tensor_maps_in: list of TensorMaps that are input names to a model
        tensor_maps_out: list of TensorMaps that are output from a model
        generator: TensorGenerator of minibatches
        minibatches: number of times to call generator and collect a minibatch
        keep_paths: also return the list of tensor files loaded

    Returns:
        A tuple of dicts mapping tensor names to big batches of numpy arrays mapping.
    """
    # saved tensors is a dictionary with pre-allocated numpy arrays
    batch_size = generator.batch_size // 2 if siamese else generator.batch_size
    nrows = batch_size * minibatches

    if siamese:
        output_tensors = ['output_siamese']
        saved_tensors = TensorMapArrayCache(2e9, tensor_maps_in, [], nrows).data
        for name in (tmap.input_name() for tmap in tensor_maps_in):
            allocated = saved_tensors.pop(name)
            saved_tensors[name + '_left'] = allocated
            saved_tensors[name + '_right'] = allocated.copy()
        input_tensors = list(saved_tensors.keys())
        saved_tensors['output_siamese'] = np.zeros((nrows, 1), dtype=np.float32)
    else:
        input_tensors = [tm.input_name() for tm in tensor_maps_in]
        output_tensors = [tm.output_name() for tm in tensor_maps_out]
        saved_tensors = TensorMapArrayCache(4e9, tensor_maps_in, tensor_maps_out, nrows).data
    paths = []

    for i in range(minibatches):
        logging.info(f'big_batch_from_minibatch {100 * i / minibatches:.2f}% done.')
        next_batch = next(generator)
        s, t = i * batch_size, (i + 1) * batch_size
        for key in input_tensors:
            saved_tensors[key][s:t] = next_batch[0][key]
        for key in output_tensors:
            saved_tensors[key][s:t] = next_batch[1][key]
        if keep_paths:
            paths.extend(next_batch[2])

    for key in saved_tensors:
        logging.info(f"Tensor '{key}' has shape {saved_tensors[key].shape}.")
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
    if len(train_paths) == 0 or len(valid_paths) == 0 or len(test_paths) == 0:
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


@contextmanager
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
    try:
        if len(balance_csvs) > 0:
            train_paths, valid_paths, test_paths = get_test_train_valid_paths_split_by_csvs(tensors, balance_csvs, valid_ratio, test_ratio, test_modulo, test_csv)
            weights = [1.0/(len(balance_csvs)+1) for _ in range(len(balance_csvs)+1)]
        else:
            train_paths, valid_paths, test_paths = get_test_train_valid_paths(tensors, valid_ratio, test_ratio, test_modulo, test_csv)
            weights = None
        generate_train = TensorGenerator(batch_size, tensor_maps_in, tensor_maps_out, train_paths, num_workers, cache_size, weights, keep_paths, mixup_alpha, name='train_worker', siamese=siamese)
        generate_valid = TensorGenerator(batch_size, tensor_maps_in, tensor_maps_out, valid_paths, num_workers // 2, cache_size, weights, keep_paths, name='validation_worker', siamese=siamese)
        generate_test = TensorGenerator(batch_size, tensor_maps_in, tensor_maps_out, test_paths, num_workers, cache_size, weights, keep_paths or keep_paths_test, name='test_worker', siamese=siamese)
        yield generate_train, generate_valid, generate_test
    finally:
        for generator in (generate_train, generate_valid, generate_test):
            if generator is not None:
                generator.kill_workers()
                del generator  # Get rid of the generator's caches


def _log_first_error(stats: Counter, tensor_path: str):
    for k in stats:
        if 'Error' in k and stats[k] == 1:
            stats[k] += 1  # Increment so we only see these messages once
            logging.debug(f"At tensor path: {tensor_path}")
            logging.debug(f"Got first error: {k}")


def _mixup_batch(in_batch: Dict[str, np.ndarray], out_batch: Dict[str, np.ndarray], alpha: float = 1.0, permute_first: bool = False):
    for k in in_batch:
        full_batch = in_batch[k].shape[0]
        half_batch = full_batch // 2
        break

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

    return mixed_ins, mixed_outs


def _make_batch_siamese(in_batch: Dict[str, np.ndarray], out_batch: Dict[str, np.ndarray]):
    for k in in_batch:
        half_batch = in_batch[k].shape[0] // 2
        break

    siamese_in = {k+'_left': np.zeros((half_batch,) + in_batch[k].shape[1:]) for k in in_batch}
    siamese_in.update({k+'_right': np.zeros((half_batch,) + in_batch[k].shape[1:]) for k in in_batch})
    siamese_out = {'output_siamese': np.zeros((half_batch, 1))}

    for i in range(half_batch):
        for k in in_batch:
            siamese_in[k+'_left'][i] = in_batch[k][i, ...]
            siamese_in[k+'_right'][i] = in_batch[k][half_batch + i, ...]
        random_task_key = np.random.choice(list(out_batch.keys()))
        siamese_out['output_siamese'][i] = 0 if np.array_equal(out_batch[random_task_key][i], out_batch[random_task_key][i+half_batch]) else 1

    return siamese_in, siamese_out
