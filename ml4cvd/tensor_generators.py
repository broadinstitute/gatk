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
import logging
import traceback
import threading
import numpy as np
from collections import Counter
from typing import List, Dict, Tuple, Generator, Optional

from ml4cvd.defines import TENSOR_EXT
from ml4cvd.TensorMap import TensorMap

np.set_printoptions(threshold=np.inf)


class TensorGenerator(object):
    """Yield minibatches of tensors given lists of I/O TensorMaps in a thread-safe way"""
    def __init__(self, batch_size, input_maps, output_maps, paths, weights=None, keep_paths=False, mixup=0.0):
        self.lock = threading.Lock()
        if weights is None:
            self.generator = multimodal_multitask_generator(batch_size, input_maps, output_maps, paths, keep_paths, mixup)
        else:
            self.generator = multimodal_multitask_weighted_generator(batch_size, input_maps, output_maps, paths, weights, keep_paths, mixup)

    def __next__(self):
        self.lock.acquire()
        try:
            return next(self.generator)
        finally:
            self.lock.release()
                    

def multimodal_multitask_generator(batch_size, input_maps, output_maps, train_paths, keep_paths, mixup_alpha):
    """Generalized data generator of input and output tensors for feed-forward networks.

    The `modes` are the different inputs, and the `tasks` are given by the outputs.
    Infinitely loops over all examples yielding batch_size input and output tensor mappings.

    Arguments:
        batch_size: number of examples in each minibatch
        input_maps: list of TensorMaps that are input names to a model
        output_maps: list of TensorMaps that are output from a model
        train_paths: list of hd5 tensors shuffled after every loop or epoch
        keep_paths: If true will also yield the paths to the tensors used in this batch
        mixup_alpha: If positive, mixup batches and use this value as shape parameter alpha

    Yields:
        A tuple of dicts for the tensor mapping tensor names to numpy arrays
        {input_1:input_tensor1, input_2:input_tensor2, ...}, {output_name1:output_tensor1, ...}
        if include_paths_in_batch is True the 3rd element of tuple will be the list of paths

    Returns:
        Never!
    """
    assert len(train_paths) > 0

    stats = Counter()
    paths_in_batch = []
    if mixup_alpha > 0:
        batch_size *= 2
    in_batch = {tm.input_name(): np.zeros((batch_size,)+tm.shape) for tm in input_maps}
    out_batch = {tm.output_name(): np.zeros((batch_size,)+tm.shape) for tm in output_maps}

    while True:
        for tp in train_paths:
            try:
                with h5py.File(tp, 'r') as hd5:
                    dependents = {}
                    for tm in input_maps:
                        in_batch[tm.input_name()][stats['batch_index']] = tm.tensor_from_file(tm, hd5, dependents)
                    
                    for tm in output_maps:
                        if tm in dependents:
                            out_batch[tm.output_name()][stats['batch_index']] = dependents[tm]
                        else:
                            out_batch[tm.output_name()][stats['batch_index']] = tm.tensor_from_file(tm, hd5)

                    paths_in_batch.append(tp)
                    stats['batch_index'] += 1
                    stats['Tensors presented'] += 1
                    if stats['batch_index'] == batch_size:
                        if mixup_alpha > 0 and keep_paths:
                            yield _mixup_batch(in_batch, out_batch, mixup_alpha) + (paths_in_batch[:batch_size//2],)
                        elif mixup_alpha > 0:
                            yield _mixup_batch(in_batch, out_batch, mixup_alpha)
                        elif keep_paths:
                            yield in_batch, out_batch, paths_in_batch
                        else:
                            yield in_batch, out_batch
                        stats['batch_index'] = 0
                        paths_in_batch = []

            except IndexError:
                stats[f"IndexError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
            except KeyError:
                stats[f"KeyError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
            except ValueError:
                stats[f"ValueError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
            except OSError:
                stats[f"OSError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
            except RuntimeError:
                stats[f"RuntimeError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
            _log_first_error(stats, tp)

        stats['epochs'] += 1
        np.random.shuffle(train_paths)
        for k in stats:
            logging.info("{}: {}".format(k, stats[k]))
        logging.info(f"Generator looped & shuffled over {len(train_paths)} tensors.")
        logging.info(f"True epoch number:{stats['epochs']} in which {int(stats['Tensors presented']/stats['epochs'])} tensors were presented.")
        if stats['Tensors presented'] == 0:
            raise ValueError(f"Completed an epoch but did not find any tensors to yield")


def multimodal_multitask_weighted_generator(batch_size, input_maps, output_maps, paths_lists, weights, keep_paths, mixup_alpha):
    """Generalized data generator of input and output tensors for feed-forward networks.

    The `modes` are the different inputs, and the `tasks` are given by the outputs.
    Infinitely loops over all examples yielding batch_size input and output tensor mappings.

    Arguments:
        batch_size: number of examples in each minibatch
        input_maps: list of TensorMaps that are input names to a model
        output_maps: list of TensorMaps that are output from a model
        paths_lists: list of lists of hd5 tensors shuffled after every loop or epoch
        weights: list of weights between (0, 1) and in total summing to 1 (i.e a distribution) for how much of each batch 
            should come from each of the lists in train_paths_lists.  Must be the same size as train_paths_lists.
        keep_paths: If true will also yield the paths to the tensors used in this batch
        mixup_alpha: If positive, mixup batches and use this value as shape parameter alpha

    Yields:
        A tuple of dicts for the tensor mapping tensor names to numpy arrays
        {input_1:input_tensor1, input_2:input_tensor2, ...}, {output_name1:output_tensor1, ...}
        if include_paths_in_batch is True the 3rd element of tuple will be the list of paths

    Returns:
        Never!
    """ 
    assert len(paths_lists) > 0 and len(paths_lists) == len(weights)

    stats = Counter()
    paths_in_batch = []
    if mixup_alpha > 0:
        batch_size *= 2
    in_batch = {tm.input_name(): np.zeros((batch_size,)+tm.shape) for tm in input_maps}
    out_batch = {tm.output_name(): np.zeros((batch_size,)+tm.shape) for tm in output_maps}
    samples = [int(w*batch_size) for w in weights]

    while True:
        for i, (tensor_list, num_samples) in enumerate(zip(paths_lists, samples)):
            for tp in np.random.choice(tensor_list, num_samples):
                try:
                    with h5py.File(tp, 'r') as hd5:
                        dependents = {}
                        for tm in input_maps:
                            in_batch[tm.input_name()][stats['batch_index']] = tm.tensor_from_file(tm, hd5, dependents)
                        
                        for tm in output_maps:
                            if tm in dependents:
                                out_batch[tm.output_name()][stats['batch_index']] = dependents[tm]
                            else:
                                out_batch[tm.output_name()][stats['batch_index']] = tm.tensor_from_file(tm, hd5)

                        paths_in_batch.append(tp)
                        stats['batch_index'] += 1
                        stats['Tensors presented from list '+str(i)] += 1
                        stats['train_paths_' + str(i)] += 1
                        if stats['batch_index'] == batch_size:
                            if mixup_alpha > 0 and keep_paths:
                                yield _mixup_batch(in_batch, out_batch, mixup_alpha, permute_first=True) + (paths_in_batch[:batch_size // 2],)
                            elif mixup_alpha > 0:
                                yield _mixup_batch(in_batch, out_batch, mixup_alpha, permute_first=True)
                            elif keep_paths:
                                yield in_batch, out_batch, paths_in_batch
                            else:
                                yield in_batch, out_batch
                            stats['batch_index'] = 0
                            paths_in_batch = []

                except IndexError as e:
                    stats['IndexError:'+str(e)] += 1
                except KeyError as e:
                    stats['Key Problems from list ' + str(i)] += 1
                    stats['KeyError:'+str(e)] += 1
                except OSError as e:
                    stats['OSError:'+str(e)] += 1
                except ValueError as e:
                    stats['ValueError:'+str(e)] += 1
                except RuntimeError:
                    stats[f"RuntimeError while attempting to generate tensor:\n{traceback.format_exc()}\n"] += 1
                _log_first_error(stats, tp)
        
        for i, tensor_list in enumerate(paths_lists):
            if len(tensor_list) <= stats['train_paths_'+str(i)]:
                stats['epochs_list_number_'+str(i)] += 1
                stats['train_paths_'+str(i)] = 0
                for k in stats:
                    logging.info(f"{k} has: {stats[k]}")
                logging.info(f"Generator looped over {len(tensor_list)} tensors from CSV group {i}.")


def big_batch_from_minibatch_generator(tensor_maps_in, tensor_maps_out, generator, minibatches, keep_paths=True):
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
    input_tensors = {tm.input_name(): [] for tm in tensor_maps_in}
    output_tensors = {tm.output_name(): [] for tm in tensor_maps_out}
    paths = []

    for _ in range(minibatches):
        next_batch = next(generator)
        for key in input_tensors:
            input_tensors[key].extend(np.copy(next_batch[0][key]))
        for key in output_tensors:
            output_tensors[key].extend(np.copy(next_batch[1][key]))
        if keep_paths:
            paths.extend(next_batch[2])
    for key in input_tensors:
        input_tensors[key] = np.array(input_tensors[key])
        logging.info("Input tensor '{}' has shape {}".format(key, input_tensors[key].shape))
    for key in output_tensors:
        output_tensors[key] = np.array(output_tensors[key])
        logging.info("Output tensor '{}' has shape {}".format(key, output_tensors[key].shape))

    if keep_paths:
        return input_tensors, output_tensors, paths
    else:
        return input_tensors, output_tensors


def get_test_train_valid_paths(tensors, valid_ratio, test_ratio, test_modulo):
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

    for root, dirs, files in os.walk(tensors):
        for name in files:
            if os.path.splitext(name)[-1].lower() != TENSOR_EXT:
                continue
            dice = np.random.rand()
            if dice < valid_ratio or (test_modulo > 1 and int(os.path.splitext(name)[0]) % test_modulo == 0):
                test_paths.append(os.path.join(root, name))
            elif dice < (valid_ratio+test_ratio):
                valid_paths.append(os.path.join(root, name))
            else:   
                train_paths.append(os.path.join(root, name))

    logging.info(f"Found {len(train_paths)} training, {len(valid_paths)} validation, and {len(test_paths)} testing tensors at: {tensors}")
    if len(train_paths) == 0 or len(valid_paths) == 0 or len(test_paths) == 0:
        raise ValueError(f"Not enough tensors at {tensors}\n")

    return train_paths, valid_paths, test_paths


def get_test_train_valid_paths_split_by_csvs(tensors, balance_csvs, valid_ratio, test_ratio, test_modulo):
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
                
            sample_id = os.path.basename(splits[0])
            group = 0
            if sample_id in sample2group:
                group = sample2group[sample_id]
            dice = np.random.rand()
            if dice < valid_ratio or (test_modulo > 1 and int(os.path.splitext(name)[0]) % test_modulo == 0):
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


def test_train_valid_tensor_generators(maps_in: List[TensorMap],
                                       maps_out: List[TensorMap],
                                       tensors: str,
                                       batch_size: int,
                                       valid_ratio: float,
                                       test_ratio: float,
                                       test_modulo: int,
                                       balance_csvs: List[str],
                                       keep_paths: bool = False,
                                       keep_paths_test: bool = True,
                                       mixup_alpha: float = -1.0) -> Tuple[
        Generator[Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Optional[List[str]]], None, None],
        Generator[Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Optional[List[str]]], None, None],
        Generator[Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], Optional[List[str]]], None, None]]:
    """ Get 3 tensor generator functions for training, validation and testing data.

    :param maps_in: list of TensorMaps that are input names to a model
    :param maps_out: list of TensorMaps that are output from a model
    :param tensors: directory containing tensors
    :param batch_size: number of examples in each mini-batch
    :param valid_ratio: rate of tensors to use for validation
    :param test_ratio: rate of tensors to use for testing
    :param test_modulo: if greater than 1, all sample ids modulo this number will be used for testing regardless of test_ratio and valid_ratio
    :param balance_csvs: if not empty, generator will provide batches balanced amongst the Sample ID in these CSVs.
    :param keep_paths: also return the list of tensor files loaded for training and validation tensors
    :param keep_paths_test:  also return the list of tensor files loaded for testing tensors
    :return: A tuple of three generators. Each yields a Tuple of dictionaries of input and output numpy arrays for training, validation and testing.
    """
    if len(balance_csvs) > 0:
        train_paths, valid_paths, test_paths = get_test_train_valid_paths_split_by_csvs(tensors, balance_csvs, valid_ratio, test_ratio, test_modulo)
        weights = [1.0/(len(balance_csvs)+1) for _ in range(len(balance_csvs)+1)]
        generate_train = TensorGenerator(batch_size, maps_in, maps_out, train_paths, weights, keep_paths, mixup_alpha)
        generate_valid = TensorGenerator(batch_size, maps_in, maps_out, valid_paths, weights, keep_paths)
        generate_test = TensorGenerator(batch_size, maps_in, maps_out, test_paths, weights, keep_paths or keep_paths_test)
    else:
        train_paths, valid_paths, test_paths = get_test_train_valid_paths(tensors, valid_ratio, test_ratio, test_modulo)
        generate_train = TensorGenerator(batch_size, maps_in, maps_out, train_paths, None, keep_paths, mixup_alpha)
        generate_valid = TensorGenerator(batch_size, maps_in, maps_out, valid_paths, None, keep_paths)
        generate_test = TensorGenerator(batch_size, maps_in, maps_out, test_paths, None, keep_paths or keep_paths_test)
    return generate_train, generate_valid, generate_test


def _log_first_error(stats: Counter, tensor_path: str):
    for k in stats:
        if 'Error' in k and stats[k] == 1:
            stats[k] += 1  # Increment so we only see these messages once
            logging.info(f"At tensor path: {tensor_path}")
            logging.info(f"Got first error: {k}")


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
