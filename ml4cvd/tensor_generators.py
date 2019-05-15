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
from typing import List
from collections import Counter

from ml4cvd.defines import TENSOR_EXT
from ml4cvd.TensorMap import TensorMap

np.set_printoptions(threshold=np.inf)


class TensorGenerator(object):
    """Yield minibatches of tensors given lists of I/O TensorMaps in a thread-safe way"""
    def __init__(self, batch_size, input_maps, output_maps, paths, weights=None, keep_paths=False):
        self.lock = threading.Lock()
        if weights is None:
            self.generator = multimodal_multitask_generator(batch_size, input_maps, output_maps, paths, keep_paths)
        else:
            self.generator = multimodal_multitask_weighted_generator(batch_size, input_maps, output_maps, paths, weights, keep_paths)

    def __next__(self):
        self.lock.acquire()
        try:
            return next(self.generator)
        finally:
            self.lock.release()
                    

def multimodal_multitask_generator(batch_size, input_maps, output_maps, train_paths, keep_paths):
    """Generalizaed data generator of input and output tensors for feed-forward networks.

    The `modes` are the different inputs, and the `tasks` are given by the outputs.
    Inifinitely loops over all examples yielding batch_size input and output tensor mappings.

    Arguments:
        batch_size: number of examples in each minibatch
        input_maps: list of TensorMaps that are input names to a model
        output_maps: list of TensorMaps that are output from a model
        train_paths: list of hd5 tensors shuffled after every loop or epoch
        keep_paths: If true will also yield the paths to the tensors used in this batch

    Yields:
        A tuple of dicts for the tensor mapping tensor names to numpy arrays
        {input_1:input_tensor1, input_2:input_tensor2, ...}, {output_name1:output_tensor1, ...}
        if include_paths_in_batch is True the 3rd element of tuple will be the list of paths

    Returns:
        Never!
    """ 
    assert(len(train_paths) > 0)

    stats = Counter()
    paths_in_batch = []
    in_batch = {tm.input_name(): np.zeros((batch_size,)+tm.shape) for tm in input_maps}
    out_batch = {tm.output_name(): np.zeros((batch_size,)+tm.shape) for tm in output_maps}

    while True:
        for tp in train_paths:
            [logging.info(f"Got first error: {k}") for k in stats if 'Error' in k and stats[k] == 1]
            try:
                with h5py.File(tp, 'r') as hd5:
                    dependents = {}
                    for tm in input_maps:
                        in_batch[tm.input_name()][stats['batch_index']] = tm.tensor_from_file(hd5, dependents)
                    
                    for tm in output_maps:
                        if tm in dependents:
                            out_batch[tm.output_name()][stats['batch_index']] = dependents[tm]
                        else:
                            out_batch[tm.output_name()][stats['batch_index']] = tm.tensor_from_file(hd5)

                    paths_in_batch.append(tp)
                    stats['batch_index'] += 1
                    stats['Tensors presented'] += 1
                    if stats['batch_index'] == batch_size:
                        if keep_paths:
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

        stats['epochs'] += 1
        np.random.shuffle(train_paths)
        for k in stats:
            logging.info("{}: {}".format(k, stats[k]))
        logging.info("Generator looped over and shuffled {} tensors.".format(len(train_paths)))


def multimodal_multitask_weighted_generator(batch_size, input_maps, output_maps, paths_lists, weights, keep_paths):
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

    Yields:
        A tuple of dicts for the tensor mapping tensor names to numpy arrays
        {input_1:input_tensor1, input_2:input_tensor2, ...}, {output_name1:output_tensor1, ...}
        if include_paths_in_batch is True the 3rd element of tuple will be the list of paths

    Returns:
        Never!
    """ 
    assert(len(paths_lists) > 0 and len(paths_lists) == len(weights))

    stats = Counter()
    paths_in_batch = []
    in_batch = {tm.input_name(): np.zeros((batch_size,)+tm.shape) for tm in input_maps}
    out_batch = {tm.output_name(): np.zeros((batch_size,)+tm.shape) for tm in output_maps}
    samples = [int(w*batch_size) for w in weights]

    while True:
        for tlist, num_samples in zip(paths_lists, samples):
            for tp in np.random.choice(tlist, num_samples):
                try:
                    with h5py.File(tp, 'r') as hd5:
                        dependents = {}
                        for tm in input_maps:
                            in_batch[tm.input_name()][stats['batch_index']] = tm.tensor_from_file(hd5, dependents)
                        
                        for tm in output_maps:
                            if tm in dependents:
                                out_batch[tm.output_name()][stats['batch_index']] = dependents[tm]
                            else:
                                out_batch[tm.output_name()][stats['batch_index']] = tm.tensor_from_file(hd5)

                        paths_in_batch.append(tp)
                        stats['batch_index'] += 1
                        stats['Tensors presented'] += 1
                        if stats['batch_index'] == batch_size:
                            if keep_paths:
                                yield in_batch, out_batch, paths_in_batch
                            else:
                                yield in_batch, out_batch
                            stats['batch_index'] = 0
                            paths_in_batch = []
                            for i, num_samples in enumerate(samples):
                                stats['train_paths_'+str(i)] += num_samples

                except IndexError as e:
                    stats['IndexError:'+str(e)] += 1
                except KeyError as e:
                    stats['KeyError:'+str(e)] += 1
                except OSError as e:
                    stats['OSError:'+str(e)] += 1
                except ValueError as e:
                    stats['ValueError:'+str(e)] += 1
        
        for i, tlist in enumerate(paths_lists):
            if len(tlist) <= stats['train_paths_'+str(i)]:
                stats['epochs_list_number_'+str(i)] += 1
                stats['train_paths_'+str(i)] = 0
                if len(tlist) > 1000 or stats['epochs_list_number_'+str(i)] % 5 == 0:
                    logging.info("Generator looped over {} tensors from ICD group: {}".format(len(tlist), i))
                    for k in stats:
                        logging.info('{} has: {}'.format(k, stats[k]))


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


def get_test_train_valid_paths(tensors, valid_ratio, test_ratio):
    """Return 3 disjoint lists of tensor paths.

    The paths are split in training, validation and testing lists
    apportioned according to valid_ratio and test_ratio

    Arguments:
        tensors: directory containing tensors
        valid_ratio: rate of tensors in validation list
        test_ratio: rate of tensors in testing list

    Returns:
        A tuple of 3 lists of hd5 tensor file paths
    """     
    test_paths = []
    train_paths = []
    valid_paths = []

    assert(valid_ratio > 0 and test_ratio > 0 and valid_ratio+test_ratio < 1.0)

    for root, dirs, files in os.walk(tensors):
        for name in files:
            if os.path.splitext(name)[-1].lower() != TENSOR_EXT:
                continue
            dice = np.random.rand()
            if dice < valid_ratio:
                valid_paths.append(os.path.join(root, name))
            elif dice < (valid_ratio+test_ratio):
                test_paths.append(os.path.join(root, name))
            else:   
                train_paths.append(os.path.join(root, name))
    
    if len(train_paths) == 0 or len(valid_paths) == 0 or len(test_paths) == 0:
        my_error = f"Not enough tensors at {tensors}\n"
        my_error += f"Found {len(train_paths)} training, {len(valid_paths)} validation, and {len(test_paths)} testing."
        raise ValueError(my_error)

    return train_paths, valid_paths, test_paths


def get_test_train_valid_paths_split_by_icds(tensors, icd_csv, icds, valid_ratio, test_ratio):
    lol = list(csv.reader(open(icd_csv, 'r'), delimiter='\t'))
    print(list(enumerate(lol[0])))
    stats = Counter()
    sample2group = {}
    for row in lol[1:]:
        sample_id = row[0]
        sample2group[sample_id] = 0
        for i,icd_index in enumerate(icds):
            if row[icd_index] == '1':
                sample2group[sample_id] = i+1 # group 0 means no ICD code
                stats['group_'+str(i+1)] += 1
    print(stats)
    
    test_paths = [[] for _ in range(len(icds)+1)] 
    train_paths = [[] for _ in range(len(icds)+1)] 
    valid_paths = [[] for _ in range(len(icds)+1)]                
    for root, dirs, files in os.walk(tensors):
        for name in files:
            splitted = os.path.splitext(name)
            if splitted[-1].lower() != TENSOR_EXT:
                continue
                
            sample_id = os.path.basename(splitted[0])
            group = 0
            if sample_id in sample2group:
                group = sample2group[sample_id]
            dice = np.random.rand()
            if dice < valid_ratio:
                valid_paths[group].append(os.path.join(root, name))
            elif dice < (valid_ratio+test_ratio):
                test_paths[group].append(os.path.join(root, name))
            else:   
                train_paths[group].append(os.path.join(root, name))

    logging.info('Found {} training {} validation and {} testing tensors without ICD codes.'.format(
        len(train_paths[0]), len(valid_paths[0]), len(test_paths[0])))
    for i, icd_group in enumerate(icds):
        logging.info('For ICD indexes:{} Found {} for training {} for validation and {} for testing'.format(
            icd_group, len(train_paths[i+1]), len(valid_paths[i+1]), len(test_paths[i+1])))
    
    return train_paths, valid_paths, test_paths


def test_train_valid_tensor_generators(maps_in: List[TensorMap],
                                       maps_out: List[TensorMap],
                                       tensors: str,
                                       batch_size: int,
                                       valid_ratio: float,
                                       test_ratio: float,
                                       icd_csv: str,
                                       balance_by_icds: List[str],
                                       keep_paths: bool=False,
                                       keep_paths_test: bool=True):
    if len(balance_by_icds) > 0:
        train_paths, valid_paths, test_paths = get_test_train_valid_paths_split_by_icds(tensors, icd_csv, balance_by_icds, valid_ratio, test_ratio)
        weights = [1.0/len(balance_by_icds) for _ in range(len(balance_by_icds)+1)]
        generate_train = TensorGenerator(batch_size, maps_in, maps_out, train_paths, weights, keep_paths)
        generate_valid = TensorGenerator(batch_size, maps_in, maps_out, valid_paths, weights, keep_paths)
        generate_test = TensorGenerator(batch_size, maps_in, maps_out, test_paths, weights, keep_paths or keep_paths_test)
    else:
        train_paths, valid_paths, test_paths = get_test_train_valid_paths(tensors, valid_ratio, test_ratio)
        generate_train = TensorGenerator(batch_size, maps_in, maps_out, train_paths, None, keep_paths)
        generate_valid = TensorGenerator(batch_size, maps_in, maps_out, valid_paths, None, keep_paths)
        generate_test = TensorGenerator(batch_size, maps_in, maps_out, test_paths, None, keep_paths or keep_paths_test)
    return generate_train, generate_valid, generate_test


