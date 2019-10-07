# hyperparameters.py

# Imports
import os
import logging
import numpy as np
from collections import Counter
from timeit import default_timer as timer

import hyperopt
from hyperopt import fmin, tpe, hp

import matplotlib
matplotlib.use('Agg') # Need this to write images from the GSA servers.  Order matters:
import matplotlib.pyplot as plt # First import matplotlib, then use Agg, then import plt

import keras.backend as K

from ml4cvd.defines import IMAGE_EXT
from ml4cvd.arguments import parse_args
from ml4cvd.tensor_maps_by_script import TMAPS
from ml4cvd.models import train_model_from_generators, make_multimodal_multitask_model
from ml4cvd.tensor_generators import test_train_valid_tensor_generators, big_batch_from_minibatch_generator

MAX_LOSS = 9e9


def run(args):
    try:
        # Keep track of elapsed execution time
        start_time = timer()

        if 'conv' == args.mode:
            optimize_conv_layers_multimodal_multitask(args)
        elif 'dense_layers' == args.mode:
            optimize_dense_layers_multimodal_multitask(args)
        elif 'lr' == args.mode:
            optimize_lr_multimodal_multitask(args)
        elif 'inputs' == args.mode:
            optimize_input_tensor_maps(args)
        else:
            raise ValueError('Unknown hyperparameter optimization mode:', args.mode)
  
    except Exception as e:
        logging.exception(e)
        
    end_time = timer()
    elapsed_time = end_time - start_time
    logging.info("Executed the '{}' operation in {:.2f} seconds".format(args.mode, elapsed_time))


def optimize_conv_layers_multimodal_multitask(args):
    stats = Counter()
    generate_train, _, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors, args.batch_size,
                                                                          args.valid_ratio, args.test_ratio, args.test_modulo, args.balance_csvs, False, False)
    test_data, test_labels = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps, False)

    dense_blocks_sets = [[16], [32], [48], [32, 16], [32, 32], [32, 24, 16], [48, 32, 24], [48, 48, 48]]
    conv_layers_sets = [[64], [48], [32], [24]]
    dense_layers_sets = [[16, 64], [8, 128], [48], [32], [24], [16]]
    pool_zs = [1]
    param_lists = {'conv_layers': conv_layers_sets, 'dense_blocks': dense_blocks_sets, 'dense_layers': dense_layers_sets}
    space = {
        'pool_x': hp.loguniform('pool_x', 0, 4),
        'conv_layers': hp.choice('conv_layers', conv_layers_sets),
        'dense_blocks': hp.choice('dense_blocks', dense_blocks_sets),      
        'dense_layers': hp.choice('dense_layers', dense_layers_sets),
    }

    def loss_from_multimodal_multitask(x):
        try:
            x['pool_x'] = int(x['pool_x'])
            set_args_from_x(args, x)
            model = make_multimodal_multitask_model(**args.__dict__)
            
            if model.count_params() > args.max_parameters:
                logging.info(f"Model too big, max parameters is:{args.max_parameters}, model has:{model.count_params()}. Return max loss.")
                del model
                return MAX_LOSS

            model = train_model_from_generators(model, generate_train, generate_test, args.training_steps, args.validation_steps, args.batch_size,
                                                args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)
            loss_and_metrics = model.evaluate(test_data, test_labels, batch_size=args.batch_size)
            stats['count'] += 1
            logging.info('Current architecture: {}'.format(string_from_arch_dict(x)))
            logging.info('Iteration {} out of maximum {}: Loss: {} Current model size: {}.'.format(stats['count'], args.max_models, loss_and_metrics[0], model.count_params()))
            del model
            return loss_and_metrics[0]
        
        except ValueError as e:
            logging.exception('ValueError trying to make a model for hyperparameter optimization. Returning max loss.')
            return MAX_LOSS
        except:
            logging.exception('Error trying hyperparameter optimization. Returning max loss.')
            return MAX_LOSS

    trials = hyperopt.Trials()
    fmin(loss_from_multimodal_multitask, space=space, algo=tpe.suggest, max_evals=args.max_models, trials=trials)
    plot_trials(trials, os.path.join(args.output_folder, args.id, 'loss_per_iteration'+IMAGE_EXT), param_lists)

    args = args_from_best_trials(args, trials, param_lists)
    _ = make_multimodal_multitask_model(**args.__dict__)


def optimize_dense_layers_multimodal_multitask(args):
    stats = Counter()
    generate_train, _, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors, args.batch_size,
                                                                          args.valid_ratio, args.test_ratio, args.test_modulo, args.icd_csv,
                                                                          args.balance_by_icds, False, False)
    test_data, test_labels = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps, False)
    space = {'num_layers': hp.uniform('num_layers', 1, 6),
             'layer_width': hp.loguniform('layer_width', 2, 7)}

    def loss_from_multimodal_multitask(x):
        try:
            args.dense_layers = [int(x['layer_width'])] * int(x['num_layers'])
            model = make_multimodal_multitask_model(**args.__dict__)

            if model.count_params() > args.max_parameters:
                logging.info(f"Model too big, max parameters is:{args.max_parameters}, model has:{model.count_params()}. Return max loss.")
                del model
                return MAX_LOSS

            model = train_model_from_generators(model, generate_train, generate_test, args.training_steps, args.validation_steps, args.batch_size,
                                                args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)
            loss_and_metrics = model.evaluate(test_data, test_labels, batch_size=args.batch_size)
            stats['count'] += 1
            logging.info('Current architecture: {}'.format(string_from_arch_dict(x)))
            logging.info('Iteration {} out of maximum {}: Loss: {} Current model size: {}.'.format(stats['count'], args.max_models, loss_and_metrics[0], model.count_params()))
            del model
            return loss_and_metrics[0]

        except ValueError as e:
            logging.exception('ValueError trying to make a model for hyperparameter optimization. Returning max loss.')
            return MAX_LOSS
        except:
            logging.exception('Error trying hyperparameter optimization. Returning max loss.')
            return MAX_LOSS

    trials = hyperopt.Trials()
    fmin(loss_from_multimodal_multitask, space=space, algo=tpe.suggest, max_evals=args.max_models, trials=trials)
    plot_trials(trials, os.path.join(args.output_folder, args.id, 'loss_per_iteration'+IMAGE_EXT))

    # Re-train the best model so it's easy to view it at the end of the logs
    args = args_from_best_trials(args, trials)
    _ = make_multimodal_multitask_model(**args.__dict__)


def optimize_lr_multimodal_multitask(args):
    stats = Counter()
    generate_train, _, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors, args.batch_size,
                                                                          args.valid_ratio, args.test_ratio, args.test_modulo, args.icd_csv,
                                                                          args.balance_by_icds, False, False)
    test_data, test_labels = big_batch_from_minibatch_generator(args.tensor_maps_in, args.tensor_maps_out, generate_test, args.test_steps, False)

    space = {'learning_rate': hp.loguniform('learning_rate', -10, -2)}

    def loss_from_multimodal_multitask(x):
        try:
            set_args_from_x(args, x)
            model = make_multimodal_multitask_model(**args.__dict__)
            if model.count_params() > args.max_parameters:
                logging.info(f"Model too big, max parameters is:{args.max_parameters}, model has:{model.count_params()}. Return max loss.")
                return MAX_LOSS
            
            logging.info('Current parameter set: {} \n'.format(string_from_arch_dict(x)))
            model = train_model_from_generators(model, generate_train, generate_test, args.training_steps, args.validation_steps, args.batch_size,
                                                args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)
            loss_and_metrics = model.evaluate(test_data, test_labels, batch_size=args.batch_size)
            stats['count'] += 1
            logging.info('Iteration {} out of maximum {}. Loss: {} Current model size {}.'.format(stats['count'], args.max_models, loss_and_metrics[0], model.count_params()))
            return loss_and_metrics[0]
        
        except ValueError as e:
            logging.exception('ValueError trying to make a model for hyperparameter optimization. Returning max loss.')
            return MAX_LOSS
    
    trials = hyperopt.Trials()
    fmin(loss_from_multimodal_multitask, space=space, algo=tpe.suggest, max_evals=args.max_models, trials=trials)
    plot_trials(trials, os.path.join(args.output_folder, args.id, 'loss_per_iteration'+IMAGE_EXT))

    # Rebuild the best model so it's easy to view it at the end of the logs
    args = args_from_best_trials(args, trials)
    _ = make_multimodal_multitask_model(**args.__dict__)


def optimize_input_tensor_maps(args):
    stats = Counter()
    generate_train, _, generate_test = test_train_valid_tensor_generators(args.tensor_maps_in, args.tensor_maps_out, args.tensors, args.batch_size,
                                                                          args.valid_ratio, args.test_ratio, args.test_modulo, args.icd_csv,
                                                                          args.balance_by_icds, False, False)
    input_tensor_map_sets = [['categorical-phenotypes-72'], ['mri-slice'], ['sax_inlinevf_zoom'], ['cine_segmented_sax_inlinevf'], ['ekg-leads']]
    param_lists = {'input_tensor_maps': input_tensor_map_sets}
    space = {'input_tensor_maps': hp.choice('input_tensor_maps', input_tensor_map_sets),}

    def loss_from_multimodal_multitask(x):
        try:
            set_args_from_x(args, x)
            model = make_multimodal_multitask_model(**args.__dict__)
            if model.count_params() > args.max_parameters:
                logging.info(f"Model too big, max parameters is:{args.max_parameters}, model has:{model.count_params()}. Return max loss.")
                return MAX_LOSS
            
            logging.info('Current parameter set: {} \n'.format(string_from_arch_dict(x)))
            model = train_model_from_generators(model, generate_train, generate_test, args.training_steps, args.validation_steps, args.batch_size,
                                                args.epochs, args.patience, args.output_folder, args.id, args.inspect_model, args.inspect_show_labels)
            loss_and_metrics = model.evaluate_generator(generate_test, steps=args.test_steps)
            stats['count'] += 1
            logging.info('Iteration {} out of maximum {}: Loss: {} Current model size: {}.'.format(stats['count'], args.max_models, loss_and_metrics[0], model.count_params()))
            return loss_and_metrics[0]
        
        except ValueError as e:
            logging.exception('ValueError trying to make a model for hyperparameter optimization. Returning max loss.')
            return MAX_LOSS
    
    trials = hyperopt.Trials()
    fmin(loss_from_multimodal_multitask, space=space, algo=tpe.suggest, max_evals=args.max_models, trials=trials)
    plot_trials(trials, os.path.join(args.output_folder, args.id, 'loss_per_iteration'+IMAGE_EXT), param_lists)

    # Rebuild the best model so it's easy to view it at the end of the logs
    args = args_from_best_trials(args, trials, param_lists)
    _ = make_multimodal_multitask_model(**args.__dict__)


def set_args_from_x(args, x):
    for k in args.__dict__:
        if k in x:
            if isinstance(args.__dict__[k], int):
                args.__dict__[k] = int(x[k])
            elif isinstance(args.__dict__[k], float):
                args.__dict__[k] = float(x[k])
            else:
                args.__dict__[k] = x[k]
    logging.info(f"Set arguments to: {args}")
    args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
    args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]


def string_from_arch_dict(x):
    s = ''
    for k in x:
        s += '\n' + k + ' = '
        s += str(x[k])     
    return s


def args_from_best_trials(args, trials, param_lists={}):
    best_trial_idx = np.argmin(trials.losses())
    x = trials.trials[best_trial_idx]['misc']['vals']
    logging.info(f"got best x {x} best model is:{string_from_trials(trials, best_trial_idx, param_lists)}")
    for k in x:
        v = x[k][0]
        if k in param_lists:
            args.__dict__[k] = param_lists[k][int(v)]
        elif k in ['conv_x', 'conv_y', 'conv_z']:
            args.__dict__[k] = int(v)
        else:
            args.__dict__[k] = v
    args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
    args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
    return args   


def string_from_trials(trials, index, param_lists={}):
    s = ''
    x = trials.trials[index]['misc']['vals']
    for k in x:
        s += '\n' + k + ' = '
        v = x[k][0]
        if k in param_lists:
            s += str(param_lists[k][int(v)])
        elif k in ['num_layers', 'layer_width']:
            s += str(int(v))
        else:
            s += str(v)
    return s


def plot_trials(trials, figure_path, param_lists={}):
    lmax = max([x for x in trials.losses() if x != MAX_LOSS]) + 1 # add to the max to distinguish real losses from max loss
    lplot = [x if x != MAX_LOSS else lmax for x in trials.losses()]
    best_loss = min(lplot)
    worst_loss = max(lplot)
    std = np.std(lplot)
    plt.figure(figsize=(64, 64))
    matplotlib.rcParams.update({'font.size': 9})
    plt.plot(lplot)
    for i in range(len(trials.trials)):
        if best_loss+std > lplot[i]:
            plt.text(i, lplot[i], string_from_trials(trials, i, param_lists))
        elif worst_loss-std < lplot[i]:
            plt.text(i, lplot[i], string_from_trials(trials, i, param_lists))

    plt.xlabel('Iterations')
    plt.ylabel('Losses')
    plt.title('Hyperparameter Optimization\n')
    if not os.path.exists(os.path.dirname(figure_path)):
        os.makedirs(os.path.dirname(figure_path))
    plt.savefig(figure_path)
    logging.info('Saved loss plot to:{}'.format(figure_path))


def limit_mem():
    try:
        K.clear_session()
        cfg = K.tf.ConfigProto()
        cfg.gpu_options.allow_growth = True
        K.set_session(K.tf.Session(config=cfg))
    except AttributeError as e:
        logging.exception('Could not clear session. Maybe you are using Theano backend?')


if __name__=='__main__':
    args = parse_args()
    run(args)  # back to the top