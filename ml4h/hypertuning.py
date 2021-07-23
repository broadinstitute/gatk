# hypertuning.py

# Imports
import logging
from timeit import default_timer as timer

import numpy as np
import tensorflow as tf
import kerastuner as kt
import kerastuner.engine.trial as trial_module
from tensorflow.keras.callbacks import EarlyStopping
from kerastuner.tuners import RandomSearch, BayesianOptimization, Hyperband

from ml4h.arguments import parse_args
from ml4h.models.layer_wrappers import NORMALIZATION_CLASSES
from ml4h.models.model_factory import block_make_multimodal_multitask_model
from ml4h.tensor_generators import test_train_valid_tensor_generators

tuner_type = 'bayes'
MAX_MODEL_SIZE = 15000000
MAX_LOSS = 9e9

def run(args):
    start_time = timer()
    if 'conv' == args.mode:
        model_builder = make_model_builder(args)
    elif 'activation' == args.mode:
        model_builder = make_model_builder_activation(args)
    elif 'normalization' == args.mode:
        model_builder = make_model_builder_normalization(args)
    else:
        raise ValueError('Unknown hyper-parameter optimization mode:', args.mode)

    if 'random' == tuner_type:
        tuner = RandomSearch(
            model_builder,
            objective='val_loss', # kt.Objective("val_pearson", direction="max"),
            max_trials=args.max_models,
            executions_per_trial=args.min_samples,
            directory=args.output_folder,
            project_name=args.id,
            seed=args.random_seed,
        )
    elif 'bayes' == tuner_type:
        tuner = BayesianSearchEdit(
            model_builder,
            objective='val_loss', #kt.Objective("val_pearson", direction="max"),
            max_trials=args.max_models,
            #max_model_size=args.max_parameters,
            executions_per_trial=args.min_samples,
            directory=args.output_folder,
            project_name=args.id,
            seed=args.random_seed,
            beta=2.3,  # Explore exploit tradeoff, higher value mean more exploration
        )
    generate_train, generate_valid, generate_test = test_train_valid_tensor_generators(**args.__dict__)
    stop_early = EarlyStopping(monitor='val_loss', patience=args.patience)
    tuner.search(generate_train, epochs=args.epochs, steps_per_epoch=args.training_steps,
                 validation_data=generate_valid, validation_steps=args.validation_steps, callbacks=[stop_early])
    logging.info(f"Tuning done best models below!")
    end_time = timer()
    tuner.search_space_summary()
    tuner.results_summary()
    for i, best_hyper in enumerate(tuner.get_best_hyperparameters(num_trials=8)):
        logging.info(f"\n\n#{i+1} best hyperparameters:\n{best_hyper.values}")
    logging.info(f"\nExecuted {args.mode} mode in {(end_time - start_time) / 60.0:.1f} minutes.")


def make_model_builder(args):
    model_count = 0

    def model_builder(hp):
        num_conv_layers = hp.Int('num_conv_layers', 0, 3)
        conv_layer_size = hp.Int('conv_layer_size', 8, 64, step=8)
        args.__dict__['conv_layers'] = [conv_layer_size] * num_conv_layers
        num_dense_blocks = hp.Int('num_dense_blocks', 1, 3)
        dense_block_size = hp.Int('dense_block_size', 8, 48, step=8)
        args.__dict__['dense_blocks'] = [dense_block_size] * num_dense_blocks
        args.__dict__['block_size'] = hp.Int('block_size', 1, 6)
        # num_dense_layers = hp.Int('num_dense_layers', 1, 4)
        # dense_layer_size = hp.Int('dense_layer_size', 16, 256, sampling='log')
        # args.__dict__['dense_layers'] = [dense_layer_size] * num_dense_layers
        args.__dict__['activation'] = hp.Choice('activation', ['leaky', 'swish', 'gelu', 'lisht', 'mish', 'relu', 'selu'])
        dense_normalize = hp.Choice('dense_normalize', list(NORMALIZATION_CLASSES.keys()) + ['None'])
        args.__dict__['dense_normalize'] = None if dense_normalize == 'None' else dense_normalize
        conv_normalize = hp.Choice('conv_normalize', list(NORMALIZATION_CLASSES.keys()) + ['None'])
        args.__dict__['conv_normalize'] = None if conv_normalize == 'None' else conv_normalize
        args.__dict__['pool_type'] = 'max' if hp.Boolean('pool_type_is_max') else 'average'
        model, _, _, _ = block_make_multimodal_multitask_model(**args.__dict__)
        nonlocal model_count
        logging.info(f'Hyper-tuner is {100.0*(model_count / (args.max_models*args.min_samples)):0.1f}% complete. '
                     f'Built model #{model_count} with {model.count_params()} parameters of a maximum of {args.max_models*args.min_samples}.')
        model_count += 1
        return model
    return model_builder


def make_model_builder_activation(args):
    def model_builder(hp):
        args.__dict__['activation'] = hp.Choice('activation', ['leaky', 'swish', 'gelu', 'lisht', 'mish', 'relu', 'selu'])
        model, _, _, _ = block_make_multimodal_multitask_model(**args.__dict__)
        return model
    return model_builder


def make_model_builder_normalization(args):
    def model_builder(hp):
        dense_normalize = hp.Choice('dense_normalize', list(NORMALIZATION_CLASSES.keys()) + ['None'])
        args.__dict__['dense_normalize'] = None if dense_normalize == 'None' else dense_normalize
        conv_normalize = hp.Choice('conv_normalize', list(NORMALIZATION_CLASSES.keys()) + ['None'])
        args.__dict__['conv_normalize'] = None if conv_normalize == 'None' else conv_normalize
        model, _, _, _ = block_make_multimodal_multitask_model(**args.__dict__)
        return model
    return model_builder


class BayesianSearchEdit(BayesianOptimization):
    """
    TO-DO: add custom max_model_size input param to class
    def __init__(self):
        pass
    """

    def on_trial_end(self, trial):
        """A hook called after each trial is run.
        # Arguments:
            trial: A `Trial` instance.
        """
        # Send status to Logger
        if self.logger:
            self.logger.report_trial_state(trial.trial_id, trial.get_state())

        if not trial.get_state().get("status") == trial_module.TrialStatus.INVALID:
            self.oracle.end_trial(trial.trial_id, trial_module.TrialStatus.COMPLETED)

        self.oracle.update_space(trial.hyperparameters)
        # Display needs the updated trial scored by the Oracle.
        self._display.on_trial_end(self.oracle.get_trial(trial.trial_id))
        self.save()

    def _build_and_fit_model(self, trial, fit_args, fit_kwargs):
        model = self.hypermodel.build(trial.hyperparameters)
        model_size = self.maybe_compute_model_size(model)
        print("Considering model with size: {}".format(model_size))

        if model_size > MAX_MODEL_SIZE:
            #self.oracle.end_trial(trial.trial_id, trial_module.TrialStatus.INVALID)

            dummy_history_obj = tf.keras.callbacks.History()
            dummy_history_obj.on_train_begin()
            dummy_history_obj.history.setdefault('val_loss', []).append(MAX_LOSS)
            return dummy_history_obj

        return model.fit(*fit_args, **fit_kwargs)

    def maybe_compute_model_size(self, model):
        """Compute the size of a given model, if it has been built."""
        if model.built:
            params = [tf.keras.backend.count_params(p) for p in model.trainable_weights]
            return int(np.sum(params))
        return 0


if __name__ == '__main__':
    args = parse_args()
    run(args)  # back to the top
