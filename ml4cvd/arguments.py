# arguments.py
#
# Command Line Arguments for Machine Learning 4 CardioVascular Disease
# Shared by recipes.py and other command-line runnable files.
# These arguments are a bit of a hodge-podge and are used promiscuously throughout these files.
# Sometimes code overwrites user-provided arguments to enforce assumptions or sanity.
#
# October 2018
# Sam Friedman 
# sam@broadinstitute.org

# Imports
import os
import sys
import logging
import argparse
import operator
import datetime
import numpy as np
import multiprocessing

from ml4cvd.logger import load_config
from ml4cvd.TensorMap import TensorMap
from ml4cvd.tensor_maps_by_hand import TMAPS
from ml4cvd.defines import IMPUTATION_RANDOM, IMPUTATION_MEAN
from ml4cvd.tensor_map_maker import generate_multi_field_continuous_tensor_map, generate_continuous_tensor_map_from_file


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mode', default='mlp', help='What would you like to do?')

    # Config arguments
    parser.add_argument("--logging_level", default='INFO', choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level. Overrides any configuration given in the logging configuration file.")

    # Tensor Map arguments
    parser.add_argument('--input_tensors', default=[], nargs='+')
    parser.add_argument('--output_tensors', default=[], nargs='+')
    parser.add_argument('--input_continuous_tensors', default=[], nargs='+', help='Continuous tensor maps to be combined.')
    parser.add_argument('--tensor_maps_in', default=[], help='Do not set this directly. Use input_tensors')
    parser.add_argument('--tensor_maps_out', default=[], help='Do not set this directly. Use output_tensors')

    # Input and Output files and directories
    parser.add_argument('--bigquery_credentials_file', default='/mnt/ml4cvd/projects/jamesp/bigquery/bigquery-viewer-credentials.json',
                        help='Path to service account credentials for looking up BigQuery tables.')
    parser.add_argument('--bigquery_dataset', default='broad-ml4cvd.ukbb7089_r10data', help='BigQuery dataset containing tables we want to query.')
    parser.add_argument('--xml_folder', default='/mnt/disks/ecg-rest-xml/', help='Path to folder of XMLs of ECG data.')
    parser.add_argument('--zip_folder', default='/mnt/disks/sax-mri-zip/', help='Path to folder of zipped dicom images.')
    parser.add_argument('--phenos_folder', default='gs://ml4cvd/phenotypes/', help='Path to folder of phenotype defining CSVs.')
    parser.add_argument('--phecode_definitions', default='/mnt/ml4cvd/projects/jamesp/data/phecode_definitions1.2.csv', help='CSV of phecode definitions')
    parser.add_argument('--dicoms', default='./dicoms/', help='Path to folder of dicoms.')
    parser.add_argument('--test_csv', default=None, help='Path to CSV with Sample IDs to reserve for testing')
    parser.add_argument('--app_csv', help='Path to file used to link sample IDs between UKBB applications 17488 and 7089')
    parser.add_argument('--tensors', help='Path to folder containing tensors, or where tensors will be written.')
    parser.add_argument('--output_folder', default='./recipes_output/', help='Path to output folder for recipes.py runs.')
    parser.add_argument('--model_file', help='Path to a saved model architecture and weights (hd5).')
    parser.add_argument('--model_files', nargs='*', default=[], help='List of paths to saved model architectures and weights (hd5).')
    parser.add_argument('--model_layers', help='Path to a model file (hd5) which will be loaded by layer, useful for transfer learning.')
    parser.add_argument('--freeze_model_layers', default=False, action='store_true', help='Whether to freeze the layers from model_layers.')
    parser.add_argument('--continuous_file', default=None, help='Path to a file containing continuous values from which a output TensorMap will be made.'
                                                               'Note that setting this argument has the effect of linking the first output_tensors'
                                                               'argument to the TensorMap made from this file.')

    # Data selection parameters
    parser.add_argument('--continuous_file_column', default=None, help='Column header in file from which a continuous TensorMap will be made.')
    parser.add_argument('--continuous_file_normalize', default=False, action='store_true', help='Whether to normalize a continuous TensorMap made from a file.')
    parser.add_argument('--categorical_field_ids', nargs='*', default=[], type=int,
        help='List of field ids from which input features will be collected.')
    parser.add_argument('--continuous_field_ids', nargs='*', default=[], type=int,
        help='List of field ids from which continuous real-valued input features will be collected.')
    parser.add_argument('--include_array', default=False, action='store_true', help='Include array idx for UKBB phenotypes.')
    parser.add_argument('--include_instance', default=False, action='store_true', help='Include instances for UKBB phenotypes.')
    parser.add_argument('--min_values', default=10, type=int, help='Per feature size minimum.')
    parser.add_argument('--min_samples', default=3, type=int, help='Min number of samples to require for calculating correlations.')
    parser.add_argument('--max_samples', type=int, default=None,
                        help='Max number of samples to use for tensor reporting -- all samples are used if not specified.')
    parser.add_argument('--mri_field_ids', default=['20208', '20209'], nargs='*', help='Field id for MR images.')
    parser.add_argument('--xml_field_ids', default=['20205', '6025'], nargs='*', help='Field id for XMLs of resting and exercise ECG data.')
    parser.add_argument('--max_patients', default=999999, type=int,  help='Maximum number of patient data to read')
    parser.add_argument('--min_sample_id', default=0, type=int, help='Minimum sample id to write to tensor.')
    parser.add_argument('--max_sample_id', default=7000000, type=int, help='Maximum sample id to write to tensor.')
    parser.add_argument('--max_slices', default=999999, type=int, help='Maximum number of dicom slices to read')
    parser.add_argument('--dicom_series', default='cine_segmented_sax_b6', help='Maximum number of dicom slices to read')
    parser.add_argument('--b_slice_force', default=None,
                        help='If set, will only load specific b slice for short axis MRI diastole systole tensor maps (i.e b0, b1, b2, ... b10).')
    parser.add_argument('--include_missing_continuous_channel', default=False, action='store_true',
                        help='Include missing channels in continuous tensors')
    parser.add_argument('--imputation_method_for_continuous_fields', default=IMPUTATION_RANDOM, help='can be random or mean',
                        choices=[IMPUTATION_RANDOM, IMPUTATION_MEAN])

    # Model Architecture Parameters
    parser.add_argument('--x', default=256, type=int, help='x tensor resolution')
    parser.add_argument('--y', default=256, type=int, help='y tensor resolution')
    parser.add_argument('--zoom_x', default=50, type=int, help='zoom_x tensor resolution')
    parser.add_argument('--zoom_y', default=35, type=int, help='zoom_y tensor resolution')
    parser.add_argument('--zoom_width', default=96, type=int, help='zoom_width tensor resolution')
    parser.add_argument('--zoom_height', default=96, type=int, help='zoom_height tensor resolution')
    parser.add_argument('--z', default=48, type=int, help='z tensor resolution')
    parser.add_argument('--t', default=48, type=int, help='Number of time slices')
    parser.add_argument('--mlp_concat', default=False, action='store_true', help='Concatenate input with every multiplayer perceptron layer.')
    parser.add_argument('--dense_layers', nargs='*', default=[16, 64], type=int, help='List of number of hidden units in neural nets dense layers.')
    parser.add_argument('--dropout', default=0.0, type=float, help='Dropout rate of dense layers must be in [0.0, 1.0].')
    parser.add_argument('--activation', default='relu',  help='Activation function for hidden units in neural nets dense layers.')
    parser.add_argument('--conv_layers', nargs='*', default=[32], type=int, help='List of number of kernels in convolutional layers.')
    parser.add_argument('--conv_x', default=3, type=int, help='X dimension of convolutional kernel.')
    parser.add_argument('--conv_y', default=3, type=int, help='Y dimension of convolutional kernel.')
    parser.add_argument('--conv_z', default=2, type=int, help='Z dimension of convolutional kernel.')
    parser.add_argument('--conv_width', default=71, type=int, help='Width of convolutional kernel for 1D CNNs.')
    parser.add_argument('--conv_dilate', default=False, action='store_true', help='Dilate the convolutional layers.')
    parser.add_argument('--conv_dropout', default=0.0, type=float, help='Dropout rate of convolutional kernels must be in [0.0, 1.0].')
    parser.add_argument('--conv_type', default='conv', choices=['conv', 'separable', 'depth'], help='Type of convolutional layer')
    parser.add_argument('--conv_normalize', default=None, choices=['', 'batch_norm'], help='Type of normalization layer for convolutions')
    parser.add_argument('--conv_regularize', default=None, choices=['dropout', 'spatial_dropout'], help='Type of regularization layer for convolutions.')
    parser.add_argument('--max_pools', nargs='*', default=[], type=int, help='List of maxpooling layers.')
    parser.add_argument('--pool_type', default='max', choices=['max', 'average'], help='Type of pooling layers.')
    parser.add_argument('--pool_x', default=2, type=int, help='Pooling size in the x-axis, if 1 no pooling will be performed.')
    parser.add_argument('--pool_y', default=2, type=int, help='Pooling size in the y-axis, if 1 no pooling will be performed.')
    parser.add_argument('--pool_z', default=1, type=int, help='Pooling size in the z-axis, if 1 no pooling will be performed.')
    parser.add_argument('--res_layers', nargs='*', default=[], type=int, help='List of residual layers.')
    parser.add_argument('--padding', default='same', help='Valid or same border padding on the convolutional layers.')
    parser.add_argument('--dense_blocks', nargs='*', default=[32, 24, 16], type=int, help='List of number of kernels in convolutional layers.')
    parser.add_argument('--block_size', default=3, type=int, help='Number of convolutional layers within a block.')
    parser.add_argument('--u_connect', default=False, action='store_true', help='Connect early convolutional layers to later ones of the same size, as in U-Net.')
    parser.add_argument('--aligned_dimension', default=16, type=int, help='Dimensionality of aligned embedded space for multi-modal alignment models.')
    parser.add_argument('--max_parameters', default=9000000, type=int,
                        help='Maximum number of trainable parameters in a model during hyperparameter optimization.')
    parser.add_argument('--hidden_layer', default='embed', help='Name of a hidden layer for inspections.')

    # Training and Hyper-Parameter Optimization Parameters
    parser.add_argument('--epochs', default=12, type=int, help='Number of training epochs.')
    parser.add_argument('--batch_size', default=16, type=int, help='Mini batch size for stochastic gradient descent algorithms.')
    parser.add_argument('--valid_ratio', default=0.2, type=float, help='Rate of training tensors to save for validation must be in [0.0, 1.0].')
    parser.add_argument('--test_ratio', default=0.1, type=float, help='Rate of training tensors to save for testing [0.0, 1.0].')
    parser.add_argument('--test_modulo', default=10, type=int,
                        help='Sample IDs modulo this number will be reserved for testing. Set to 1 to only reserve test_ratio for testing.')
    parser.add_argument('--test_steps', default=32, type=int, help='Number of batches to use for testing.')
    parser.add_argument('--training_steps', default=400, type=int, help='Number of training batches to examine in an epoch.')
    parser.add_argument('--validation_steps', default=40, type=int, help='Number of validation batches to examine in an epoch validation.')
    parser.add_argument('--learning_rate', default=0.0002, type=float, help='Learning rate during training.')
    parser.add_argument('--mixup_alpha', default=0, type=float, help='If positive apply mixup and sample from a Beta with this value as shape parameter alpha.')
    parser.add_argument('--label_weights', nargs='*', type=float,
                        help='List of per-label weights for weighted categorical cross entropy. If provided, must map 1:1 to number of labels.')
    parser.add_argument('--patience', default=8, type=int,
                        help='Early Stopping parameter: Maximum number of epochs to run without validation loss improvements.')
    parser.add_argument('--max_models', default=16, type=int,
                        help='Maximum number of models for the hyper-parameter optimizer to evaluate before returning.')
    parser.add_argument('--balance_csvs', default=[], nargs='*', help='Balances batches with representation from sample IDs in this list of CSVs')
    parser.add_argument('--optimizer', default='radam', type=str, help='Optimizer for model training')

    # Run specific and debugging arguments
    parser.add_argument('--id', default='no_id', help='Identifier for this run, user-defined string to keep experiments organized.')
    parser.add_argument('--random_seed', default=12878, type=int, help='Random seed to use throughout run.  Always use np.random.')
    parser.add_argument('--write_pngs', default=False, action='store_true', help='Write pngs of slices.')
    parser.add_argument('--debug', default=False, action='store_true', help='Run in debug mode.')
    parser.add_argument('--inspect_model', default=False, action='store_true', help='Plot model architecture, measure inference and training speeds.')
    parser.add_argument('--inspect_show_labels', default=True, action='store_true', help='Plot model architecture with labels for each layer.')
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha transparency for t-SNE plots must in [0.0-1.0].')

    # Training optimization options
    parser.add_argument('--num_workers', default=multiprocessing.cpu_count(), type=int, help="Number of workers to use for every tensor generator.")
    parser.add_argument('--cache_size', default=3.5e9/multiprocessing.cpu_count(), type=float, help="Tensor map cache size per worker.")

    args = parser.parse_args()
    _process_args(args)
    return args


def _get_tmap(name: str) -> TensorMap:
    """
    This allows tensor_maps_by_script to only be imported if necessary, because it's slow.
    """
    if name in TMAPS:
        return TMAPS[name]
    from ml4cvd.tensor_maps_by_script import TMAPS as SCRIPT_TMAPS
    TMAPS.update(SCRIPT_TMAPS)
    return TMAPS[name]


def _process_args(args):
    now_string = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    args_file = os.path.join(args.output_folder, args.id, 'arguments_' + now_string + '.txt')
    command_line = f"\n\n./scripts/tf.sh {' '.join(sys.argv)}\n\n\n"
    if not os.path.exists(os.path.dirname(args_file)):
        os.makedirs(os.path.dirname(args_file))
    with open(args_file, 'w') as f:
        f.write(command_line)
        for k, v in sorted(args.__dict__.items(), key=operator.itemgetter(0)):
            f.write(k + ' = ' + str(v) + '\n')
    load_config(args.logging_level, os.path.join(args.output_folder, args.id), 'log_' + now_string, args.min_sample_id)
    args.tensor_maps_in = [_get_tmap(it) for it in args.input_tensors]
    if len(args.input_continuous_tensors) > 0:
        multi_field_tensor_map = [generate_multi_field_continuous_tensor_map(args.input_continuous_tensors, args.include_missing_continuous_channel,
                                                                             args.imputation_method_for_continuous_fields)]
        args.tensor_maps_in.extend(multi_field_tensor_map)

    args.tensor_maps_out = []
    if args.continuous_file is not None:
        # Continuous TensorMap generated from file is given the name specified by the first output_tensors argument
        args.tensor_maps_out.append(generate_continuous_tensor_map_from_file(args.continuous_file, args.continuous_file_column,
                                                                             args.output_tensors.pop(0), args.continuous_file_normalize))
    args.tensor_maps_out.extend([_get_tmap(ot) for ot in args.output_tensors])

    np.random.seed(args.random_seed)

    logging.info(f"Command Line was:{command_line}")
    logging.info(f"Total TensorMaps:{len(TMAPS)} Arguments are {args}")
