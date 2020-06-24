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
from typing import Set, Dict, List, Optional
from collections import defaultdict

from ml4cvd.logger import load_config
from ml4cvd.TensorMap import TensorMap
from ml4cvd.models import parent_sort, BottleneckType, check_no_bottleneck
from ml4cvd.tensor_maps_by_hand import TMAPS
from ml4cvd.defines import IMPUTATION_RANDOM, IMPUTATION_MEAN
from ml4cvd.tensor_maps_partners_ecg import build_partners_tensor_maps, build_cardiac_surgery_tensor_maps, build_partners_time_series_tensor_maps
from ml4cvd.tensor_map_maker import generate_continuous_tensor_map_from_file


BOTTLENECK_STR_TO_ENUM = {
    'flatten_restructure': BottleneckType.FlattenRestructure,
    'global_average_pool': BottleneckType.GlobalAveragePoolStructured,
    'variational': BottleneckType.Variational,
    'no_bottleneck': BottleneckType.NoBottleNeck,
}


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mode', default='mlp', help='What would you like to do?')

    # Config arguments
    parser.add_argument(
        "--logging_level", default='INFO', choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level. Overrides any configuration given in the logging configuration file.",
    )

    # Tensor Map arguments
    parser.add_argument('--input_tensors', default=[], nargs='+')
    parser.add_argument('--output_tensors', default=[], nargs='+')
    parser.add_argument('--sample_weight', default=None,  help='TensorMap key for sample weight in training.')
    parser.add_argument('--tensor_maps_in', default=[], help='Do not set this directly. Use input_tensors')
    parser.add_argument('--tensor_maps_out', default=[], help='Do not set this directly. Use output_tensors')

    # Input and Output files and directories
    parser.add_argument(
        '--bigquery_credentials_file', default='/mnt/ml4cvd/projects/jamesp/bigquery/bigquery-viewer-credentials.json',
        help='Path to service account credentials for looking up BigQuery tables.',
    )
    parser.add_argument('--bigquery_dataset', default='broad-ml4cvd.ukbb7089_r10data', help='BigQuery dataset containing tables we want to query.')
    parser.add_argument('--xml_folder', default='/mnt/disks/ecg-rest-xml/', help='Path to folder of XMLs of ECG data.')
    parser.add_argument('--zip_folder', default='/mnt/disks/sax-mri-zip/', help='Path to folder of zipped dicom images.')
    parser.add_argument('--phenos_folder', default='gs://ml4cvd/phenotypes/', help='Path to folder of phenotype defining CSVs.')
    parser.add_argument('--phecode_definitions', default='/mnt/ml4cvd/projects/jamesp/data/phecode_definitions1.2.csv', help='CSV of phecode definitions')
    parser.add_argument('--dicoms', default='./dicoms/', help='Path to folder of dicoms.')
    parser.add_argument('--sample_csv', default=None, help='Path to CSV with Sample IDs to restrict tensor paths')
    parser.add_argument('--tsv_style', default='standard', choices=['standard', 'genetics'], help='Format choice for the TSV file produced in output by infer and explore modes.')
    parser.add_argument('--app_csv', help='Path to file used to link sample IDs between UKBB applications 17488 and 7089')
    parser.add_argument('--tensors', help='Path to folder containing tensors, or where tensors will be written.')
    parser.add_argument('--output_folder', default='./recipes_output/', help='Path to output folder for recipes.py runs.')
    parser.add_argument('--model_file', help='Path to a saved model architecture and weights (hd5).')
    parser.add_argument('--model_files', nargs='*', default=[], help='List of paths to saved model architectures and weights (hd5).')
    parser.add_argument('--model_layers', help='Path to a model file (hd5) which will be loaded by layer, useful for transfer learning.')
    parser.add_argument('--freeze_model_layers', default=False, action='store_true', help='Whether to freeze the layers from model_layers.')
    parser.add_argument(
        '--continuous_file', default=None, help='Path to a file containing continuous values from which a output TensorMap will be made.'
        'Note that setting this argument has the effect of linking the first output_tensors'
        'argument to the TensorMap made from this file.',
    )

    # Data selection parameters
    parser.add_argument('--continuous_file_column', default=None, help='Column header in file from which a continuous TensorMap will be made.')
    parser.add_argument('--continuous_file_normalize', default=False, action='store_true', help='Whether to normalize a continuous TensorMap made from a file.')
    parser.add_argument(
        '--continuous_file_discretization_bounds', default=[], nargs='*', type=float,
        help='Bin boundaries to use to discretize a continuous TensorMap read from a file.',
    )
    parser.add_argument(
        '--categorical_field_ids', nargs='*', default=[], type=int,
        help='List of field ids from which input features will be collected.',
    )
    parser.add_argument(
        '--continuous_field_ids', nargs='*', default=[], type=int,
        help='List of field ids from which continuous real-valued input features will be collected.',
    )
    parser.add_argument('--include_array', default=False, action='store_true', help='Include array idx for UKBB phenotypes.')
    parser.add_argument('--include_instance', default=False, action='store_true', help='Include instances for UKBB phenotypes.')
    parser.add_argument('--min_values', default=10, type=int, help='Per feature size minimum.')
    parser.add_argument('--min_samples', default=3, type=int, help='Min number of samples to require for calculating correlations.')
    parser.add_argument(
        '--max_samples', type=int, default=None,
        help='Max number of samples to use for tensor reporting -- all samples are used if not specified.',
    )
    parser.add_argument('--mri_field_ids', default=['20208', '20209'], nargs='*', help='Field id for MR images.')
    parser.add_argument('--xml_field_ids', default=['20205', '6025'], nargs='*', help='Field id for XMLs of resting and exercise ECG data.')
    parser.add_argument('--max_patients', default=999999, type=int,  help='Maximum number of patient data to read')
    parser.add_argument('--min_sample_id', default=0, type=int, help='Minimum sample id to write to tensor.')
    parser.add_argument('--max_sample_id', default=7000000, type=int, help='Maximum sample id to write to tensor.')
    parser.add_argument('--max_slices', default=999999, type=int, help='Maximum number of dicom slices to read')
    parser.add_argument('--dicom_series', default='cine_segmented_sax_b6', help='Maximum number of dicom slices to read')
    parser.add_argument(
        '--b_slice_force', default=None,
        help='If set, will only load specific b slice for short axis MRI diastole systole tensor maps (i.e b0, b1, b2, ... b10).',
    )
    parser.add_argument(
        '--include_missing_continuous_channel', default=False, action='store_true',
        help='Include missing channels in continuous tensors',
    )
    parser.add_argument(
        '--imputation_method_for_continuous_fields', default=IMPUTATION_RANDOM, help='can be random or mean',
        choices=[IMPUTATION_RANDOM, IMPUTATION_MEAN],
    )

    # Model Architecture Parameters
    parser.add_argument('--x', default=256, type=int, help='x tensor resolution')
    parser.add_argument('--y', default=256, type=int, help='y tensor resolution')
    parser.add_argument('--zoom_x', default=50, type=int, help='zoom_x tensor resolution')
    parser.add_argument('--zoom_y', default=35, type=int, help='zoom_y tensor resolution')
    parser.add_argument('--zoom_width', default=96, type=int, help='zoom_width tensor resolution')
    parser.add_argument('--zoom_height', default=96, type=int, help='zoom_height tensor resolution')
    parser.add_argument('--z', default=48, type=int, help='z tensor resolution')
    parser.add_argument('--t', default=48, type=int, help='Number of time slices')
    parser.add_argument('--mlp_concat', default=False, action='store_true', help='Concatenate input with every multiplayer perceptron layer.')  # TODO: should be the same style as u_connect
    parser.add_argument('--dense_layers', nargs='*', default=[16, 64], type=int, help='List of number of hidden units in neural nets dense layers.')
    parser.add_argument('--dropout', default=0.0, type=float, help='Dropout rate of dense layers must be in [0.0, 1.0].')
    parser.add_argument('--activation', default='relu',  help='Activation function for hidden units in neural nets dense layers.')
    parser.add_argument('--conv_layers', nargs='*', default=[32], type=int, help='List of number of kernels in convolutional layers.')
    parser.add_argument('--conv_x', default=[3], nargs='*', type=int, help='X dimension of convolutional kernel. Filter sizes are specified per layer given by conv_layers and per block given by dense_blocks. Filter sizes are repeated if there are less than the number of layers/blocks.')
    parser.add_argument('--conv_y', default=[3], nargs='*', type=int, help='Y dimension of convolutional kernel. Filter sizes are specified per layer given by conv_layers and per block given by dense_blocks. Filter sizes are repeated if there are less than the number of layers/blocks.')
    parser.add_argument('--conv_z', default=[2], nargs='*', type=int, help='Z dimension of convolutional kernel. Filter sizes are specified per layer given by conv_layers and per block given by dense_blocks. Filter sizes are repeated if there are less than the number of layers/blocks.')
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
    parser.add_argument('--padding', default='same', help='Valid or same border padding on the convolutional layers.')
    parser.add_argument('--dense_blocks', nargs='*', default=[32, 24, 16], type=int, help='List of number of kernels in convolutional layers.')
    parser.add_argument('--block_size', default=3, type=int, help='Number of convolutional layers within a block.')
    parser.add_argument(
        '--u_connect', nargs=2, action='append',
        help='U-Net connect first TensorMap to second TensorMap. They must be the same shape except for number of channels. Can be provided multiple times.',
    )
    parser.add_argument('--aligned_dimension', default=16, type=int, help='Dimensionality of aligned embedded space for multi-modal alignment models.')
    parser.add_argument(
        '--max_parameters', default=9000000, type=int,
        help='Maximum number of trainable parameters in a model during hyperparameter optimization.',
    )
    parser.add_argument('--bottleneck_type', type=str, default=list(BOTTLENECK_STR_TO_ENUM)[0], choices=list(BOTTLENECK_STR_TO_ENUM))
    parser.add_argument('--hidden_layer', default='embed', help='Name of a hidden layer for inspections.')
    parser.add_argument('--language_layer', default='ecg_rest_text', help='Name of TensorMap for learning language models (eg train_char_model).')
    parser.add_argument('--language_prefix', default='ukb_ecg_rest', help='Path prefix for a TensorMap to learn language models (eg train_char_model)')

    # Training and Hyper-Parameter Optimization Parameters
    parser.add_argument('--epochs', default=12, type=int, help='Number of training epochs.')
    parser.add_argument('--batch_size', default=16, type=int, help='Mini batch size for stochastic gradient descent algorithms.')
    parser.add_argument('--train_csv', help='Path to CSV with Sample IDs to reserve for training.')
    parser.add_argument('--valid_csv', help='Path to CSV with Sample IDs to reserve for validation. Takes precedence over valid_ratio.')
    parser.add_argument('--test_csv', help='Path to CSV with Sample IDs to reserve for testing. Takes precedence over test_ratio.')
    parser.add_argument(
        '--valid_ratio', default=0.2, type=float,
        help='Rate of training tensors to save for validation must be in [0.0, 1.0]. '
             'If any of train/valid/test csv is specified, split by ratio is applied on the remaining tensors after reserving tensors given by csvs. '
             'If not specified, default 0.2 is used. If default ratios are used with train_csv, some tensors may be ignored because ratios do not sum to 1.',
    )
    parser.add_argument(
        '--test_ratio', default=0.1, type=float,
        help='Rate of training tensors to save for testing must be in [0.0, 1.0]. '
             'If any of train/valid/test csv is specified, split by ratio is applied on the remaining tensors after reserving tensors given by csvs. '
             'If not specified, default 0.1 is used. If default ratios are used with train_csv, some tensors may be ignored because ratios do not sum to 1.',
    )
    parser.add_argument('--test_steps', default=32, type=int, help='Number of batches to use for testing.')
    parser.add_argument('--training_steps', default=400, type=int, help='Number of training batches to examine in an epoch.')
    parser.add_argument('--validation_steps', default=40, type=int, help='Number of validation batches to examine in an epoch validation.')
    parser.add_argument('--learning_rate', default=0.0002, type=float, help='Learning rate during training.')
    parser.add_argument('--mixup_alpha', default=0, type=float, help='If positive apply mixup and sample from a Beta with this value as shape parameter alpha.')
    parser.add_argument(
        '--label_weights', nargs='*', type=float,
        help='List of per-label weights for weighted categorical cross entropy. If provided, must map 1:1 to number of labels.',
    )
    parser.add_argument(
        '--patience', default=8, type=int,
        help='Early Stopping parameter: Maximum number of epochs to run without validation loss improvements.',
    )
    parser.add_argument(
        '--max_models', default=16, type=int,
        help='Maximum number of models for the hyper-parameter optimizer to evaluate before returning.',
    )
    parser.add_argument('--balance_csvs', default=[], nargs='*', help='Balances batches with representation from sample IDs in this list of CSVs')
    parser.add_argument('--optimizer', default='radam', type=str, help='Optimizer for model training')
    parser.add_argument('--learning_rate_schedule', default=None, type=str, choices=['triangular', 'triangular2'], help='Adjusts learning rate during training.')
    parser.add_argument('--anneal_rate', default=0., type=float, help='Annealing rate in epochs of loss terms during training')
    parser.add_argument('--anneal_shift', default=0., type=float, help='Annealing offset in epochs of loss terms during training')
    parser.add_argument('--anneal_max', default=2.0, type=float, help='Annealing maximum value')

    # Run specific and debugging arguments
    parser.add_argument('--id', default='no_id', help='Identifier for this run, user-defined string to keep experiments organized.')
    parser.add_argument('--random_seed', default=12878, type=int, help='Random seed to use throughout run.  Always use np.random.')
    parser.add_argument('--write_pngs', default=False, action='store_true', help='Write pngs of slices.')
    parser.add_argument('--debug', default=False, action='store_true', help='Run in debug mode.')
    parser.add_argument('--eager', default=False, action='store_true', help='Run tensorflow functions in eager execution mode (helpful for debugging).')
    parser.add_argument('--inspect_model', default=False, action='store_true', help='Plot model architecture, measure inference and training speeds.')
    parser.add_argument('--inspect_show_labels', default=True, action='store_true', help='Plot model architecture with labels for each layer.')
    parser.add_argument('--alpha', default=0.5, type=float, help='Alpha transparency for t-SNE plots must in [0.0-1.0].')
    parser.add_argument('--plot_mode', default='clinical', choices=['clinical', 'full'], help='ECG view to plot for partners ECGs.')
    parser.add_argument("--embed_visualization", help="Method to visualize embed layer. Options: None, tsne, or umap")
    parser.add_argument("--explore_export_errors", default=False, action="store_true", help="Export error_type columns in tensors_all*.csv generated by explore.")
    parser.add_argument('--plot_hist', default=True, help='Plot histograms of continuous tensors in explore mode.')

    # Training optimization options
    parser.add_argument('--num_workers', default=multiprocessing.cpu_count(), type=int, help="Number of workers to use for every tensor generator.")
    parser.add_argument('--cache_size', default=3.5e9/multiprocessing.cpu_count(), type=float, help="Tensor map cache size per worker.")

    # Cross reference arguments
    parser.add_argument(
        '--tensors_source',
        help='Either a csv or directory of hd5 containing a source dataset.',
    )
    parser.add_argument(
        '--tensors_name', default='Tensors',
        help='Name of dataset at tensors, e.g. ECG. '
             'Adds contextual detail to summary CSV and plots.',
    )
    parser.add_argument(
        '--join_tensors', default=['partners_ecg_patientid_clean'], nargs='+',
        help='TensorMap or column name in csv of value in tensors used in join with reference. '
             'Can be more than 1 join value.',
    )
    parser.add_argument(
        '--time_tensor', default='partners_ecg_datetime',
        help='TensorMap or column name in csv of value in tensors to perform time cross-ref on. '
             'Time cross referencing is optional.',
    )
    parser.add_argument(
        '--reference_tensors',
        help='Either a csv or directory of hd5 containing a reference dataset.',
    )
    parser.add_argument(
        '--reference_name', default='Reference',
        help='Name of dataset at reference, e.g. STS. '
             'Adds contextual detail to summary CSV and plots.',
    )
    parser.add_argument(
        '--reference_join_tensors', nargs='+',
        help='TensorMap or column name in csv of value in reference used in join in tensors. '
             'Can be more than 1 join value.',
    )
    parser.add_argument(
        '--reference_start_time_tensor', action='append', nargs='+',
        help='TensorMap or column name in csv of start of time window in reference. '
             'Define multiple time windows by using this argument more than once. '
             'The number of time windows must match across all time window arguments. '
             'An integer can be provided as a second argument to specify an offset to the start time. '
             'e.g. tStart -30',
    )
    parser.add_argument(
        '--reference_end_time_tensor', action='append', nargs='+',
        help='TensorMap or column name in csv of end of time window in reference. '
             'Define multiple time windows by using this argument more than once. '
             'The number of time windows must match across all time window arguments. '
             'An integer can be provided as a second argument to specify an offset to the end time. '
             'e.g. tEnd 30',
    )
    parser.add_argument(
        '--window_name', action='append',
        help='Name of time window. By default, the name of the window is the index of the window. '
             'Define multiple time windows by using this argument more than once. '
             'The number of time windows must match across all time window arguments.',
    )
    parser.add_argument(
        '--order_in_window', action='append', choices=['newest', 'oldest', 'random'],
        help='If specified, exactly --number_in_window rows with join tensor are used in time window. '
             'Defines which source tensors in a time series to use in time window. '
             'Define multiple time windows by using this argument more than once. '
             'The number of time windows must match across all time window arguments.',
    )
    parser.add_argument(
        '--number_per_window', type=int, default=1,
        help='Minimum number of rows with join tensor to use in each time window. '
             'By default, 1 tensor is used for each window.',
    )
    parser.add_argument(
        '--match_any_window', action='store_true', default=False,
        help='If specified, join tensor does not need to be found in every time window. '
             'Join tensor needs only be found in at least 1 time window. '
             'Default only use rows with join tensor that appears across all time windows.',
    )
    parser.add_argument(
        '--reference_labels', nargs='+',
        help='TensorMap or column name of values in csv to report distribution on, e.g. mortality. '
             'Label distribution reporting is optional. Can list multiple labels to report.',
    )

    args = parser.parse_args()
    _process_args(args)
    return args


def _get_tmap(name: str, needed_tensor_maps: List[str]) -> TensorMap:
    """
    This allows tensor_maps_by_script to only be imported if necessary, because it's slow.
    """
    if name in TMAPS:
        return TMAPS[name]

    TMAPS.update(build_partners_tensor_maps(needed_tensor_maps))
    if name in TMAPS:
        return TMAPS[name]

    TMAPS.update(build_cardiac_surgery_tensor_maps(needed_tensor_maps))
    if name in TMAPS:
        return TMAPS[name]

    TMAPS.update(build_partners_time_series_tensor_maps(needed_tensor_maps))
    if name in TMAPS:
        return TMAPS[name]

    from ml4cvd.tensor_maps_partners_ecg import TMAPS as partners_tmaps
    TMAPS.update(partners_tmaps)
    if name in TMAPS:
        return TMAPS[name]

    from ml4cvd.tensor_maps_partners_ecg_labels import TMAPS as partners_label_tmaps
    TMAPS.update(partners_label_tmaps)
    if name in TMAPS:
        return TMAPS[name]

    from ml4cvd.tensor_maps_by_script import TMAPS as script_tmaps
    TMAPS.update(script_tmaps)

    from ml4cvd.tensor_maps_by_script import TMAPS as script_tmaps
    TMAPS.update(script_tmaps)

    return TMAPS[name]


def _process_u_connect_args(u_connect: Optional[List[List]]) -> Dict[TensorMap, Set[TensorMap]]:
    u_connect = u_connect or []
    new_u_connect = defaultdict(set)
    for connect_pair in u_connect:
        tmap_key_in, tmap_key_out = connect_pair[0], connect_pair[1]
        tmap_in, tmap_out = _get_tmap(tmap_key_in, []), _get_tmap(tmap_key_out, [])
        if tmap_in.shape[:-1] != tmap_out.shape[:-1]:
            raise TypeError(f'u_connect of {tmap_in} {tmap_out} requires matching shapes besides channel dimension.')
        if tmap_in.axes() < 2 or tmap_out.axes() < 2:
            raise TypeError(f'Cannot u_connect 1d TensorMaps ({tmap_in} {tmap_out}).')
        new_u_connect[tmap_in].add(tmap_out)
    return new_u_connect


def _process_args(args):
    now_string = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    args_file = os.path.join(args.output_folder, args.id, 'arguments_' + now_string + '.txt')
    command_line = f"\n./scripts/tf.sh {' '.join(sys.argv)}\n"
    if not os.path.exists(os.path.dirname(args_file)):
        os.makedirs(os.path.dirname(args_file))
    with open(args_file, 'w') as f:
        f.write(command_line)
        for k, v in sorted(args.__dict__.items(), key=operator.itemgetter(0)):
            f.write(k + ' = ' + str(v) + '\n')
    load_config(args.logging_level, os.path.join(args.output_folder, args.id), 'log_' + now_string, args.min_sample_id)
    args.u_connect = _process_u_connect_args(args.u_connect)
    needed_tensor_maps = args.input_tensors + args.output_tensors + [args.sample_weight] if args.sample_weight else args.input_tensors + args.output_tensors
    args.tensor_maps_in = [_get_tmap(it, needed_tensor_maps) for it in args.input_tensors]
    args.sample_weight = _get_tmap(args.sample_weight, needed_tensor_maps) if args.sample_weight else None
    if args.sample_weight:
        assert args.sample_weight.shape == (1,)

    args.tensor_maps_out = []
    if args.continuous_file is not None:
        # Continuous TensorMap generated from file is given the name specified by the first output_tensors argument
        args.tensor_maps_out.append(
            generate_continuous_tensor_map_from_file(
                args.continuous_file,
                args.continuous_file_column,
                args.output_tensors.pop(0),
                args.continuous_file_normalize,
                args.continuous_file_discretization_bounds,
            ),
        )
    args.tensor_maps_out.extend([_get_tmap(ot, needed_tensor_maps) for ot in args.output_tensors])
    args.tensor_maps_out = parent_sort(args.tensor_maps_out)

    args.bottleneck_type = BOTTLENECK_STR_TO_ENUM[args.bottleneck_type]
    if args.bottleneck_type == BottleneckType.NoBottleNeck:
        check_no_bottleneck(args.u_connect, args.tensor_maps_out)

    if args.learning_rate_schedule is not None and args.patience < args.epochs:
        raise ValueError(f'learning_rate_schedule is not compatible with ReduceLROnPlateau. Set patience > epochs.')

    np.random.seed(args.random_seed)

    logging.info(f"Command Line was: {command_line}")
    logging.info(f"Total TensorMaps: {len(TMAPS)} Arguments are {args}")

    if args.eager:
        import tensorflow as tf
        tf.config.experimental_run_functions_eagerly(True)
