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
import logging
import argparse
import operator
import datetime
import numpy as np

from ml4cvd.logger import load_config
from ml4cvd.tensor_map_maker import generate_multi_field_continuous_tensor_map
from ml4cvd.tensor_maps_by_script import TMAPS


CATEGORICAL_PHENOTYPES = [54, 924, 943, 971, 981, 1011, 1100, 1239, 1249, 1259, 1329, 1339, 1349, 1359, 1369, 1379, 1389, 1408, 1418, 1428, 1448, 1468, 1478, 1508, 1518, 1528, 1538, 1548, 1558, 1618, 1628, 1647, 1677, 1687, 1697, 1707, 1717, 1727, 1747, 1757, 1767, 1777, 1787, 1797, 1835, 2178, 2188, 2207, 2247, 2316, 2306, 2415, 2443, 2453, 2463, 2473, 2674, 2694, 2724, 2784, 2814, 2877, 3079, 3616, 3637, 3773, 3799, 4717, 4825, 4935, 4957, 4968, 4979, 4990, 5001, 5012, 6015, 6017, 6148, 6149, 6150, 6152, 6153, 6154, 6155, 6157, 6159, 6162, 6164, 6177, 6179, 20001, 20003, 20004, 20116, 22001, 22609, 22610, 22611, 22612, 22613, 22614, 22615, 22616, 22650]
CONTINUOUS_PHENOTYPES = [34, 48, 49, 50, 102, 137, 864, 874, 884, 894, 904, 914, 1070, 1080, 1090, 1309, 1438, 1478, 1488, 1498, 1568, 4237, 4239, 4288, 4407, 20015, 20023, 21001, 21002, 22200, 22602, 22603, 22003, 22009, 23098, 23099, 23104, 23106, 23120, 23128, 30000, 30010, 30020, 30030, 30040, 30050, 30060, 30070 ]


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mode', default='mlp', help='What would you like to do?')

    # Config arguments
    parser.add_argument("--logging_level", default='INFO',
                        help="Logging level. Overrides any configuration given in the logging configuration file.",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])

    # Tensor Map arguments
    parser.add_argument('--input_tensors', default=[], nargs='+')
    parser.add_argument('--output_tensors', default=[], nargs='+')
    parser.add_argument('--input_continuous_tensors', default=[], nargs='+', help='Continuous tensor maps to be combined.')
    parser.add_argument('--tensor_maps_in', default=[], help='Do not set this directly. Use input_tensors')
    parser.add_argument('--tensor_maps_out', default=[], help='Do not set this directly. Use output_tensors')

    # Input and Output files and directories
    parser.add_argument('--db', default='/mnt/disks/data/raw/sql/ukbb7089.r10data.db',
                        help='Path to sqlite3 data base file.')
    parser.add_argument('--bigquery_credentials_file',
                        default='/mnt/ml4cvd/projects/jamesp/bigquery/bigquery-viewer-credentials.json',
                        help='Path to service account credentials for looking up BigQuery tables.')
    parser.add_argument('--bigquery_dataset',
                        default='broad-ml4cvd.ukbb7089_r10data',
                        help='BigQuery dataset containing tables we want to query.')
    parser.add_argument('--xml_folder', default='/mnt/disks/data/raw/ecgs/',
        help='Path to folder of XMLs of ECG data.')
    parser.add_argument('--zip_folder', default='/mnt/disks/data/raw/mris/cardiac/',
        help='Path to folder of zipped dicom images.')
    parser.add_argument('--phenos_folder', default='/mnt/disks/data/raw/phenotypes/',
        help='Path to folder of phenotype defining CSVs.')    
    parser.add_argument('--dicoms', default='./dicoms/',
        help='Path to folder of dicoms ( dicoms/labels/sample_id/field_id/*dcm.')
    parser.add_argument('--icd_csv', default='/mnt/disks/data/raw/tsvs/modified.zmerge.prs.full.csv.20190219',
        help='Path to CSV with ICD status for UKBB Sample IDs')
    parser.add_argument('--volume_csv', default='/mnt/disks/data/raw/tsvs/round2_4sdfixedpoint_parsed_lvedv_lvesv_lvef.tsv',
        help='Path to left ventricle volumes')
    parser.add_argument('--app_csv', default='/mnt/disks/data/raw/tsvs/ukb_app17488_app7089_link.csv',
        help='Path to file used to link sample IDs between UKBB applications 17488 and 7089')
    parser.add_argument('--tensors', default='/mnt/disks/data/generated/tensors/test/2019-03-21/',
        help='Path to folder containing tensors, or where tensors will be written.')
    parser.add_argument('--output_folder', default='./recipes_output/',
        help='Path to output folder for recipes.py runs.')
    parser.add_argument('--model_file',
        help='Path to a saved model architecture and weights (hd5).')
    parser.add_argument('--model_files', nargs='*', default=[],
        help='List of paths to saved model architectures and weights (hd5).')
    parser.add_argument('--model_layers',
        help='Path to a model file (hd5) which will be loaded by layer, useful for transfer learning.')
    parser.add_argument('--model_freeze',
        help='Like the model_layers argument, except all loaded layers will also be frozen.')
    parser.add_argument('--lv_mass_csv', default='/mnt/ml4cvd/projects/jamesp/data/returned_lv_mass.tsv',
        help='Path to left ventricular mass and other cardiac MRI readouts on ~5000 people returned from app 2964')

    # Data selection parameters
    parser.add_argument('--categorical_field_ids', nargs='*', default=CATEGORICAL_PHENOTYPES, type=int,
        help='List of field ids from which input features will be collected.')
    parser.add_argument('--continuous_field_ids', nargs='*', default=CONTINUOUS_PHENOTYPES, type=int,
        help='List of field ids from which continuous real-valued input features will be collected.')
    parser.add_argument('--include_array', default=False, action='store_true',
        help='Include array idx for UKBB phenotypes.')
    parser.add_argument('--include_instance', default=False, action='store_true',
        help='Include instances for UKBB phenotypes.')
    parser.add_argument('--min_values', default=10, type=int,
        help='Per feature size minimum.')
    parser.add_argument('--mri_field_ids', default=['20208', '20209'], nargs='*',
        help='Field id for MR images.')
    parser.add_argument('--xml_field_ids', default=['20205', '6025'], nargs='*',
        help='Field id for XMLs of resting and exercise ECG data.')
    parser.add_argument('--max_patients', default=999999, type=int,
        help='Maximum number of patient data to read')
    parser.add_argument('--min_sample_id', default=0, type=int,
        help='Minimum sample id to write to tensor.')
    parser.add_argument('--max_sample_id', default=7000000, type=int,
        help='Maximum sample id to write to tensor.')
    parser.add_argument('--max_slices', default=999999, type=int,
        help='Maximum number of dicom slices to read')
    parser.add_argument('--b_slice_force', default=None,
        help='If set, will only load specific b slice for short axis MRI diastole systole tensor maps (i.e b0, b1, b2, ... b10).')
    parser.add_argument('--include_heart_zoom', default=False, action='store_true',
        help='Include the heart zoom')

    # Model Architecture Parameters
    parser.add_argument('--x', default=256, type=int,
        help='x tensor resolution')
    parser.add_argument('--y', default=256, type=int,
        help='y tensor resolution')
    parser.add_argument('--zoom_x', default=50, type=int,
        help='zoom_x tensor resolution')
    parser.add_argument('--zoom_y', default=35, type=int,
        help='zoom_y tensor resolution')
    parser.add_argument('--zoom_width', default=96, type=int,
        help='zoom_width tensor resolution')
    parser.add_argument('--zoom_height', default=96, type=int,
        help='zoom_height tensor resolution')
    parser.add_argument('--z', default=48, type=int,
        help='z tensor resolution')
    parser.add_argument('--t', default=48, type=int,
        help='Number of time slices')
    parser.add_argument('--mlp_concat', default=False, action='store_true',
        help='Concatenate input with every multiplayer perceptron layer.')
    parser.add_argument('--channels_last', default=True, dest='channels_last', action='store_true',
        help='Store the channels in the last axis of tensors, tensorflow->true, theano->false')
    parser.add_argument('--channels_first', dest='channels_last', action='store_false',
        help='Store the channels in the first axis of tensors, tensorflow->false, theano->true')
    parser.add_argument('--dense_layers', nargs='*', default=[16, 64], type=int,
        help='List of number of hidden units in neural nets dense layers.')
    parser.add_argument('--dropout', default=0.0, type=float,
        help='Dropout rate of dense layers must be in [0.0, 1.0].')
    parser.add_argument('--activation', default='relu',
        help='Activation function for hidden units in neural nets dense layers.')
    parser.add_argument('--conv_layers', nargs='*', default=[32], type=int,
        help='List of number of kernels in convolutional layers.')
    parser.add_argument('--conv_x', default=3, type=int,
        help='X dimension of convolutional kernel.')
    parser.add_argument('--conv_y', default=3, type=int,
        help='Y dimension of convolutional kernel.')
    parser.add_argument('--conv_z', default=2, type=int,
        help='Z dimension of convolutional kernel.')
    parser.add_argument('--conv_width', default=71, type=int,
        help='Width of convolutional kernel for 1D CNNs.')
    parser.add_argument('--conv_bn', default=False, action='store_true',
        help='Batch normalize convolutional layers.')
    parser.add_argument('--conv_dropout', default=0.0, type=float,
        help='Dropout rate of convolutional kernels must be in [0.0, 1.0].')
    parser.add_argument('--max_pools', nargs='*', default=[], type=int,
        help='List of maxpooling layers.')
    parser.add_argument('--pool_x', default=2, type=int,
        help='Pooling size in the x-axis, if 1 no pooling will be performed.')
    parser.add_argument('--pool_y', default=2, type=int,
        help='Pooling size in the y-axis, if 1 no pooling will be performed.')
    parser.add_argument('--pool_z', default=1, type=int,
        help='Pooling size in the z-axis, if 1 no pooling will be performed.')
    parser.add_argument('--res_layers', nargs='*', default=[], type=int,
        help='List of residual layers.')
    parser.add_argument('--padding', default='same',
        help='Valid or same border padding on the convolutional layers.')
    parser.add_argument('--dense_blocks', nargs='*', default=[32, 24, 16], type=int,
        help='List of number of kernels in convolutional layers.')
    parser.add_argument('--block_size', default=3, type=int,
        help='Number of convolutional layers within a block.')
    parser.add_argument('--u_connect', default=False, action='store_true',
        help='Connect early convolutional layers to later ones of the same size, as in U-Net.')
    parser.add_argument('--aligned_dimension', default=16, type=int,
        help='Dimensionality of aligned embedded space for multi-modal alignment models.')
    parser.add_argument('--max_parameters', default=9000000, type=int,
        help='Maximum number of trainable parameters in a model during hyperparameter optimization.')

    # Training and Hyper-Parameter Optimization Parameters
    parser.add_argument('--epochs', default=12, type=int,
        help='Number of training epochs.')
    parser.add_argument('--batch_size', default=16, type=int,
        help='Mini batch size for stochastic gradient descent algorithms.')
    parser.add_argument('--valid_ratio', default=0.2, type=float,
        help='Rate of training tensors to save for validation must be in [0.0, 1.0].')
    parser.add_argument('--test_ratio', default=0.1, type=float,
        help='Rate of training tensors to save for testing [0.0, 1.0].')
    parser.add_argument('--test_steps', default=32, type=int,
        help='Number of batches to use for testing.')
    parser.add_argument('--training_steps', default=400, type=int,
        help='Number of training batches to examine in an epoch.')
    parser.add_argument('--validation_steps', default=40, type=int,
        help='Number of validation batches to examine in an epoch validation.')
    parser.add_argument('--learning_rate', default=0.001, type=float,
        help='Learning rate during training.')
    parser.add_argument('--label_weights', nargs='*', type=float,
        help='List of per-label weights for weighted categorical cross entropy. If provided, must map 1:1 to number of labels.')
    parser.add_argument('--patience', default=8, type=int,
        help='Early Stopping parameter: Maximum number of epochs to run without validation loss improvements.')
    parser.add_argument('--max_models', default=16, type=int,
        help='Maximum number of models for the hyper-parameter optimizer to evaluate before returning.')
    parser.add_argument('--balance_by_icds', default=[], nargs='*', type=int, #[135, 125], [92, 116, 154, 174]
        help='Balances minibatches with samples with different ICD status')

    # Run specific and debugging arguments
    parser.add_argument('--id', default='no_id',
        help='Identifier for this run, user-defined string to keep experiments organized.')
    parser.add_argument('--random_seed', default=12878, type=int,
        help='Random seed to use throughout run.  Always use np.random.')
    parser.add_argument('--write_pngs', default=False, action='store_true',
        help='Write pngs of slices.')
    parser.add_argument('--debug', default=False, action='store_true',
        help='Run in debug mode.')
    parser.add_argument('--inspect_model', default=False, action='store_true',
        help='Plot model architecture, measure inference and training speeds.')
    parser.add_argument('--inspect_show_labels', default=True, action='store_true',
        help='Plot model architecture with labels for each layer.')
    parser.add_argument('--num_samples', type=int, default=None,
        help='Max number of samples to use for tensor reporting -- all samples are used if not specified.')

    args = parser.parse_args()

    _process_args(args)

    return args


def _process_args(args):
    if len(args.input_continuous_tensors) > 0:
        multi_field_tensor_map = generate_multi_field_continuous_tensor_map(args.input_continuous_tensors)
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors] + [multi_field_tensor_map]
    else:
        args.tensor_maps_in = [TMAPS[it] for it in args.input_tensors]
    args.tensor_maps_out = [TMAPS[ot] for ot in args.output_tensors]
    np.random.seed(args.random_seed)

    now_string = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M')
    args_file = os.path.join(args.output_folder, args.id, 'arguments_'+now_string+'.txt')
    if not os.path.exists(os.path.dirname(args_file)):
        os.makedirs(os.path.dirname(args_file))
    with open(args_file, 'w') as f:
        for k, v in sorted(args.__dict__.items(), key=operator.itemgetter(0)):
            f.write(k + ' = ' + str(v) + '\n')

    load_config(args.logging_level, os.path.join(args.output_folder, args.id), 'log_'+now_string, args.min_sample_id)
    logging.info('Total TensorMaps:{} Arguments are {}'.format(len(TMAPS), args))

