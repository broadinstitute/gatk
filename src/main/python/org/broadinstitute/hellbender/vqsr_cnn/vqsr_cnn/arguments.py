import argparse
import numpy as np

import keras.backend as K

from . import defines

def parse_args():
    """Parse command line arguments.

    The args namespace is used promiscuously in this module.
    Its fields control the tensor definition, dataset generation, training, file I/O and evaluation.
    Some of the fields are typically dicts or lists that are not actually set on the command line,
    but via a companion argument also in the namespace.
    For example, input_symbols is set via the input_symbol_set string
    and, annotations is set via the annotation_set string.
    Here we also seed the random number generator.
    The keras image data format is set here as well via the channels_last or channels_first arguments.

    Returns:
        namespace: The args namespace that is used throughout this module.
    """
    parser = argparse.ArgumentParser()

    # Tensor defining arguments
    parser.add_argument('--tensor_name', default='read_tensor', choices=defines.TENSOR_MAPS_1D+defines.TENSOR_MAPS_2D,
                        help='String key which identifies the map from tensor channels to their meaning.')
    parser.add_argument('--labels', default=defines.SNP_INDEL_LABELS,
                        help='Dict mapping label names to their index within label tensors.')
    parser.add_argument('--input_symbol_set', default='dna_indel', choices=defines.INPUT_SYMBOLS.keys(),
                        help='Key which maps to an input symbol to index mapping.')
    parser.add_argument('--input_symbols', help='Dict mapping input symbols to their index within input tensors, '
                        + 'initialised via input_symbols_set argument')

    parser.add_argument('--batch_size', default=32, type=int,
                        help='Mini batch size for stochastic gradient descent algorithms.')
    parser.add_argument('--read_limit', default=128, type=int,
                        help='Maximum number of reads to load.')
    parser.add_argument('--window_size', default=128, type=int,
                        help='Size of sequence window to use as input, typically centered at a variant.')
    parser.add_argument('--base_quality_mode', default='phot', choices=['phot', 'phred', '1hot'],
                        help='How to treat base qualities, must be in [phot, phred, 1hot]')
    parser.add_argument('--channels_last', default=True, dest='channels_last', action='store_true',
                        help='Store the channels in the last axis of tensors, tensorflow->true, theano->false')
    parser.add_argument('--channels_first', dest='channels_last', action='store_false',
                        help='Store the channels in the first axis of tensors, tensorflow->false, theano->true')

    # Annotation arguments
    parser.add_argument('--annotations', help='Array of annotation names, initialised via annotation_set argument')
    parser.add_argument('--annotation_set', default='best_practices', choices=defines.ANNOTATIONS.keys(),
                        help='Key which maps to an annotations list (or _ to ignore annotations).')

    # Dataset generation related arguments
    parser.add_argument('--samples', default=500, type=int,
                        help='Maximum number of data samples to write or load.')
    parser.add_argument('--downsample_snps', default=1.0, type=float,
                        help='Rate of SNP examples that are kept must be in [0.0, 1.0].')
    parser.add_argument('--downsample_indels', default=1.0, type=float,
                        help='Rate of INDEL examples that are kept must be in [0.0, 1.0].')
    parser.add_argument('--downsample_not_snps', default=1.0, type=float,
                        help='Rate of NOT_SNP examples that are kept must be in [0.0, 1.0].')
    parser.add_argument('--downsample_not_indels', default=1.0, type=float,
                        help='Rate of NOT_INDEL examples that are kept must be in [0.0, 1.0].')
    parser.add_argument('--downsample_reference', default=0.001, type=float,
                        help='Rate of reference genotype examples that are kept must be in [0.0, 1.0].')
    parser.add_argument('--downsample_homozygous', default=0.001, type=float,
                        help='Rate of homozygous genotypes that are kept must be in [0.0, 1.0].')
    parser.add_argument('--start_pos', default=0, type=int,
                        help='Genomic position start for parallel tensor writing.')
    parser.add_argument('--end_pos', default=0, type=int,
                        help='Genomic position end for parallel tensor writing.')
    parser.add_argument('--skip_positive_class', default=False, action='store_true',
                        help='Whether to skip positive examples when writing tensors.')
    parser.add_argument('--chrom', help='Chromosome to load for parallel tensor writing.')

    # I/O files and directories: vcfs, bams, beds, hd5, fasta
    parser.add_argument('--output_dir', default='./weights/', help='Directory to write models or other data out.')
    parser.add_argument('--image_dir', default=None, help='Directory to write images and plots to.')
    parser.add_argument('--reference_fasta', help='The reference FASTA file (e.g. HG19 or HG38).')
    parser.add_argument('--weights_hd5', default='',
                        help='A hd5 file of weights to initialize a model, will use all layers with names that match.')
    parser.add_argument('--architecture', default='',
                        help='A json file specifying semantics and architecture of a neural net.')
    parser.add_argument('--bam_file',
                        help='Path to a BAM file to train from or generate tensors with.')
    parser.add_argument('--train_vcf',
                        help='Path to a VCF that has verified true calls from NIST, platinum genomes, etc.')
    parser.add_argument('--input_vcf',
                        help='Haplotype Caller or VQSR generated VCF with raw annotation values [and quality scores].')
    parser.add_argument('--output_vcf', default=None,
                        help='Optional VCF to write to.')
    parser.add_argument('--bed_file',
                        help='Bed file specifying high confidence intervals associated with args.train_vcf.')
    parser.add_argument('--data_dir',
                        help='Directory of tensors, must be split into test/valid/train directories'
                            +'with subdirectories for each label.')

    # Training and optimization related arguments
    parser.add_argument('--epochs', default=25, type=int,
                        help='Number of epochs, typically passes through the entire dataset, not always well-defined.')
    parser.add_argument('--batch_normalization', default=False, action='store_true',
                        help='Mini batch normalization layers after convolutions.')
    parser.add_argument('--patience', default=4, type=int,
                        help='Maximum number of epochs to run without validation loss improvements (Early Stopping).')
    parser.add_argument('--training_steps', default=80, type=int,
                        help='Number of training batches to examine in an epoch.')
    parser.add_argument('--validation_steps', default=40, type=int,
                        help='Number of validation batches to examine in an epoch validation.')
    parser.add_argument('--iterations', default=5, type=int,
                        help='Generic iteration limit for hyperparameter optimization, animation, and other counts.')
    parser.add_argument('--tensor_board', default=False, action='store_true',
                        help='Add the tensor board callback.')

    # Evaluation related arguments
    parser.add_argument('--score_keys', nargs='+', default=['VQSLOD'],
                        help='List of variant score keys for performance comparisons.')
    parser.add_argument('--tranches', nargs='+', default=[100, 99.9, 99, 95, 90], type=float,
                        help='List of variant score keys for performance comparisons.')

    # Run specific arguments
    parser.add_argument('--mode', help='High level recipe: write tensors, train, test or evaluate models.')
    parser.add_argument('--id', default='no_id',
                        help='Identifier for this run, user-defined string to keep experiments organized.')
    parser.add_argument('--random_seed', default=12878, type=int,
                        help='Random seed to use throughout run.  Always use np.random.')

    # Parse, print, set annotations and seed
    args = parser.parse_args()
    args.annotations = annotations_from_args(args)
    args.input_symbols = input_symbols_from_args(args)
    np.random.seed(args.random_seed)

    if args.channels_last:
        K.set_image_data_format('channels_last')
    else:
        K.set_image_data_format('channels_first')

    print('Arguments are', args)
    return args


def annotations_from_args(args):
    """Get list of annotations corresponding to the args.annotation_set.

    The annotation_set argument allows us to name commonly used groups of annotations
    without having to specify each annotation individually.

    Arguments:
        args.annotation_set: The key for which annotation set to use.

    Returns:
        list: Annotation strings as they appear in a VCF info/format field or None.
    """
    if args.annotation_set and args.annotation_set in defines.ANNOTATIONS:
        return defines.ANNOTATIONS[args.annotation_set]
    return None


def input_symbols_from_args(args):
    """Get dictionary mapping input data symbols to indices in the input tensor.

    Arguments:
        args.input_symbol_set: The key for the symbol set to use.

    Returns:
        dict: if there is a input symbol dict otherwise None
    """
    if args.input_symbol_set and args.input_symbol_set in defines.INPUT_SYMBOLS:
        return defines.INPUT_SYMBOLS[args.input_symbol_set]
    return None



def weight_path_from_args(args):
    """Create a weight file name from the command line arguments.

    Arguments:
        args.output_dir: The directory where the file will be saved
        args.id: The name of the file is this run's id with tensor suffix as file extension
    """
    save_weight_hd5 = args.output_dir + args.id + defines.TENSOR_SUFFIX
    return save_weight_hd5

