# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Tensor Maps ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from . import defines
import keras.backend as K

def get_tensor_channel_map_from_args(args):
    '''Return tensor mapping dict given args.tensor_name'''
    if not args.tensor_name:
        return None

    if 'read_tensor' == args.tensor_name:
        return get_read_tensor_channel_map()
    elif 'reference' == args.tensor_name:
        return get_tensor_channel_map_1d_dna()
    else:
        raise ValueError('Unknown tensor mapping mode:', args.tensor_name)


def get_tensor_channel_map_1d_dna():
    '''1D Reference tensor with 4 channel DNA encoding.'''
    tensor_map = {}
    for k in defines.DNA_SYMBOLS.keys():
        tensor_map[k] = defines.DNA_SYMBOLS[k]

    return tensor_map


def get_tensor_channel_map_reference_reads():
    '''Read and reference tensor with 4 channel DNA encoding.
    Plus insertions and deletions.
    '''
    tensor_map = {}
    for k in defines.INPUTS_INDEL.keys():
        tensor_map['read_'+k] = defines.INPUTS_INDEL[k]
    for k in defines.INPUTS_INDEL.keys():
        tensor_map['reference_'+k] = len(defines.INPUTS_INDEL) + defines.INPUTS_INDEL[k]

    return tensor_map


def get_read_tensor_channel_map():
    '''Read and reference tensor with 4 channel DNA encoding.
    Also includes read flags.
    '''
    tensor_map = {}
    for k in defines.INPUTS_INDEL.keys():
        tensor_map['read_'+k] = defines.INPUTS_INDEL[k]
    for k in defines.INPUTS_INDEL.keys():
        tensor_map['reference_'+k] = len(defines.INPUTS_INDEL) + defines.INPUTS_INDEL[k]
    tensor_map['flag_bit_4'] = 10
    tensor_map['flag_bit_5'] = 11
    tensor_map['flag_bit_6'] = 12
    tensor_map['flag_bit_7'] = 13
    tensor_map['mapping_quality'] = 14
    return tensor_map


def tensor_shape_from_args(args):
    in_channels = len(get_tensor_channel_map_from_args(args))
    if K.image_data_format() == 'channels_last':
        tensor_shape = (args.read_limit, args.window_size, in_channels)
    else:
        tensor_shape = (in_channels, args.read_limit, args.window_size)
    return tensor_shape

def total_input_channels_from_args(args):
    '''Get the number of channels in the tensor map'''
    return len(get_tensor_channel_map_from_args(args))

