# Imports
import os
import math
import h5py
import numpy as np
from collections import namedtuple
from typing import List, Tuple, Dict, TextIO

from gatktool import tool

# Keras Imports
import keras
import keras.backend as K

# Package Imports
from . import defines
from . import tensor_maps

Variant = namedtuple("Variant", "contig pos ref alt type")
Read = namedtuple("Read", "seq qual cigar reverse mate_reverse first mapping_quality reference_start")

READ_BASES_FIFO_INDEX = 0
READ_QUAL_FIFO_INDEX = 1
READ_CIGAR_FIFO_INDEX = 2
READ_REVERSE_FIFO_INDEX = 3
READ_MATE_REVERSE_FIFO_INDEX = 4
READ_FIRST_IN_PAIR_FIFO_INDEX = 5
READ_MQ_FIFO_INDEX = 6
READ_REF_START_FIFO_INDEX = 7
READ_ELEMENTS = 8  # The number of fields of the namedtuple defined above

CONTIG_FIFO_INDEX = 0
POS_FIFO_INDEX = 1
REF_FIFO_INDEX = 2
ALT_FIFO_INDEX = 3
REF_STRING_FIFO_INDEX = 4
ANNOTATION_FIFO_INDEX = 5
VARIANT_TYPE_FIFO_INDEX = 6
VARIANT_FIFO_FIELDS = 7



CIGAR_CODES_TO_COUNT = [
    defines.CIGAR_CODE['M'], defines.CIGAR_CODE['I'], defines.CIGAR_CODE['S'], defines.CIGAR_CODE['D']
]

p_lut = np.zeros((256,))
not_p_lut = np.zeros((256,))

for i in range(256):
    exponent = float(-i) / 10.0
    p_lut[i] = 1.0 - (10.0**exponent)
    not_p_lut[i] = (1.0 - p_lut[i]) / 3.0


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Inference ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def score_and_write_batch(model: keras.Model,
                          file_out: TextIO,
                          batch_size: int,
                          python_batch_size: int,
                          tensor_type: str,
                          annotation_set: str,
                          window_size: int,
                          read_limit: int,
                          tensor_dir: str = '') -> None:
    """Score a batch of variants with a CNN model. Write tab delimited temp file with scores.

    This function is tightly coupled with the CNNScoreVariants.java
    It requires data written to the fifo in the order given by transferToPythonViaFifo

    Arguments
        model: a keras model
        file_out: The temporary VCF-like file where variants scores will be written
        batch_size: The total number of variants available in the fifo
        python_batch_size: the number of variants to process in each inference
        tensor_type: The name for the type of tensor to make
        annotation_set: The name for the set of annotations to use
        window_size: The size of the context window of genomic bases, i.e the width of the tensor
        read_limit: The maximum number of reads to encode in a tensor, i.e. the height of the tensor
        tensor_dir : If this path exists write hd5 files for each tensor (optional for debugging)
    """
    annotation_batch = []
    reference_batch = []
    variant_types = []
    variant_data = []
    read_batch = []
    for _ in range(batch_size):
        fifo_line = tool.readDataFIFO()
        fifo_data = fifo_line.split(defines.DATA_TYPE_SEPARATOR)

        variant_data.append(fifo_data[CONTIG_FIFO_INDEX] + defines.DATA_TYPE_SEPARATOR
                            + fifo_data[POS_FIFO_INDEX] + defines.DATA_TYPE_SEPARATOR
                            + fifo_data[REF_FIFO_INDEX] + defines.DATA_TYPE_SEPARATOR + fifo_data[ALT_FIFO_INDEX])
        reference_batch.append(reference_string_to_tensor(fifo_data[REF_STRING_FIFO_INDEX]))
        annotation_batch.append(annotation_string_to_tensor(annotation_set, fifo_data[ANNOTATION_FIFO_INDEX]))
        variant_types.append(fifo_data[VARIANT_TYPE_FIFO_INDEX].strip())

        fifo_idx = VARIANT_FIFO_FIELDS
        if tensor_type in defines.TENSOR_MAPS_2D and len(fifo_data) > fifo_idx:
            read_tuples = []
            var = Variant(fifo_data[CONTIG_FIFO_INDEX], int(fifo_data[POS_FIFO_INDEX]), fifo_data[POS_FIFO_INDEX],
                          fifo_data[ALT_FIFO_INDEX], fifo_data[VARIANT_TYPE_FIFO_INDEX])
            while fifo_idx+READ_ELEMENTS <= len(fifo_data):
                read_tuples.append(
                    Read(fifo_data[fifo_idx + READ_BASES_FIFO_INDEX],
                         list(map(int, fifo_data[fifo_idx+READ_QUAL_FIFO_INDEX].split(defines.DATA_VALUE_SEPARATOR))),
                         fifo_data[fifo_idx+READ_CIGAR_FIFO_INDEX],
                         bool_from_java(fifo_data[fifo_idx+READ_REVERSE_FIFO_INDEX]),
                         bool_from_java(fifo_data[fifo_idx+READ_MATE_REVERSE_FIFO_INDEX]),
                         bool_from_java(fifo_data[fifo_idx+READ_FIRST_IN_PAIR_FIFO_INDEX]),
                         int(fifo_data[fifo_idx+READ_MQ_FIFO_INDEX]),
                         int(fifo_data[fifo_idx+READ_REF_START_FIFO_INDEX])))
                fifo_idx += READ_ELEMENTS
            _, ref_start, _ = get_variant_window(window_size, var)
            insert_dict = get_inserts(read_tuples, var, window_size)
            tensor = read_tuples_to_tensor(read_tuples, ref_start, insert_dict, tensor_type, window_size, read_limit)
            reference_sequence_into_tensor(fifo_data[4], tensor, insert_dict, window_size, read_limit)
            if os.path.exists(tensor_dir):
                _write_tensor_to_hd5(tensor, annotation_batch[-1], fifo_data[0], fifo_data[1], fifo_data[6],
                                     tensor_type, annotation_set, tensor_dir)
            read_batch.append(tensor)

    if tensor_type in defines.TENSOR_MAPS_1D:
        predictions = model.predict([np.array(reference_batch), np.array(annotation_batch)],
                                    batch_size=python_batch_size)
    elif tensor_type in defines.TENSOR_MAPS_2D:
        predictions = model.predict(
            [np.array(read_batch), np.array(annotation_batch)], batch_size=python_batch_size)
    else:
        raise ValueError('Unknown tensor mapping.  Check architecture file.', tensor_type)

    indel_scores = predictions_to_indel_scores(predictions)
    snp_scores = predictions_to_snp_scores(predictions)

    for i in range(batch_size):
        if 'SNP' == variant_types[i]:
            file_out.write(variant_data[i] + defines.DATA_TYPE_SEPARATOR + '{0:.3f}'.format(snp_scores[i]) + '\n')
        elif 'INDEL' == variant_types[i]:
            file_out.write(variant_data[i] + defines.DATA_TYPE_SEPARATOR + '{0:.3f}'.format(indel_scores[i]) + '\n')
        else:
            file_out.write(variant_data[i] + defines.DATA_TYPE_SEPARATOR
                           + '{0:.3f}'.format(max(snp_scores[i], indel_scores[i])) + '\n')


def reference_string_to_tensor(reference: str) -> np.ndarray:
    dna_data = np.zeros((len(reference), len(defines.DNA_SYMBOLS)))
    for i,b in enumerate(reference):
        if b in defines.DNA_SYMBOLS:
            dna_data[i, defines.DNA_SYMBOLS[b]] = 1.0
        elif b in defines.AMBIGUITY_CODES:
            dna_data[i] = defines.AMBIGUITY_CODES[b]
        elif b == '\x00':
            break
        else:
            raise ValueError('Error! Unknown code:', b)
    return dna_data


def annotation_string_to_tensor(annotation_set: str, annotation_string: str) -> np.ndarray:
    name_val_pairs = annotation_string.split(defines.ANNOTATION_SEPARATOR)
    annotation_names = annotation_set.split(defines.DATA_VALUE_SEPARATOR)
    name_val_arrays = [p.split(defines.ANNOTATION_SET_STRING) for p in name_val_pairs]
    annotation_map = {str(p[0]).strip(): p[1] for p in name_val_arrays if len(p) > 1}
    annotation_data = np.zeros((len(annotation_names),))
    for ii, a in enumerate(annotation_names):
        if a in annotation_map and not math.isnan(float(annotation_map[a])):
            annotation_data[ii] = annotation_map[a]

    return annotation_data


def get_inserts(read_tuples: List[Read], variant: Variant, window_size: int, sort_by: str='base') -> Dict:
    """A dictionary mapping insertions to reference positions.

    Ignores artificial haplotype read group.
    Relies on pysam's cigartuples structure see: http://pysam.readthedocs.io/en/latest/api.html
    Match, M -> 0
    Insert, I -> 1
    Deletion, D -> 2
    Ref Skip, N -> 3
    Soft Clip, S -> 4

    Arguments:
        read_tuples: list of aligned read tuples to find insertions within
        variant: the variant around which reads will load
        window_size: The size of the context window of genomic bases, i.e the width of the tensor
        sort_by: sort reads at the variant by base or refernce start

    Returns:
        insert_dict: a dict mapping read indices to max insertions at that point
    """
    insert_dict = {}

    idx_offset, ref_start, ref_end = get_variant_window(window_size, variant)

    for read in read_tuples:
        index_dif = ref_start - read.reference_start
        if abs(index_dif) >= window_size:
            continue

        if 'I' in read.cigar:
            cur_idx = 0
            for t in cigar_string_to_tuples(read.cigar):
                if t[0] == defines.CIGAR_CODE['I']:
                    insert_idx = cur_idx - index_dif
                    if insert_idx not in insert_dict:
                        insert_dict[insert_idx] = t[1]
                    elif insert_dict[insert_idx] < t[1]:
                        insert_dict[insert_idx] = t[1]

                if t[0] in CIGAR_CODES_TO_COUNT:
                    cur_idx += t[1]

    read_tuples.sort(key=lambda r: r.reference_start)
    if sort_by == 'base':
        read_tuples.sort(key=lambda r: get_base_to_sort_by(r, variant))

    return insert_dict


def get_base_to_sort_by(read: Read, variant: Variant) -> str:
    if len(read.seq) > 0:
        max_idx = len(read.seq)-1
    else:
        return 'Z'

    if variant.type == 'SNP':
        return read.seq[clamp((variant.pos-read.reference_start), 0, max_idx)]
    else:
        var_idx = (variant.pos-read.reference_start)+1
        cur_idx = 0
        for cur_op, length in cigar_string_to_tuples(read.cigar):
            cur_idx += length
            if cur_idx > var_idx:
                if cur_op == defines.CIGAR_CODE['M']:
                    return read.seq[clamp(var_idx, 0, max_idx)]
                else:
                    return defines.CODE2CIGAR[cur_op]
        return 'Y'


def cigar_string_to_tuples(cigar: str) -> List[Tuple]:
    if not cigar or len(cigar) == 0:
        return []
    parts = defines.CIGAR_REGEX.findall(cigar)
    # reverse order
    return [(defines.CIGAR2CODE[y], int(x)) for x,y in parts]


def get_variant_window(window_size: int, variant: Variant) -> Tuple:
    index_offset = (window_size//2)
    reference_start = variant.pos-index_offset
    reference_end = variant.pos + index_offset + (window_size % 2)
    return index_offset, reference_start, reference_end


def bool_from_java(val: str) -> bool:
    return val == 'true'


def clamp(n: int, minn: int, maxn: int) -> int:
    return max(min(maxn, n), minn)


def read_tuples_to_tensor(read_tuples: List[Read],
                          ref_start: int,
                          insert_dict: Dict,
                          tensor_type: str,
                          window_size: int,
                          read_limit: int,
                          base_quality_mode: str='phot') -> np.ndarray:
    """Create a read tensor based on a tensor channel map.

    Assumes read pairs have the same name.
    Only loads reads that might align inside the tensor.

    Arguments:
        read_tuples: list of reads to make into a tensor
        ref_start: the beginning of the window in reference coordinates
        insert_dict: a dict mapping read indices to max insertions at that point.
        tensor_type: The name for the type of tensor to make
        window_size: The size of the context window of genomic bases, i.e the width of the tensor
        read_limit: The maximum number of reads to encode in a tensor, i.e. the height of the tensor
        base_quality_mode: How to encode qualities in the tensor (phot, 1hot or phred)

    Returns:
        tensor: 3D read tensor.
    """
    channel_map = tensor_maps.get_tensor_channel_map_from_tensor_type(tensor_type)
    tensor = np.zeros(tensor_maps.tensor_shape_from_tensor_type(tensor_type, window_size, read_limit))

    if len(read_tuples) > read_limit:
        read_tuples_idx = np.random.choice(range(len(read_tuples)), size=read_limit, replace=False)
        read_tuples = [read_tuples[ii] for ii in read_tuples_idx]

    for j, read in enumerate(read_tuples):
        rseq, rqual = sequence_and_qualities_from_read(read, ref_start, insert_dict, window_size)
        flag_start = -1
        flag_end = 0

        for ii, b in enumerate(rseq):

            if ii == window_size:
                break

            if b == defines.SKIP_CHAR:
                continue
            elif flag_start == -1:
                flag_start = ii
            else:
                flag_end = ii

            if b in defines.INPUTS_INDEL:
                if b == defines.INDEL_CHAR:
                    if K.image_data_format() == 'channels_last':
                        tensor[j, ii, defines.INPUTS_INDEL[b]] = 1.0
                    else:
                        tensor[defines.INPUTS_INDEL[b], j, ii] = 1.0
                else:
                    hot_array = quality_from_mode(rqual[ii], b, defines.INPUTS_INDEL, base_quality_mode)
                    if K.image_data_format() == 'channels_last':
                        tensor[j, ii, :4] = hot_array
                    else:
                        tensor[:4, j, ii] = hot_array
            elif b in defines.AMBIGUITY_CODES:
                if K.image_data_format() == 'channels_last':
                    tensor[j, ii, :4] = defines.AMBIGUITY_CODES[b]
                else:
                    tensor[:4, j, ii] = defines.AMBIGUITY_CODES[b]
            else:
                raise ValueError('Unknown symbol in seq block:', b)

        if K.image_data_format() == 'channels_last':
            tensor[j, flag_start:flag_end, channel_map['flag_bit_4']] = 1.0 if read.reverse else 0.0
            tensor[j, flag_start:flag_end, channel_map['flag_bit_5']] = 1.0 if read.mate_reverse else 0.0
            tensor[j, flag_start:flag_end, channel_map['flag_bit_6']] = 1.0 if read.first else 0.0
            tensor[j, flag_start:flag_end, channel_map['flag_bit_7']] = 0.0 if read.first else 1.0
        else:
            tensor[channel_map['flag_bit_4'], j, flag_start:flag_end] = 1.0 if read.reverse else 0.0
            tensor[channel_map['flag_bit_5'], j, flag_start:flag_end] = 1.0 if read.mate_reverse else 0.0
            tensor[channel_map['flag_bit_6'], j, flag_start:flag_end] = 1.0 if read.first else 0.0
            tensor[channel_map['flag_bit_7'], j, flag_start:flag_end] = 0.0 if read.first else 1.0

        if 'mapping_quality' in channel_map:
            mq = float(read.mapping_quality) / defines.MAPPING_QUALITY_MAX
            if K.image_data_format() == 'channels_last':
                tensor[j, flag_start:flag_end, channel_map['mapping_quality']] = mq
            else:
                tensor[channel_map['mapping_quality'], j, flag_start:flag_end] = mq

    return tensor


def sequence_and_qualities_from_read(read: Read, ref_start: int, insert_dict: Dict, window_size: int) -> Tuple:
    cur_idx = 0
    my_indel_dict = {}
    no_qual_filler = 0

    index_dif = ref_start - read.reference_start
    for t in cigar_string_to_tuples(read.cigar):
        my_ref_idx = cur_idx - index_dif
        if t[0] == defines.CIGAR_CODE['I'] and my_ref_idx in insert_dict:
            my_indel_dict[my_ref_idx] = insert_dict[my_ref_idx] - t[1]
        elif t[0] == defines.CIGAR_CODE['D']:
            my_indel_dict[my_ref_idx] = t[1]
        if t[0] in CIGAR_CODES_TO_COUNT:
            cur_idx += t[1]

    for k in insert_dict.keys():
        if k not in my_indel_dict:
            my_indel_dict[k] = insert_dict[k]

    rseq = read.seq[:window_size]
    rqual = read.qual[:window_size]

    if index_dif > 0:
        rseq = rseq[index_dif:]
        rqual = rqual[index_dif:]
    elif index_dif < 0:
        rseq = defines.SKIP_CHAR * (-index_dif) + rseq
        rqual = [no_qual_filler]*(-index_dif) + rqual

    for j in sorted(my_indel_dict.keys(), key=int, reverse=True):
        if j < 1:
            rseq = (defines.INDEL_CHAR * my_indel_dict[j]) + rseq
            rqual = ([no_qual_filler]*my_indel_dict[j]) + rqual
        else:
            rseq = rseq[:j] + (defines.INDEL_CHAR * my_indel_dict[j]) + rseq[j:]
            rqual = rqual[:j] + ([no_qual_filler]*my_indel_dict[j]) + rqual[j:]

    return rseq, rqual


def reference_sequence_into_tensor(reference_seq: str,
                                   tensor: np.ndarray,
                                   insert_dict: Dict,
                                   window_size: int,
                                   read_limit: int):
    ref_offset = len(defines.INPUTS_INDEL)

    for ii in sorted(insert_dict.keys(), key=int, reverse=True):
        if ii < 0:
            reference_seq = defines.INDEL_CHAR*insert_dict[ii] + reference_seq
        else:
            reference_seq = reference_seq[:ii] + defines.INDEL_CHAR*insert_dict[ii] + reference_seq[ii:]

    for ii,b in enumerate(reference_seq):
        if ii == window_size:
            break

        if b in defines.INPUTS_INDEL:
            if K.image_data_format() == 'channels_last':
                tensor[:, ii, ref_offset+defines.INPUTS_INDEL[b]] = 1.0
            else:
                tensor[ref_offset+defines.INPUTS_INDEL[b], :, ii] = 1.0
        elif b in defines.AMBIGUITY_CODES:
            if K.image_data_format() == 'channels_last':
                tensor[:, ii, ref_offset:ref_offset+4] = np.tile(defines.AMBIGUITY_CODES[b], (read_limit, 1))
            else:
                tensor[ref_offset:ref_offset+4, :, ii] = np.transpose(
                    np.tile(defines.AMBIGUITY_CODES[b], (read_limit, 1)))


def base_quality_to_phred_array(base_quality: int, base: str, base_dict: Dict) -> np.ndarray:
    phred = np.zeros((4,))
    exponent = float(-base_quality) / 10.0
    p = 1.0-(10.0**exponent)  # Convert to probability
    not_p = (1.0-p) / 3.0  # Error could be any of the other 3 bases
    not_base_quality = -10 * np.log10(not_p)  # Back to Phred

    for b in base_dict.keys():
        if b == defines.INDEL_CHAR:
            continue
        elif b == base:
            phred[base_dict[b]] = base_quality
        else:
            phred[base_dict[b]] = not_base_quality
    return phred


def base_quality_to_p_hot_array(base_quality: int, base: str, base_dict: Dict) -> np.ndarray:
    not_p = not_p_lut[base_quality]
    phot = [not_p, not_p, not_p, not_p]
    phot[base_dict[base]] = p_lut[base_quality]

    return phot


def quality_from_mode(base_quality: int, base: str, base_dict: Dict, base_quality_mode: str) -> np.ndarray:
    if base_quality_mode == 'phot':
        return base_quality_to_p_hot_array(base_quality, base, base_dict)
    elif base_quality_mode == 'phred':
        return base_quality_to_phred_array(base_quality, base, base_dict)
    elif base_quality_mode == '1hot':
        one_hot = np.zeros((4,))
        one_hot[base_dict[base]] = 1.0
        return one_hot
    else:
        raise ValueError('Unknown base quality mode:', base_quality_mode)


def predictions_to_snp_scores(predictions: np.ndarray, eps: float=1e-7) -> np.ndarray:
    snp = predictions[:, defines.SNP_INDEL_LABELS['SNP']]
    not_snp = predictions[:, defines.SNP_INDEL_LABELS['NOT_SNP']]
    return np.log(eps + snp / (not_snp + eps))


def predictions_to_indel_scores(predictions: np.ndarray, eps: float=1e-7) -> np.ndarray:
    indel = predictions[:, defines.SNP_INDEL_LABELS['INDEL']]
    not_indel = predictions[:, defines.SNP_INDEL_LABELS['NOT_INDEL']]
    return np.log(eps + indel / (not_indel + eps))


def predictions_to_snp_indel_scores(predictions: np.ndarray) -> Tuple:
    snp_dict = predictions_to_snp_scores(predictions)
    indel_dict = predictions_to_indel_scores(predictions)
    return snp_dict, indel_dict


def _write_tensor_to_hd5(tensor: np.ndarray,
                         annotations: np.ndarray,
                         contig: str,
                         pos: str,
                         variant_type: str,
                         tensor_type: str,
                         annotation_set: str,
                         output_dir: str,) -> None:
    tensor_path = os.path.join(output_dir, 'inference_tensor_'+contig+pos+variant_type+defines.TENSOR_SUFFIX)
    if not os.path.exists(os.path.dirname(tensor_path)):
        os.makedirs(os.path.dirname(tensor_path))
    with h5py.File(tensor_path, 'w') as hf:
        hf.create_dataset(tensor_type, data=tensor, compression='gzip')
        hf.create_dataset(annotation_set, data=annotations, compression='gzip')
