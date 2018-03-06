# Imports
import numpy as np
from collections import Counter, defaultdict, namedtuple

# Keras Imports
import keras.backend as K

# Package Imports
from . import defines
from . import tensor_maps

READ_ELEMENTS = 8
Read = namedtuple("Read", "seq qual cigar reverse mate_reverse first mapping_quality reference_start")
Variant = namedtuple("Variant", "contig pos ref alt type")

CIGAR_CODES_TO_COUNT = [
    defines.CIGAR_CODE['M'], defines.CIGAR_CODE['I'], defines.CIGAR_CODE['S'], defines.CIGAR_CODE['D']
]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Inference ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def score_and_write_batch(args, model, file_out, fifo, batch_size, python_batch_size):
    '''Score a batch of variants with a CNN model. Write tab delimited temp file with scores.

    This function is tightly coupled with the CNNVariantScore.java
    It requires data written to the fifo in the order given by transferToPythonViaFifo

    Arguments
        args: Namespace with command line or configuration file set arguments
        model: a keras model
        file_out: The VCF file where variants scores are written
        fifo: The fifo opened by GATK Streaming executor
        batch_size: The total number of variants available in the fifo
        python_batch_size: the number of variants to process in each inference
        use_reads : Use read data from a BAM for 2D CNNs or just reference for 1D CNNs
    '''
    annotation_batch = []
    reference_batch = []
    variant_types = []
    variant_data = []
    read_batch = []

    for i in range(batch_size):
        fifo_line = fifo.readline()
        fifo_data = fifo_line.split(defines.SEPARATOR_CHAR)

        variant_data.append(fifo_data[0] + '\t' + fifo_data[1] + '\t' + fifo_data[2] + '\t' + fifo_data[3])
        reference_batch.append(reference_string_to_tensor(fifo_data[4]))
        annotation_batch.append(annotation_string_to_tensor(args, fifo_data[5]))
        variant_types.append(fifo_data[6])

        fidx = 7 # 7 Because above we parsed: contig pos ref alt reference_string annotation variant_type
        if args.tensor_name in defines.TENSOR_MAPS_2D and len(fifo_data) > fidx:
            read_tuples = []
            var = Variant(fifo_data[0], int(fifo_data[1]), fifo_data[2], fifo_data[3], fifo_data[6])
            while fidx+7 < len(fifo_data):
                read_tuples.append( Read(fifo_data[fidx],
                                         [int(q) for q in fifo_data[fidx+1].split(',')],
                                         fifo_data[fidx+2],
                                         bool_from_java(fifo_data[fidx+3]),
                                         bool_from_java(fifo_data[fidx+4]),
                                         bool_from_java(fifo_data[fidx+5]),
                                         int(fifo_data[fidx+6]),
                                         int(fifo_data[fidx+7])))
                fidx += READ_ELEMENTS
            _, ref_start, _ = get_variant_window(args, var)
            insert_dict = get_inserts(args, read_tuples, var)
            read_batch.append(read_tuples_to_read_tensor(args, read_tuples, ref_start, insert_dict))

    if args.tensor_name in defines.TENSOR_MAPS_1D:
        predictions = model.predict([np.array(reference_batch), np.array(annotation_batch)],
                                    batch_size=python_batch_size)
    elif args.tensor_name in defines.TENSOR_MAPS_2D:
        if len(read_batch) > 0:
            predictions = model.predict(
                {'read_tensor':np.array(read_batch), 'annotations':np.array(annotation_batch)},
                batch_size=python_batch_size)
    else:
        raise ValueError('Unknown tensor mapping.  Check architecture file.', args.tensor_name)

    indel_scores = predictions_to_indel_scores(predictions)
    snp_scores = predictions_to_snp_scores(predictions)

    for i in range(batch_size):
        if 'SNP' == variant_types[i]:
            file_out.write(variant_data[i]+'\t{0:.3f}'.format(snp_scores[i])+'\n')
        elif 'INDEL' == variant_types[i]:
            file_out.write(variant_data[i]+'\t{0:.3f}'.format(indel_scores[i])+'\n')
        else:
            file_out.write(variant_data[i]+'\t{0:.3f}'.format(snp_scores[i])+'\n')


def reference_string_to_tensor(reference):
    dna_data = np.zeros((len(reference), len(defines.DNA_SYMBOLS)))
    for i,b in enumerate(reference):
        if b in defines.DNA_SYMBOLS:
            dna_data[i, defines.DNA_SYMBOLS[b]] = 1.0
        elif b in defines.AMBIGUITY_CODES:
            dna_data[i] = defines.AMBIGUITY_CODES[b]
        else:
            raise ValueError('Error! Unknown code:', b)

    return dna_data


def annotation_string_to_tensor(args, annotation_string):
    name_val_pairs = annotation_string.split(';')
    name_val_arrays = [p.split('=') for p in name_val_pairs]
    annotation_map = {str(p[0]).strip() : p[1] for p in name_val_arrays if len(p) > 1}
    annotation_data = np.zeros(( len(defines.ANNOTATIONS[args.annotation_set]),))

    for i,a in enumerate(defines.ANNOTATIONS[args.annotation_set]):
        if a in annotation_map:
            annotation_data[i] = annotation_map[a]

    return annotation_data


def get_inserts(args, read_tuples, variant, sort_by='base'):
    '''A dictionary mapping insertions to reference positions.

    Ignores artificial haplotype read group.
    Relies on pysam's cigartuples structure see: http://pysam.readthedocs.io/en/latest/api.html
    Match, M -> 0
    Insert, I -> 1
    Deletion, D -> 2
    Ref Skip, N -> 3
    Soft Clip, S -> 4

    Arguments:
        args.read_limit: maximum number of reads to return
        samfile: the BAM (or BAMout) file
        variant: the variant around which reads will load

    Returns:
        insert_dict: a dict mapping read indices to max insertions at that point
    '''
    insert_dict = {}

    idx_offset, ref_start, ref_end = get_variant_window(args, variant)

    for read in read_tuples:
        index_dif = ref_start - read.reference_start
        if abs(index_dif) >= args.window_size:
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


    read_tuples.sort(key=lambda read: read.reference_start)
    if sort_by == 'base':
        read_tuples.sort(key=lambda read: get_base_to_sort_by(read, variant))

    return insert_dict


def get_base_to_sort_by(read, variant):
    if len(read.seq) > 0:
        max_idx = len(read.seq)-1
    else:
        return 'Z'

    if variant.type != 'INDEL':
        return read.seq[clamp((variant.pos-read.reference_start)-1, 0, max_idx)]
    else:
        var_idx = variant.pos-read.reference_start
        cur_idx = 0
        for cur_op, length in cigar_string_to_tuples(read.cigar):
            cur_idx += length
            if cur_idx > var_idx:
                if cur_op == defines.CIGAR_CODE['M']:
                    return read.seq[clamp(var_idx, 0, max_idx)]
                else:
                    return defines.CODE2CIGAR[cur_op]
        return 'Y'


def cigar_string_to_tuples(cigar):
    if not cigar or len(cigar) == 0:
        return []
    parts = defines.CIGAR_REGEX.findall(cigar)
    # reverse order
    return [(defines.CIGAR2CODE[y], int(x)) for x,y in parts]


def get_variant_window(args, variant):
    index_offset = (args.window_size//2)
    reference_start = variant.pos-(index_offset+1)
    reference_end = variant.pos+index_offset
    return index_offset, reference_start, reference_end

def bool_from_java(val):
    return val == 'true'

def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)


def read_tuples_to_read_tensor(args, read_tuples, ref_start, insert_dict):
    '''Create a read tensor based on a tensor channel map.

    Assumes read pairs have the same name.
    Only loads reads that might align inside the tensor.

    Arguments:
        args.read_limit: maximum number of reads to return
        good_reads: list of reads to make arrays from
        ref_start: the beginning of the window in reference coordinates
        insert_dict: a dict mapping read indices to max insertions at that point.

    Returns:
        tensor: 3D read tensor.
    '''
    channel_map = tensor_maps.get_tensor_channel_map_from_args(args)
    tensor = np.zeros( tensor_maps.tensor_shape_from_args(args) )

    for j,read in enumerate(read_tuples):
        rseq, rqual = sequence_and_qualities_from_read(args, read, ref_start, insert_dict)
        flag_start = -1
        flag_end = 0

        for i,b in enumerate(rseq):

            if i == args.window_size:
                break
            if b == defines.SKIP_CHAR:
                continue
            elif flag_start == -1:
                flag_start = i
            else:
                flag_end = i

            if b in args.input_symbols:
                if b == defines.INDEL_CHAR:
                    if K.image_data_format() == 'channels_last':
                        tensor[j, i, args.input_symbols[b]] = 1.0
                    else:
                        tensor[args.input_symbols[b], j, i] = 1.0
                else:
                    hot_array = quality_from_mode(args, rqual[i], b, args.input_symbols)
                    if K.image_data_format() == 'channels_last':
                        tensor[j, i, :4] = hot_array
                    else:
                        tensor[:4, j, i] = hot_array
            elif b in defines.AMBIGUITY_CODES:
                if K.image_data_format() == 'channels_last':
                    tensor[j, i, :4] = defines.AMBIGUITY_CODES[b]
                else:
                    tensor[:4, j, i] = defines.AMBIGUITY_CODES[b]
            else:
                raise ValueError('Unknown symbol in seq block:', b)

        if K.image_data_format() == 'channels_last':
            tensor[j, flag_start:flag_end, channel_map['flag_bit_4']] = 1.0 if read.reverse else 0.0
            tensor[j, flag_start:flag_end, channel_map['flag_bit_5']] = 1.0 if read.mate_reverse else 0.0
            tensor[j, flag_start:flag_end, channel_map['flag_bit_6']] = 1.0 if read.first else 0.0
            tensor[j, flag_start:flag_end, channel_map['flag_bit_7']] = 0.0 if read.first else 1.0
        else:
            tensor[channel_map['flag_bit_4'], j,  flag_start:flag_end] = 1.0 if read.reverse else 0.0
            tensor[channel_map['flag_bit_5'], j,  flag_start:flag_end] = 1.0 if read.mate_reverse else 0.0
            tensor[channel_map['flag_bit_6'], j, flag_start:flag_end] = 1.0 if read.first else 0.0
            tensor[channel_map['flag_bit_7'], j, flag_start:flag_end] = 0.0 if read.first else 1.0

        if 'mapping_quality' in channel_map:
            mq = float(read.mapping_quality) / defines.MAPPING_QUALITY_MAX
            if K.image_data_format() == 'channels_last':
                tensor[j, flag_start:flag_end, channel_map['mapping_quality']] = mq
            else:
                tensor[channel_map['mapping_quality'], j, flag_start:flag_end] = mq

    return tensor


def sequence_and_qualities_from_read(args, read, ref_start, insert_dict):
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

    rseq = read.seq[:args.window_size]
    rqual = read.qual[:args.window_size]

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


def base_quality_to_phred_array(base_quality, base, base_dict):
    phred = np.zeros((4,))
    exponent = float(-base_quality) / 10.0
    p = 1.0-(10.0**exponent) # Convert to probability
    not_p = (1.0-p) / 3.0 # Error could be any of the other 3 bases
    not_base_quality = -10 * np.log10(not_p) # Back to Phred

    for b in base_dict.keys():
        if b == defines.INDEL_CHAR:
            continue
        elif b == base:
            phred[base_dict[b]] = base_quality
        else:
            phred[base_dict[b]] = not_base_quality
    return phred


def base_quality_to_p_hot_array(base_quality, base, base_dict):
    phot = np.zeros((4,))
    exponent = float(-base_quality) / 10.0
    p = 1.0-(10.0**exponent)
    not_p = (1.0-p)/3.0

    for b in base_dict.keys():
        if b == base:
            phot[base_dict[b]] = p
        elif b == defines.INDEL_CHAR:
            continue
        else:
            phot[base_dict[b]] = not_p

    return phot


def quality_from_mode(args, base_quality, base, base_dict):
    if args.base_quality_mode == 'phot':
        return base_quality_to_p_hot_array(base_quality, base, base_dict)
    elif args.base_quality_mode == 'phred':
        return base_quality_to_phred_array(base_quality, base, base_dict)
    elif args.base_quality_mode == '1hot':
        one_hot = np.zeros((4,))
        one_hot[base_dict[base]] = 1.0
        return one_hot
    else:
        raise ValueError('Unknown base quality mode:', args.base_quality_mode)


def predictions_to_snp_scores(predictions, eps=1e-7):
    snp = predictions[:, defines.SNP_INDEL_LABELS['SNP']]
    not_snp = predictions[:, defines.SNP_INDEL_LABELS['NOT_SNP']]
    return np.log(eps + snp / (not_snp + eps))


def predictions_to_indel_scores(predictions, eps=1e-7):
    indel = predictions[:, defines.SNP_INDEL_LABELS['INDEL']]
    not_indel = predictions[:, defines.SNP_INDEL_LABELS['NOT_INDEL']]
    return np.log(eps + indel / (not_indel + eps))


def predictions_to_snp_indel_scores(predictions):
    snp_dict = predictions_to_snp_scores(predictions)
    indel_dict = predictions_to_indel_scores(predictions)
    return snp_dict, indel_dict


