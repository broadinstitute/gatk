import numpy as np
import abc
from pysam import VariantRecord
from typing import Dict, List
from scorevariants.utilities import variant_is_snp
from scorevariants.random_generator import Random
from enum import Enum


BEST_PRACTICES = ["MQ", "DP", "SOR", "FS", "QD", "MQRankSum", "ReadPosRankSum"]

BASE_MAP = {"A": 0, "C": 1, "G": 2, "T": 3}


# copied from gatk's vqsr_cnn.defines
CODE2CIGAR = "MIDNSHP=XB"
CIGAR_CODE = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4}
CODES_TO_COUNT = {
    CIGAR_CODE["M"],
    CIGAR_CODE["I"],
    CIGAR_CODE["S"],
    CIGAR_CODE["D"],
}
INPUTS_INDEL = {"A": 0, "C": 1, "G": 2, "T": 3, "*": 4}
AMBIGUITY_CODES = {
    "K": [0, 0, 0.5, 0.5],
    "M": [0.5, 0.5, 0, 0],
    "R": [0.5, 0, 0, 0.5],
    "Y": [0, 0.5, 0.5, 0],
    "S": [0, 0.5, 0, 0.5],
    "W": [0.5, 0, 0.5, 0],
    "B": [0, 0.333, 0.333, 0.334],
    "V": [0.333, 0.333, 0, 0.334],
    "H": [0.333, 0.333, 0.334, 0],
    "D": [0.333, 0, 0.333, 0.334],
    "X": [0.25, 0.25, 0.25, 0.25],
    "N": [0.25, 0.25, 0.25, 0.25],
}

AMBIGUITY_CODES_convertIUPACtoN = {
    "K": [0.25, 0.25, 0.25, 0.25],
    "M": [0.25, 0.25, 0.25, 0.25],
    "R": [0.25, 0.25, 0.25, 0.25],
    "Y": [0.25, 0.25, 0.25, 0.25],
    "S": [0.25, 0.25, 0.25, 0.25],
    "W": [0.25, 0.25, 0.25, 0.25],
    "B": [0.25, 0.25, 0.25, 0.25],
    "V": [0.25, 0.25, 0.25, 0.25],
    "H": [0.25, 0.25, 0.25, 0.25],
    "D": [0.25, 0.25, 0.25, 0.25],
    "X": [0.25, 0.25, 0.25, 0.25],
    "N": [0.25, 0.25, 0.25, 0.25],
}

SKIP_CHAR = "~"
INDEL_CHAR = "*"
MAPPING_QUALITY_MAX = (
    60.0  # Mapping qualities from BWA are typically capped at 60
)
# end copied material

GATK_RANDOM_SEED = np.int64(47382911)
Random.multiplier = np.int64(0x5DEECE66D)
Random.addend = np.int64(0xB)
Random.mask = np.int64((1 << 48) - 1)
mRandomGenerator = Random(GATK_RANDOM_SEED)

def get_start(read):
    if read.is_unmapped:
        return -1
    return read.reference_start - read.query_alignment_start


def clamp(n: int, minn: int, maxn: int) -> int:
    # this is copied/edited from vqsr_cnn/training.py
    return max(min(maxn, n), minn)

def get_base_to_sort_by(read, variant):
    if len(read.query_sequence) > 0:
        max_idx = len(read.query_sequence)-1
    else:
        return 'Z'

    if variant_is_snp(variant):
        return read.query_sequence[clamp((variant.pos-get_start(read))-1, 0, max_idx)]
    else:
        var_idx = (variant.pos-get_start(read))
        cur_idx = 0
        for cur_op, length in read.cigartuples or []:
            cur_idx += length
            if cur_idx > var_idx:
                if cur_op == CIGAR_CODE['M']:
                    return read.seq[clamp(var_idx, 0, max_idx)]
                else:
                    return CODE2CIGAR[cur_op]
        return 'Y'



def base_quality_to_p_hot_array(
    base_quality: int, base: str, base_dict: Dict
) -> np.ndarray:
    # this is copied/edited from vqsr_cnn/inference.py
    phot = np.zeros((4,))
    exponent = float(-base_quality) / 10.0
    p = 1.0 - (10.0 ** exponent)  # Convert to probability
    not_p = (1.0 - p) / 3.0  # Error could be any of the other 3 bases

    for b in base_dict.keys():
        if b == INDEL_CHAR:
            continue
        elif b == base:
            phot[base_dict[b]] = p
        else:
            phot[base_dict[b]] = not_p
    return phot


class Interval:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class AlignmentReads:
    def __init__(self, reads, insertions):
        self.reads = reads
        self.insertions = insertions


class VariantLabel(Enum):
    NONE = -1
    NOT_SNP = 0
    NOT_INDEL = 1
    SNP = 2
    INDEL = 3



class VariantType(Enum):
    OTHER = -1
    SNP = 0
    INDEL = 1


variant_label_map = {i.name: i.value for i in VariantLabel}
variant_type_map = {i.name: i.value for i in VariantType}


def get_variant_window(window_size, position, ref_size, offset=0):
    index_offset = window_size // 2
    reference_start = position - index_offset
    reference_end = position + index_offset + (window_size % 2)
    reference_end += min(max(-reference_start, 0), ref_size - 1)
    return reference_start + offset, reference_end + offset


def encode_annotations(variant: VariantRecord, annotation_list: List[str]):
    encoding = np.zeros(len(annotation_list))
    for i, key in enumerate(annotation_list):
        encoding[i] = variant.info.get(key, 0.0)
    return encoding


class Encoder:
    @abc.abstractmethod
    def __call__(self, sample):
        raise NotImplementedError("Encoding method not implemented.")


class AlignmentReads:
    def __init__(self, reads, insertions):
        self.reads = reads
        self.insertions = insertions

    def __iter__(self):
        for read in self.reads:
            yield read


class _DefaultMapEncoder(Encoder):
    def __init__(self, map_=None, default=-1):
        if map_ is None:
            map_ = self._get_default_map()
        self.default = default
        self.map_ = map_
        self.inverse_map_ = {v: k for k, v, in self.map_.items()}

    @abc.abstractmethod
    def _get_default_map(self):
        raise NotImplementedError("Default map Not implemented.")

    def __call__(self, item):
        return self.map_.get(item, self.default)

    def inverse(self, key):
        return self.inverse_map_.get(key)


class VariantLabelEncoder(_DefaultMapEncoder):
    @staticmethod
    def _get_default_map():
        return variant_label_map


class VariantTypeEncoder(_DefaultMapEncoder):
    @staticmethod
    def _get_default_map():
        return variant_type_map


class ReadTensorEncoder(Encoder):
    
    def __init__(
        self, alignment_file, reference, window=128, offset=0, read_limit=128
    ):
        """
        Encodes a pileup as a tensor of a window centered
        around a position.

        Args:
            alignment_file: pysam.AlignmentFile containing alignments
                for the pileup
            reference: dict of strings containing the reference genome
                for the alignment_file.
            window:
                Optional int. Determines with width of the tensor.
            offset:
                Optional int. When offset is 0, the variant POS
                is located at middle. negative shifts it left
                and positive shifts it right, with respect
                to the output tensor.
            read_limit:
                Optional int. Maximum number of reads to include in the
                encoding.
        """
        self.alignment_file = alignment_file
        self.reference = reference
        for chrom in reference:
            reference[chrom].seq = reference[chrom].seq.upper()
        self.window_size = window
        self.offset = offset
        self.read_limit = read_limit
        self.channel_map = self._create_channel_map()
        self.n_channels = max(self.channel_map.values()) + 1
        self.encode_read_base_and_quality = base_quality_to_p_hot_array
        self._read_base_channels = np.zeros((4,), dtype=int)
        for base, pos in BASE_MAP.items():
            self._read_base_channels[pos] = self.channel_map[f"read_{base}"]

    def __call__(self, variant: VariantRecord):
    
        """
        variant should be 0-based index
        """
        interval = Interval(
            *get_variant_window(
                self.window_size,
                variant.pos - 1,
                len(variant.ref),
                self.offset,
            )
        )

        tensor = self.initialize_representation()
        reads = self.get_reads(variant, interval)
        reference = self.get_reference(
            self.reference,
            variant.contig,
            interval,
            reads.insertions,
        )
        # the interval potentially needs to be shifted if the interval starts
        # to the left of the start of the reference sequence
        self.insert_reference(reference, tensor)
        self.insert_reads_and_flags(reads, interval, tensor)
        return tensor

    def initialize_representation(self):
        shape = (self.read_limit, self.window_size, self.n_channels)
        tensor = np.zeros(shape)
        return tensor

    def _add_insertions(self, reference_sequence, insertions, fill, beginning, padded):
        for i in sorted(insertions, reverse=True):
            # if the insertion happens before the intervals opens
            insertion_sequence = fill * insertions[i]
            if i < beginning:
                reference_sequence = insertion_sequence + reference_sequence
            elif (not padded) or i <= len(reference_sequence):
                reference_sequence = (
                    reference_sequence[:i]
                    + insertion_sequence
                    + reference_sequence[i:]
                )
            elif i < self.window_size:
                reference_sequence = (
                    reference_sequence
                    + '\00' * (i-len(reference_sequence))
                    + insertion_sequence
                )
        return reference_sequence

    def get_reference(self, reference, contig, interval, insertions):
        # this will get changed, so this relies on slice-by-copy
        reference_sequence = reference[contig][max(0, interval.start) : interval.end + (1 if interval.start < 0 else 0)]
        # this is reversed so the indices of the insertions
        # do not change as we modify the reference
        reference_sequence = self._add_insertions(
            reference_sequence, insertions, INDEL_CHAR, 0, padded=True
        )
        return reference_sequence

    def insert_reference(self, reference, tensor):
        base_channels = [
            self.channel_map[f"reference_{base}"] for base in "ACGT"
        ]

        for i, base in enumerate(reference):
            if i == self.window_size:
                break
            if base in INPUTS_INDEL:
                tensor[:, i, self.channel_map[f"reference_{base}"]] = 1.0
            elif base in AMBIGUITY_CODES_convertIUPACtoN:
                ambiguous_vector = np.tile(
                    AMBIGUITY_CODES_convertIUPACtoN[base], (self.read_limit, 1)
                )
                tensor[:, i, base_channels] = ambiguous_vector

    @staticmethod
    def _create_channel_map():
        tensor_map = {}
        for k in INPUTS_INDEL.keys():
            tensor_map["read_" + k] = INPUTS_INDEL[k]
        for k in INPUTS_INDEL.keys():
            tensor_map["reference_" + k] = len(INPUTS_INDEL) + INPUTS_INDEL[k]
        tensor_map["flag_bit_4"] = 10
        tensor_map["flag_bit_5"] = 11
        tensor_map["flag_bit_6"] = 12
        tensor_map["flag_bit_7"] = 13
        tensor_map["mapping_quality"] = 14
        return tensor_map

    def _get_sequence_specific_indel_dict(self, read, offset, insert_dict):
        my_indel_dict = dict()
        cur_idx = 0
        # get sequence specific indel dict
        for (code, length) in read.cigartuples or []:
            my_ref_idx = cur_idx - offset

            # this read's bases will be present in the tensor,
            # so do not add *'s for the insertions covered by this read
            if code == CIGAR_CODE["I"] and my_ref_idx in insert_dict:
                my_indel_dict[my_ref_idx] = insert_dict[my_ref_idx] - length
            # add *'s for bases any deletions
            elif code == CIGAR_CODE["D"]:
                my_indel_dict[my_ref_idx] = length

            if code in CODES_TO_COUNT:
                cur_idx += length

        # only add additional insertions from all reads if
        # they are not in this reads insertions
        for key in insert_dict:
            if key not in my_indel_dict:
                my_indel_dict[key] = insert_dict[key]

        return my_indel_dict

    def _get_sequence_and_quality(self, read, insertions, interval):
        ref_start = interval.start

        no_qual_filler = 0

        index_dif = ref_start - get_start(read)

        my_indel_dict = self._get_sequence_specific_indel_dict(
            read, index_dif, insertions
        )

        sequence = read.query_sequence[: self.window_size]
        quality = read.query_qualities[: self.window_size].tolist()

        # pad the sequence and quality scores for the window
        # chop off the left side of the sequence if it is outside
        # the window
        if index_dif > 0:
            sequence = sequence[index_dif:]
            quality = quality[index_dif:]
        # pad the left side if needed
        
        elif index_dif < 0:
            sequence = SKIP_CHAR * (-index_dif) + sequence
            quality = [no_qual_filler] * (-index_dif) + quality

        # pad the internal bits where there are insertions in the pileup
        # that are not present in this read
        sequence = self._add_insertions(sequence, my_indel_dict, INDEL_CHAR, 1, padded=False)
        quality = self._add_insertions(
            quality, my_indel_dict, [no_qual_filler], 1, padded=False
        )

        return sequence, quality

    def insert_reads_and_flags(self, alignment_reads, interval, tensor):
        for alignment_number, read in enumerate(alignment_reads):
            sequence, quality = self._get_sequence_and_quality(
                read, alignment_reads.insertions, interval
            )
            self._fill_rows(alignment_number, read, sequence, quality, tensor)

    def _fill_rows(self, alignment_number, read, sequence, quality, tensor):
        flag_start = -1
        flag_end = 0
        for position, (base, qual) in enumerate(zip(sequence, quality)):
            # this is probably not great, for longer reads we will prep a bunch
            # a bunch of extra sequence
            if position == self.window_size:
                break
            # skip beginning of the sequence (do not fill)
            if base == SKIP_CHAR:
                continue
            # if this is the first base, indicate we have started
            elif flag_start == -1:
                flag_start = position
            else:
                flag_end = position

            if base in INPUTS_INDEL:
                if base == INDEL_CHAR:
                    channel = self.channel_map[f"read_{base}"]
                    tensor[alignment_number, position, channel] = 1.0
                else:
                    base_array = self.encode_read_base_and_quality(
                        qual,
                        base,
                        BASE_MAP,
                    )
                    # the :4 indexing could change if the tensor channels change
                    tensor[
                        alignment_number, position, self._read_base_channels
                    ] = base_array

            elif base in AMBIGUITY_CODES:
                tensor[
                    alignment_number, position, self._read_base_channels
                ] = AMBIGUITY_CODES[base]

        # this matches vqsr_cnn/inference.py, the training.py had a more
        # general/extensible approach, but following this for sake of matching
        # the inference code
        channel_map = self.channel_map
        tensor[
            alignment_number, flag_start:flag_end, channel_map["flag_bit_4"]
        ] = (1.0 if read.is_reverse else 0.0)
        tensor[
            alignment_number, flag_start:flag_end, channel_map["flag_bit_5"]
        ] = (1.0 if read.mate_is_reverse else 0.0)
        tensor[
            alignment_number, flag_start:flag_end, channel_map["flag_bit_6"]
        ] = (1.0 if read.is_read1 else 0.0)
        tensor[
            alignment_number, flag_start:flag_end, channel_map["flag_bit_7"]
        ] = (1.0 if read.is_read2 else 0.0)

        mq = float(read.mapping_quality) / MAPPING_QUALITY_MAX
        tensor[
            alignment_number,
            flag_start:flag_end,
            channel_map["mapping_quality"],
        ] = mq

    def filter_read(self, read, window):
        if not read:
            True

        if (not read.is_unmapped) and (read.query_alignment_start < 0):
            True

        if (not read.is_unmapped) and (read.query_alignment_length + 1 < 0):
            True

        if len(read.query_qualities) != len(read.query_sequence):
            True

        if (not read.is_unmapped) and (len(read.cigarstring) != len(read.query_sequence)):
            True

        if len(read.query_sequence) <= 0:
            True

        if read.cigarstring and "N" in read.cigarstring:
            True

        # filters HaplotypeCaller artificial haplotypes
        try:
            read_group = read.get_tag("RG")
        except KeyError:
            True
        if "artificial" in read_group.lower():
            return True

        return False

    @staticmethod
    def get_insertions(reads, interval, window_size):
        insert_dict = {}
        for read in reads:
            index_dif = interval.start - get_start(read)
            if abs(index_dif) >= window_size:
                continue

            if not read.cigarstring:
                continue

            if "I" in read.cigarstring:
                cur_idx = 0
                for t in read.cigartuples:
                    if t[0] == CIGAR_CODE["I"]:
                        # this is the index with respect to the tensor
                        insert_idx = cur_idx - index_dif
                        if insert_idx not in insert_dict:
                            insert_dict[insert_idx] = t[1]
                        elif insert_dict[insert_idx] < t[1]:
                            insert_dict[insert_idx] = t[1]

                    if t[0] in CODES_TO_COUNT:
                        cur_idx += t[1]
        return insert_dict

    def _sort_reads(self, reads, variant, sort_by="base"):
        reads.sort(key=lambda x: get_start(x))
        if sort_by == "base":
            reads.sort(key=lambda read: get_base_to_sort_by(read, variant))
            
    def get_reads(self, variant, interval):
        reads = []
        insertions = dict()

        start = variant.start
        stop = variant.stop
        alignment_number = 0
        for read in self.alignment_file.fetch(
            variant.contig,
            start,
            stop,
            multiple_iterators=False,
        ):
            if self.filter_read(read, interval):
                continue
            alignment_number += 1
            if alignment_number <= self.read_limit:
                reads.append(read)
            else:
                randomSlot = mRandomGenerator.nextInt(alignment_number)
                if randomSlot < self.read_limit:
                    reads[randomSlot] = read
        insertions = self.get_insertions(reads, interval, self.window_size)

        self._sort_reads(reads, variant, sort_by="base")
        
        return AlignmentReads(reads, insertions)


class AnnotationEncoder(Encoder):
    def __init__(self, annotation_list=None):
        if annotation_list is None:
            annotation_list = BEST_PRACTICES
        self.annotation_list = annotation_list

    def __call__(self, variant: VariantRecord):
        return encode_annotations(variant, self.annotation_list)


class DNAEncoder(Encoder):
    def __init__(self, symbol_map=None):
        if symbol_map is None:
            symbol_map = BASE_MAP
        self.symbol_map = symbol_map

    def __call__(self, base):
        return self.symbol_map.get(base)


class DNAVectorEncoder(Encoder):
    def __init__(self, ambiguity_getter="N"):
        self.symbol_map = BASE_MAP
        if ambiguity_getter == "N":
            ambiguity_getter = lambda x: np.array([0.25, 0.25, 0.25, 0.25])
        elif ambiguity_getter == "code":
            ambiguity_getter = lambda x: AMBIGUITY_CODES.get(x)
        elif not callable(ambiguity_getter):
            raise ValueError(
                f"Invalid choice for `ambiguity_getter`: {ambiguity_getter}"
            )
        self.ambiguity_getter = ambiguity_getter

    def __call__(self, base):
        encoding = np.zeros(4)
        if base in self.symbol_map:
            encoding[self.symbol_map[base]] = 1.0
        elif base in AMBIGUITY_CODES:
            encoding[:] = self.ambiguity_getter(base)
        return encoding


class _EdgeStrategy:
    @abc.abstractmethod
    def get_position(*args, **kwargs):
        raise NotImplementedError()

    @abc.abstractmethod
    def get_right_side(*args, **kwargs):
        raise NotImplementedError()


class _LeftStrategy(_EdgeStrategy):
    """
    Strategy for determining reference window used by gatk
    """

    @staticmethod
    def get_position(i, reference_start, interval_start, *args, **kwargs):
        return i - interval_start

    @staticmethod
    def get_interval_end(reference_end, contig, reference_start, **kwargs):
        return min(reference_end, len(contig))  + (
            1 if reference_start < 0 else 0
        )


class _PadStrategy(_EdgeStrategy):
    @staticmethod
    def get_position(i, reference_start, interval_start, *args, **kwargs):
        return i - reference_start

    @staticmethod
    def get_interval_end(reference_end, contig, *args, **kwargs):
        return min(reference_end, len(contig))


class ReferenceEncoder(Encoder):
    def __init__(
        self, window=128, offset=0, base_encoder=None, handle_start="left"
    ):
        """
        Encodes a string as a tensor of a window centered
        around a position.

        Args:
            window:
                Optional int. Determines with width of the tensor.
            offset:
                Optional int. When offset is 0, the variant POS
                is located at middle. negative shifts it left
                and positive shifts it right, with respect
                to the output tensor.
            base_encoder:
                Optional Encoder. Encodes a base from the
                string as an int representing the channel
                in the output tensor.
            handle_start:
                Optional str or callable. Strategy for handling
                the start of the tensor when the window extends
                beyond the left side of the reference.

                Options:
                * 'left': the first base of the reference is
                    the first base in the tensor
                * 'pad': the reference's start is padded such
                    that the variant POS remains unchanged
                    with respect to the output tensor.

        """
        self.window_size = window
        self.offset = offset
        if base_encoder is None:
            base_encoder = DNAVectorEncoder()
        self.base_encoder = base_encoder

        if handle_start == "left":
            self._edge_strategy = _LeftStrategy()
        elif handle_start == "pad":
            self._edge_strategy = _PadStrategy()
        elif isinstance(handle_start, _EdgeStrategy):
            self._edge_strategy = handle_start
        else:
            raise ValueError(
                f"Invalid choice for `handle_start`: {handle_start}"
            )

    def __call__(self, reference: Dict, variant: VariantRecord):
        reference_start, reference_end = get_variant_window(
            self.window_size,
            variant.pos - 1,
            len(variant.ref),
            self.offset,
        )
        representation = np.zeros((self.window_size, 4))

        contig = reference[variant.contig]
        interval_start = max(reference_start, 0)
        interval_end = self._edge_strategy.get_interval_end(
            reference_end, contig, reference_start
        )
        for i in range(interval_start, interval_end):
            base = contig[i].upper()
            base_values = self.base_encoder(base)
            position = self._edge_strategy.get_position(
                i, reference_start, interval_start
            )
            representation[position, :] = base_values
        return representation
