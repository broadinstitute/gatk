# Imports
import os
import sys
import vcf
import math
import h5py
import time
import pysam
import vqsr_cnn
import numpy as np
from Bio import Seq, SeqIO
from collections import Counter

# Keras Imports
import keras.backend as K

def run():
    args = vqsr_cnn.parse_args()
    if 'write_reference_and_annotation_tensors' == args.mode:
        write_reference_and_annotation_tensors(args)
    elif 'write_read_and_annotation_tensors' == args.mode:
        write_read_and_annotation_tensors(args)
    elif 'train_on_reference_tensors_and_annotations' == args.mode:
        train_on_reference_tensors_and_annotations(args)
    elif 'train_on_read_tensors_and_annotations' == args.mode:
        train_on_read_tensors_and_annotations(args)
    elif 'train_tiny_model_on_read_tensors_and_annotations' == args.mode:
        train_tiny_model_on_read_tensors_and_annotations(args)
    elif 'train_small_model_on_read_tensors_and_annotations' == args.mode:
        train_small_model_on_read_tensors_and_annotations(args)
    else:
        raise ValueError('Unknown training mode:', args.mode)


def write_reference_and_annotation_tensors(args, include_dna=True, include_annotations=True):
    if not args.tensor_name in vqsr_cnn.TENSOR_MAPS_1D:
        raise ValueError('Unknown tensor name:', args.tensor_name, '1d maps must be in:', str(vqsr_cnn.TENSOR_MAPS_1D))

    record_dict = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, "fasta"))
    if os.path.splitext(args.input_vcf)[-1].lower() == '.gz':
        vcf_reader = vcf.Reader(open(args.input_vcf, 'rb'))
    else:
        vcf_reader = vcf.Reader(open(args.input_vcf, 'r'))

    if os.path.splitext(args.train_vcf)[-1].lower() == '.gz':
        vcf_ram = vcf.Reader(open(args.train_vcf, 'rb'))
    else:
        vcf_ram = vcf.Reader(open(args.train_vcf, 'r'))

    bed_dict = bed_file_to_dict(args.bed_file)
    stats = Counter()

    if args.chrom:
        variants  = vcf_reader.fetch(args.chrom, args.start_pos, args.end_pos)
    else:
        variants = vcf_reader

    for variant in variants:
        for allele_idx, allele in enumerate(variant.ALT):
            idx_offset, ref_start, ref_end = get_variant_window(args, variant)
            contig = record_dict[variant.CHROM]
            record = contig[variant.POS-idx_offset: variant.POS+idx_offset]

            cur_label_key = get_true_label(allele, variant, bed_dict, vcf_ram, stats)
            if not cur_label_key or downsample(args, cur_label_key, stats):
                continue

            if include_annotations:
                if all(map(
                        lambda x: x not in variant.INFO and x not in variant.FORMAT and x != "QUAL", args.annotations)):
                    stats['Missing ALL annotations'] += 1
                    continue # Require at least 1 annotation...
                annotation_data = get_annotation_data(args, variant, stats)

            if include_dna:
                dna_data = np.zeros( (args.window_size, len(vqsr_cnn.DNA_SYMBOLS)) )
                for i,b in enumerate(record.seq):
                    if b in vqsr_cnn.DNA_SYMBOLS:
                        dna_data[i, vqsr_cnn.DNA_SYMBOLS[b]] = 1.0
                    elif b in vqsr_cnn.AMBIGUITY_CODES:
                        dna_data[i] = vqsr_cnn.AMBIGUITY_CODES[b]
                    else:
                        raise ValueError('Error! Unknown code:', b)

            tp = get_path_to_train_valid_or_test(args.data_dir)
            tp += cur_label_key +'/'+ plain_name(args.input_vcf) +'_'+ plain_name(args.train_vcf)
            tp += '_allele_' + str(allele_idx) +'-'+ variant.CHROM +'_'+ str(variant.POS) + vqsr_cnn.TENSOR_SUFFIX
            if not os.path.exists(os.path.dirname(tp)):
                os.makedirs(os.path.dirname(tp))

            with h5py.File(tp, 'w') as hf:
                if include_annotations:
                    hf.create_dataset(args.annotation_set, data=annotation_data, compression='gzip')
                if include_dna:
                    hf.create_dataset(args.tensor_name, data=dna_data, compression='gzip')

            stats[cur_label_key] += 1
            stats['count'] += 1
            if stats['count']%500==0:
                print('Wrote', stats['count'], 'out of:', args.samples, 'Last variant:', variant)
        if args.samples <= stats['count']:
            break

    print('Done Writing 1D Tensors. Tensor Map:', args.tensor_name, ' Annotation set:', args.annotation_set)
    for k in stats.keys():
        print(k, ' has:', stats[k])



def write_read_and_annotation_tensors(args, include_annotations=True, pileup=False):
    '''Create tensors structured as tensor map of reads organized by labels in the data directory.

    Defines true variants as those in the args.train_vcf, defines false variants as
    those called in args.input_vcf and in the args.bed_file high confidence intervals,
    but not in args.train_vcf.

    Arguments
        args.data_dir: directory where tensors will live. Created here and filled with
            subdirectories of test, valid and train, each containing
            subdirectories for each label with tensors stored as hd5 files.
        args.bam_file: BAM or BAMout file where the aligned reads are stored
        args.input_vcf: VCF file with annotation values from Haplotype caller or VQSR
        args.train_vcf: VCF file with true variant (from NIST or Platinum genomes, etc.)
        args.bed_file: High confidence intervals for the calls in args.train_vcf
        args.window_size: Size of sequence window around variant (width of the tensor)
        args.read_limit: Maximum number of reads to include (height of the tensor)
        args.chrom: Only write tensors from this chromosome (optional, used for parallelization)
        args.start_pos: Only write tensors after this position (optional, used for parallelization)
        args.end_pos: Only write tensors before this position (optional, used for parallelization)
    '''
    print('Writing tensors with:', args.tensor_name, 'channel map.')
    stats = Counter()

    samfile = pysam.AlignmentFile(args.bam_file, "rb")
    bed_dict = bed_file_to_dict(args.bed_file)
    record_dict = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, "fasta"))
    vcf_reader = vcf.Reader(open(args.input_vcf, 'r'))
    vcf_ram = vcf.Reader(open(args.train_vcf, 'rb'))

    if args.chrom:
        variants = vcf_reader.fetch(args.chrom, args.start_pos, args.end_pos)
    else:
        variants = vcf_reader

    for variant in variants:
        for allele_idx, allele in enumerate(variant.ALT):
            idx_offset, ref_start, ref_end = get_variant_window(args, variant)
            contig = record_dict[variant.CHROM]
            record = contig[ ref_start : ref_end ]

            cur_label_key = get_true_label(allele, variant, bed_dict, vcf_ram, stats)
            if not cur_label_key or downsample(args, cur_label_key, stats):
                continue

            if include_annotations:
                if all(map(
                        lambda x: x not in variant.INFO and x not in variant.FORMAT and x != "QUAL", args.annotations)):
                    stats['Missing ALL annotations'] += 1
                    continue # Require at least 1 annotation...
                annotation_data = get_annotation_data(args, variant, stats)

            good_reads, insert_dict = get_good_reads(args, samfile, variant)
            if len(good_reads) >= args.read_limit:
                stats['More reads than read_limit'] += 1
            if len(good_reads) == 0:
                stats['No reads aligned'] += 1
                continue

            reference_seq = record.seq
            for i in sorted(insert_dict.keys(), key=int, reverse=True):
                if i < 0:
                    reference_seq = vqsr_cnn.INDEL_CHAR*insert_dict[i] + reference_seq
                else:
                    reference_seq = reference_seq[:i] + vqsr_cnn.INDEL_CHAR*insert_dict[i] + reference_seq[i:]

            read_tensor = good_reads_to_tensor(args, good_reads, ref_start, insert_dict)
            reference_sequence_into_tensor(args, reference_seq, read_tensor)

            tensor_path = get_path_to_train_valid_or_test(args.data_dir)
            tensor_prefix = plain_name(args.input_vcf) +'_'+ plain_name(args.train_vcf)
            tensor_prefix += '_allele_' + str(allele_idx) + '-' + cur_label_key
            tensor_path += cur_label_key + '/' + tensor_prefix + '-' + variant.CHROM
            tensor_path += '_' + str(variant.POS) + vqsr_cnn.TENSOR_SUFFIX
            stats[cur_label_key] += 1

            if not os.path.exists(os.path.dirname(tensor_path)):
                os.makedirs(os.path.dirname(tensor_path))
            with h5py.File(tensor_path, 'w') as hf:
                if pileup:
                    pileup_tensor = read_tensor_to_pileup(args, read_tensor)
                    hf.create_dataset('pileup_tensor', data=pileup_tensor, compression='gzip')
                hf.create_dataset(args.tensor_name, data=read_tensor, compression='gzip')
                if include_annotations:
                    hf.create_dataset(args.annotation_set, data=annotation_data, compression='gzip')

            stats['count'] += 1
            if stats['count']%100 == 0:
                print('Wrote', stats['count'], 'tensors out of', args.samples, ' last variant:', str(variant))
        if stats['count'] >= args.samples:
            break

    for s in stats.keys():
        print(s, 'has:', stats[s])
    if variant:
        print('Done generating tensors. Last variant:', str(variant), 'from vcf:', args.input_vcf)


def train_on_reference_tensors_and_annotations(args):
    '''Train a 1D Convolution plus reference tracks and MLP Annotation architecture.

    Arguments:
        args.data_dir: must be set to an appropriate directory with
            subdirectories of test, valid and train, each containing
            subdirectories for each label with tensors stored as hd5 files.

    Reference and Annotation tensors must be generated by calling
    write_reference_and_annotation_tensors() before this function is used.
    Performance curves for CNN are plotted on the test dataset.
    '''
    train_paths, valid_paths, test_paths = get_train_valid_test_paths(args)

    generate_train = dna_annotation_generator(args, train_paths)
    generate_valid = dna_annotation_generator(args, valid_paths)

    weight_path = vqsr_cnn.weight_path_from_args(args)
    model = vqsr_cnn.build_reference_annotation_model(args)
    model = vqsr_cnn.train_model_from_generators(args, model, generate_train, generate_valid, weight_path)

    test = load_dna_annotations_positions_from_class_dirs(args, test_paths, per_class_max=args.samples)
    if args.image_dir:
        vqsr_cnn.plot_roc_per_class(model, [test[0], test[1]], test[2], args.labels, args.id, prefix=args.image_dir)



def train_on_read_tensors_and_annotations(args):
    '''Trains a reference, read, and annotation CNN architecture on tensors at the supplied data directory.

    This architecture looks at reads, read flags, reference sequence, and variant annotations.
    Tensors must be generated by calling write_read_and_annotation_tensors() before this function is used.
    After training with early stopping performance curves are plotted on the test dataset.

    Arguments:
        args.data_dir: must be set to an appropriate directory with
            subdirectories of test, valid and train, each containing
            subdirectories for each label with tensors stored as hd5 files.

    '''
    train_paths, valid_paths, test_paths = get_train_valid_test_paths(args)

    generate_train = tensor_generator_from_label_dirs_and_args(args, train_paths)
    generate_valid = tensor_generator_from_label_dirs_and_args(args, valid_paths)

    weight_path = vqsr_cnn.weight_path_from_args(args)
    model = vqsr_cnn.build_read_tensor_2d_and_annotations_model(args)
    model = vqsr_cnn.train_model_from_generators(args, model, generate_train, generate_valid, weight_path)

    test = load_tensors_and_annotations_from_class_dirs(args, test_paths, per_class_max=args.samples)
    if args.image_dir:
        vqsr_cnn.plot_roc_per_class(model, [test[0], test[1]], test[2], args.labels, args.id,
                                    prefix=args.image_dir, batch_size=args.batch_size)


def train_tiny_model_on_read_tensors_and_annotations(args):
    '''Trains a reference, read, and annotation CNN architecture on tensors at the supplied data directory.

    This architecture looks at reads, read flags, reference sequence, and variant annotations.
    Tensors must be generated by calling write_read_and_annotation_tensors() before this function is used.
    After training with early stopping performance curves are plotted on the test dataset.

    Arguments:
        args.data_dir: must be set to an appropriate directory with
            subdirectories of test, valid and train, each containing
            subdirectories for each label with tensors stored as hd5 files.

    '''
    train_paths, valid_paths, test_paths = get_train_valid_test_paths(args)

    generate_train = tensor_generator_from_label_dirs_and_args(args, train_paths)
    generate_valid = tensor_generator_from_label_dirs_and_args(args, valid_paths)

    weight_path = vqsr_cnn.weight_path_from_args(args)
    model = vqsr_cnn.build_tiny_2d_annotation_model(args)
    model = vqsr_cnn.train_model_from_generators(args, model, generate_train, generate_valid, weight_path)

    test = load_tensors_and_annotations_from_class_dirs(args, test_paths, per_class_max=args.samples)
    if args.image_dir:
        vqsr_cnn.plot_roc_per_class(model, [test[0], test[1]], test[2], args.labels, args.id,
                                    prefix=args.image_dir, batch_size=args.batch_size)


def train_small_model_on_read_tensors_and_annotations(args):
    '''Trains a reference, read, and annotation CNN architecture on tensors at the supplied data directory.

    This architecture looks at reads, read flags, reference sequence, and variant annotations.
    Tensors must be generated by calling write_read_and_annotation_tensors() before this function is used.
    After training with early stopping performance curves are plotted on the test dataset.

    Arguments:
        args.data_dir: must be set to an appropriate directory with
            subdirectories of test, valid and train, each containing
            subdirectories for each label with tensors stored as hd5 files.

    '''
    train_paths, valid_paths, test_paths = get_train_valid_test_paths(args)

    generate_train = tensor_generator_from_label_dirs_and_args(args, train_paths)
    generate_valid = tensor_generator_from_label_dirs_and_args(args, valid_paths)

    weight_path = vqsr_cnn.weight_path_from_args(args)
    model = vqsr_cnn.build_small_2d_annotation_model(args)
    model = vqsr_cnn.train_model_from_generators(args, model, generate_train, generate_valid, weight_path)

    test = load_tensors_and_annotations_from_class_dirs(args, test_paths, per_class_max=args.samples)
    if args.image_dir:
        vqsr_cnn.plot_roc_per_class(model, [test[0], test[1]], test[2], args.labels, args.id,
                                    prefix=args.image_dir, batch_size=args.batch_size)


def get_annotation_data(args, annotation_variant, stats):
    '''Return an array annotation data about the variant.

    Arguments:
        args.annotations: List of variant annotations to use
        annotation_variant: the variant with annotation
        stats: Counter of run statistics

    Returns:
        annotation_data: numpy array of annotation values
    '''
    annotation_data = np.zeros((len(args.annotations),))

    for i, a in enumerate(args.annotations):
        if a == 'QUAL':
            annotation_data[i] = annotation_variant.QUAL
        elif a == 'AF':
            annotation_data[i] = annotation_variant.INFO[a][0]
        elif a in annotation_variant.INFO and not math.isnan(annotation_variant.INFO[a]):
            annotation_data[i] = annotation_variant.INFO[a]
        elif a == 'MBQ':
            call = annotation_variant.genotype(args.sample_name)
            annotation_data[i] = call.data.MBQ
        elif a == 'MPOS':
            call = annotation_variant.genotype(args.sample_name)
            annotation_data[i] = call.data.MPOS
        elif a == 'MMQ':
            call = annotation_variant.genotype(args.sample_name)
            annotation_data[i] = call.data.MMQ
        elif a == 'MFRL_0':
            call = annotation_variant.genotype(args.sample_name)
            annotation_data[i] = call.data.MFRL[0]
        elif a == 'MFRL_1':
            call = annotation_variant.genotype(args.sample_name)
            annotation_data[i] = call.data.MFRL[1]
        elif a == 'AD_0':
            call = annotation_variant.genotype(args.sample_name)
            annotation_data[i] = call.data.AD[0]
        elif a == 'AD_1':
            call = annotation_variant.genotype(args.sample_name)
            annotation_data[i] = call.data.AD[1]
        else:
            stats['Could not find annotation:' + a] += 1

    return annotation_data


def get_good_reads(args, samfile, variant, sort_by='base'):
    '''Return an array of usable reads centered at the variant.

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
        good_reads: array of usable reads sorted by reference start position
        insert_dict: a dict mapping read indices to max insertions at that point
    '''
    good_reads = []
    insert_dict = {}

    idx_offset, ref_start, ref_end = get_variant_window(args, variant)

    for read in samfile.fetch(variant.CHROM, variant.POS-1, variant.POS+1):

        if not read or not hasattr(read, 'cigarstring') or read.cigarstring is None:
            continue

        read_group = read.get_tag('RG')
        if 'artificial' in read_group.lower():
            continue

        index_dif = ref_start - read.reference_start
        if abs(index_dif) >= args.window_size:
            continue

        if 'I' in read.cigarstring:
            cur_idx = 0
            for t in read.cigartuples:
                if t[0] == vqsr_cnn.CIGAR_CODE['I']:
                    insert_idx = cur_idx - index_dif
                    if insert_idx not in insert_dict:
                        insert_dict[insert_idx] = t[1]
                    elif insert_dict[insert_idx] < t[1]:
                        insert_dict[insert_idx] = t[1]

                if t[0] in [vqsr_cnn.CIGAR_CODE['M'], vqsr_cnn.CIGAR_CODE['I'],
                            vqsr_cnn.CIGAR_CODE['S'], vqsr_cnn.CIGAR_CODE['D']]:
                    cur_idx += t[1]

        good_reads.append(read)

    if len(good_reads) > args.read_limit:
        good_reads = np.random.choice(good_reads, size=args.read_limit, replace=False).tolist()

    good_reads.sort(key=lambda x: x.reference_start + x.query_alignment_start)
    if sort_by == 'base':
        good_reads.sort(key=lambda read: get_base_to_sort_by(read, variant))

    return good_reads, insert_dict


def get_base_to_sort_by(read, variant):
    if len(read.query_alignment_sequence) > 0:
        max_idx = len(read.query_alignment_sequence)-1
    else:
        return 'Z'

    if variant.is_snp:
        return read.query_alignment_sequence[clamp((variant.POS-read.reference_start)-1, 0, max_idx)]
    elif variant.is_indel:
        var_idx = variant.POS-read.reference_start
        cur_idx = 0
        for cur_op, length in read.cigartuples:
            cur_idx += length
            if cur_idx > var_idx:
                if cur_op == vqsr_cnn.CIGAR_CODE['M']:
                    return read.query_alignment_sequence[clamp(var_idx, 0, max_idx)]
                else:
                    return vqsr_cnn.CODE2CIGAR[cur_op]
        return 'Y'


def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)


def good_reads_to_tensor(args, good_reads, ref_start, insert_dict):
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
    channel_map = vqsr_cnn.get_tensor_channel_map_from_args(args)
    tensor = np.zeros( vqsr_cnn.tensor_shape_from_args(args) )

    for j,read in enumerate(good_reads):

        rseq, rqual = sequence_and_qualities_from_read(args, read, ref_start, insert_dict)
        flag_start = -1
        flag_end = 0

        for i,b in enumerate(rseq):

            if i == args.window_size:
                break

            if b == vqsr_cnn.SKIP_CHAR:
                continue
            elif flag_start == -1:
                flag_start = i
            else:
                flag_end = i

            if b in args.input_symbols:
                if b == vqsr_cnn.INDEL_CHAR:
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

            elif b in vqsr_cnn.AMBIGUITY_CODES:
                if K.image_data_format() == 'channels_last':
                    tensor[j, i, :4] = vqsr_cnn.AMBIGUITY_CODES[b]
                else:
                    tensor[:4, j, i] = vqsr_cnn.AMBIGUITY_CODES[b]

            else:
                print('Error! Unknown symbol in seq block:', b)
                return

        flags = flag_to_array(read.flag)
        for i in range(vqsr_cnn.READ_FLAGS):
            flag_str = 'flag_bit_'+ str(i)

            if flags[i] and flag_str in channel_map:
                if K.image_data_format() == 'channels_last':
                    tensor[j, flag_start:flag_end, channel_map[flag_str]] = 1.0
                else:
                    tensor[channel_map[flag_str], j,  flag_start:flag_end] = 1.0

        if 'mapping_quality' in channel_map:
            mq = float(read.mapping_quality)/vqsr_cnn.MAPPING_QUALITY_MAX
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
    for t in read.cigartuples:
        my_ref_idx = cur_idx - index_dif
        if t[0] == vqsr_cnn.CIGAR_CODE['I'] and my_ref_idx in insert_dict:
            my_indel_dict[my_ref_idx] = insert_dict[my_ref_idx] - t[1]
        elif t[0] == vqsr_cnn.CIGAR_CODE['D']:
            my_indel_dict[my_ref_idx] = t[1]
        if t[0] in [vqsr_cnn.CIGAR_CODE['M'], vqsr_cnn.CIGAR_CODE['I'],
                    vqsr_cnn.CIGAR_CODE['S'], vqsr_cnn.CIGAR_CODE['D']]:
            cur_idx += t[1]

    for k in insert_dict.keys():
        if k not in my_indel_dict:
            my_indel_dict[k] = insert_dict[k]

    rseq = read.query_alignment_sequence[:args.window_size]
    rqual = read.query_alignment_qualities[:args.window_size].tolist()

    if index_dif > 0:
        rseq = rseq[index_dif:]
        rqual = rqual[index_dif:]
    elif index_dif < 0:
        rseq = vqsr_cnn.SKIP_CHAR*(-index_dif) + rseq
        rqual = [no_qual_filler]*(-index_dif) + rqual

    for j in sorted(my_indel_dict.keys(), key=int, reverse=True):
        if j < 1:
            rseq = (vqsr_cnn.INDEL_CHAR*my_indel_dict[j]) + rseq
            rqual = ([no_qual_filler]*my_indel_dict[j]) + rqual
        else:
            rseq = rseq[:j] + (vqsr_cnn.INDEL_CHAR*my_indel_dict[j]) + rseq[j:]
            rqual = rqual[:j] + ([no_qual_filler]*my_indel_dict[j]) + rqual[j:]

    return rseq, rqual


def read_tensor_to_pileup(args, read_tensor):
    tensor_map = vqsr_cnn.get_tensor_channel_map_from_args(args)
    channels = vqsr_cnn.get_reference_and_read_channels(args)
    pileup_tensor = np.zeros((args.window_size, channels))

    for i in range(args.window_size):
        for key in tensor_map:
            if 'read' not in key and 'reference' not in key:
                continue

            if 'read' in key and K.image_data_format() == 'channels_last':
                pileup_tensor[i, tensor_map[key]] = np.sum(read_tensor[:, i, tensor_map[key]]) / args.window_size
            elif 'read' in key:
                pileup_tensor[i, tensor_map[key]] = np.sum(read_tensor[tensor_map[key], :, i]) / args.window_size
            elif 'reference' in key and K.image_data_format() == 'channels_last':
                pileup_tensor[i, tensor_map[key]] = np.amax(read_tensor[:, i, tensor_map[key]])
            elif 'reference' in key:
                pileup_tensor[i, tensor_map[key]] = np.amax(read_tensor[tensor_map[key], :, i])
            else:
                raise ValueError('Error unexpected key:'+key)

    return pileup_tensor


def reference_sequence_into_tensor(args, reference_seq, tensor):
    ref_offset = len(set(args.input_symbols.values()))
    for i,b in enumerate(reference_seq):
        if i == args.window_size:
            break
        if b in args.input_symbols:
            if K.image_data_format() == 'channels_last':
                tensor[:, i, ref_offset+args.input_symbols[b]] = 1.0
            else:
                tensor[ref_offset+args.input_symbols[b], :, i] = 1.0
        elif b in vqsr_cnn.AMBIGUITY_CODES:
            ambiguous_vector = np.tile(vqsr_cnn.AMBIGUITY_CODES[b], (args.read_limit, 1))
            if K.image_data_format() == 'channels_last':
                tensor[:, i, ref_offset:ref_offset+4] = ambiguous_vector
            else:
                tensor[ref_offset:ref_offset+4, :, i] = np.transpose(ambiguous_vector)


def flag_to_array(flag):
    flags = []

    for i in range(vqsr_cnn.READ_FLAGS):
        flags.append((flag>>i)&1)

    return np.array(flags)


def add_flags_to_read_tensor(args, tensor, tensor_channel_map, flags):
    for k in tensor_channel_map.keys():
        if 'flag' in k:
            flag_bit = int(k.split('_')[-1])
            for read_idx in range(flags.shape[1]):
                if K.image_data_format() == 'channels_last':
                    tensor[read_idx, :, tensor_channel_map[k]] = flags[flag_bit, read_idx]
                else:
                    tensor[tensor_channel_map[k], read_idx, :] = flags[flag_bit, read_idx]


def add_mq_to_read_tensor(args, tensor, tensor_channel_map, mapping_qualities):
    if not 'mapping_quality' in tensor_channel_map:
        return

    for read_idx, mq in enumerate(mapping_qualities):
        if K.image_data_format() == 'channels_last':
            tensor[read_idx, :, tensor_channel_map['mapping_quality']] = float(mq) / vqsr_cnn.MAPPING_QUALITY_MAX
        else:
            tensor[tensor_channel_map['mapping_quality'], read_idx, :] = float(mq) / vqsr_cnn.MAPPING_QUALITY_MAX


def base_quality_to_phred_array(base_quality, base, base_dict):
    phred = np.zeros((4,))
    exponent = float(-base_quality) / 10.0
    p = 1.0-(10.0**exponent) # Convert to probability
    not_p = (1.0-p) / 3.0 # Error could be any of the other 3 bases
    not_base_quality = -10 * np.log10(not_p) # Back to Phred

    for b in base_dict.keys():
        if b == vqsr_cnn.INDEL_CHAR:
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
        elif b == vqsr_cnn.INDEL_CHAR:
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
        raise ValueError('Error! Unknown base quality mode:', args.base_quality_mode)



def get_true_label(allele, variant, bed_dict, truth_vcf, stats):
    '''Defines the truth status of a variant allele given a truth vcf and confident region.

    Arguments:
        allele: The allele to check
        variant: the variant whose allele we will check
        bed_dict: confident region dict defined by intervals e.g. from bed_file_to_dict()
        truth_vcf: vcf of validated variants
        stats: Counter dict used to keep track of the label distribution, etc.

    Returns:
        None if outside the confident region
        Otherwise a label string:
            SNP if variant is snp and in truth vcf
            INDEL if variant is indel and in truth vcf
            NOT_SNP if variant is snp and not in truth vcf
            NOT_INDEL if variant is indel and not in truth vcf
    '''
    in_bed = in_bed_file(bed_dict, variant.CHROM, variant.POS)
    if allele_in_vcf(allele, variant, truth_vcf) and in_bed:
        class_prefix = ''
    elif in_bed:
        class_prefix = 'NOT_'
    else:
        stats['Variant outside confident bed file'] += 1
        return None

    if variant.is_snp:
        cur_label_key = class_prefix + 'SNP'
    elif variant.is_indel:
        cur_label_key = class_prefix + 'INDEL'
    else:
        stats['Not SNP or INDEL'] += 1
        return None

    return cur_label_key


def downsample(args, cur_label_key, stats):
    '''Indicates whether or not to downsample a variant.

    Arguments:
        args.skip_positive_class: Skip all positive examples
        args.downsample_snps: fraction of SNPs to keep
        args.downsample_indels: fraction of INDELs to keep
        cur_label_key: truth label from get_true_label()
        stats: Counter dict used to keep track of a run

    Returns:
        Boolean: should we downsample this variant or not.
    '''
    if args.skip_positive_class and cur_label_key in ['SNP', 'INDEL']:
        stats['Downsampled positive examples'] += 1
        return True

    if args.downsample_snps < 1.0 and cur_label_key == 'SNP':
        dice = np.random.rand()
        if dice > args.downsample_snps:
            stats['Downsampled SNPs'] += 1
            return True
    elif args.downsample_indels < 1.0 and cur_label_key == 'INDEL':
        dice = np.random.rand()
        if dice > args.downsample_indels:
            stats['Downsampled INDELs'] += 1
            return True
    if args.downsample_not_snps < 1.0 and cur_label_key == 'NOT_SNP':
        dice = np.random.rand()
        if dice > args.downsample_not_snps:
            stats['Downsampled NOT_SNPs'] += 1
            return True
    elif args.downsample_not_indels < 1.0 and cur_label_key == 'NOT_INDEL':
        dice = np.random.rand()
        if dice > args.downsample_not_indels:
            stats['Downsampled NOT_INDELs'] += 1
            return True

    return False


def interval_file_to_dict(interval_file, shift1=0, skip=['@']):
    ''' Create a dict to store intervals from a interval list file.

    Arguments:
        interval_file: the file to load either a bed file -> shift1 should be 1
            or a picard style interval_list file -> shift1 should be 0
        shift1: Shift the intervals 1 position over to align with 1-indexed VCFs
        skip: Comment character to ignore
    Returns:
        intervals: dict where keys in the dict are contig ids
            values are a tuple of arrays the first array
            in the tuple contains the start positions
            the second array contains the end positions.
    '''
    intervals = {}

    with open(interval_file) as f:
        for line in f:
            if line[0] in skip:
                continue

            parts = line.split()
            contig = parts[0]
            lower = int(parts[1])+shift1
            upper = int(parts[2])+shift1

            if contig not in intervals:
                intervals[contig] = ([], [])

            intervals[contig][0].append(lower)
            intervals[contig][1].append(upper)

    for k in intervals.keys():
        intervals[k] = (np.array(intervals[k][0]), np.array(intervals[k][1]))

    return intervals


def bed_file_to_dict(bed_file):
    return interval_file_to_dict(bed_file, shift1=1)


def in_bed_file(bed_dict, contig, pos):
    # Exclusive
    lows = bed_dict[contig][0]
    ups = bed_dict[contig][1]
    return np.any((lows <= pos) & (pos < ups))


def allele_in_vcf(allele, variant, vcf_ram):
    ''' Check if variant's allele is in a VCF file.

    Arguments
        allele: the allele from the provided variant that we are checking
        variant: the variant whose allele we are looking for
        vcf_ram: the VCF we look in, must have an index (tbi, or idx)

    Returns
        variant if it is found otherwise None
    '''
    variants = vcf_ram.fetch(variant.CHROM, variant.POS-1, variant.POS)

    for v in variants:
        if v.CHROM == variant.CHROM and v.POS == variant.POS and allele in v.ALT:
            return v

    return None


def get_variant_window(args, variant):
    index_offset = (args.window_size//2)
    reference_start = variant.POS-(index_offset+1)
    reference_end = variant.POS+index_offset

    return index_offset, reference_start, reference_end


def dna_annotation_generator(args, train_paths):
    """Data generator of DNA and annotation tensors.

    Assumes train paths contains example in labelled directories.
    Loops over all examples sampling args.batch_size examples
    uniformly from each label.

    Arguments:
        args: args object needed for batch_size, labels, and annotations
        train_paths: array of label directories with hd5 tensors within each

    Returns:
        A tuple with a dict of the input tensors
        and a 1-Hot matrix (2D numpy array) of the labels.
    """
    per_batch_per_label = (args.batch_size // len(args.labels))
    tensor_counts = Counter()
    tensors = {}

    if args.window_size > 0:
        channel_map = vqsr_cnn.get_tensor_channel_map_from_args(args)
        tensor = np.zeros((args.batch_size, args.window_size, len(channel_map)))

    annotation_data = np.zeros((args.batch_size, len(args.annotations)))
    label_matrix = np.zeros((args.batch_size, len(args.labels)))


    for tp in train_paths:
        label_key = os.path.basename(tp)
        if label_key not in args.labels:
            print('Skipping label directory:', label_key, ' which is not in args label set:', args.labels.keys())
            continue
        label = args.labels[label_key]

        tensors[label] = [os.path.join(tp, t) for t in os.listdir(tp)
                          if os.path.splitext(t)[1] == vqsr_cnn.TENSOR_SUFFIX]
        tensor_counts[label] = 0
        print('Found ', len(tensors[label]), 'examples of label:', label, 'in:', tp)

    while True:
        cur_example = 0
        for label in tensors.keys():
            for i in range(per_batch_per_label):
                tensor_path = tensors[label][tensor_counts[label]]
                label_matrix[cur_example, label] = 1.0
                with h5py.File(tensor_path,'r') as hf:
                    annotation_data[cur_example,:] = np.array(hf.get(args.annotation_set))
                    if args.window_size > 0:
                        tensor[cur_example,:,:] = np.array(hf.get(args.tensor_name))

                tensor_counts[label] += 1
                if tensor_counts[label] == len(tensors[label]):
                    np.random.shuffle(tensors[label])
                    print('\nGenerator shuffled & looped over:', tensor_counts[label],
                          'examples of label:',label, '\nLast tensor was:', tensor_path)
                    tensor_counts[label] = 0
                cur_example += 1
                if cur_example == args.batch_size:
                    break

        if args.window_size > 0:
            yield ({args.tensor_name:tensor, args.annotation_set:annotation_data}, label_matrix)
        else:
            yield (annotation_data, label_matrix)



def tensor_generator_from_label_dirs_and_args(args, train_paths, with_positions=False):
    """Data generator of tensors with reads, and annotations.

    Assumes train paths contains example in labelled directories.
    Loops over all examples sampling args.batch_size examples
    uniformly from each label.

    Arguments:
        args: args object needed for batch_size, labels, and annotations
        train_paths: array of label directories with hd5 tensors within each
        with_positions: boolean if True will include a position string
            (i.e. "1_1234_0" for tensor from contig one base 1234 and first allele)
            as the last element in each tensor tuple.
    Returns:
        A tuple with a dict of the input tensors
        and a 1-Hot matrix (2D numpy array) of the labels.
    """
    batch = {}
    tensors = {}
    tensor_counts = Counter()
    per_batch_per_label = (args.batch_size // len(args.labels) )

    tm = vqsr_cnn.get_tensor_channel_map_from_args(args)
    if tm:
        tensor_shape = vqsr_cnn.tensor_shape_from_args(args)
        batch[args.tensor_name] = np.zeros(((args.batch_size,)+tensor_shape))

    if vqsr_cnn.annotations_from_args(args):
        batch[args.annotation_set] = np.zeros((args.batch_size, len(args.annotations)))

    if with_positions:
        positions = []

    label_matrix = np.zeros((args.batch_size, len(args.labels)))

    for tp in train_paths:
        label_key = os.path.basename(tp)
        if label_key not in args.labels:
            print('Skipping label directory:', label_key, ' which is not in args label set:', args.labels.keys())
            continue
        label = args.labels[label_key]
        tensors[label] = [os.path.join(tp, t) for t in os.listdir(tp)
                          if os.path.splitext(t)[1] == vqsr_cnn.TENSOR_SUFFIX]
        tensor_counts[label] = 0
        print('Found ', len(tensors[label]), 'examples of label:', label, 'in:', tp)

    while True:
        cur_example = 0
        for label in tensors.keys():
            for i in range(per_batch_per_label):
                tensor_path = tensors[label][tensor_counts[label]]

                with h5py.File(tensor_path, 'r') as hf:
                    for key in batch.keys():
                        batch[key][cur_example] = np.array(hf.get(key))

                label_matrix[cur_example, label] = 1.0
                tensor_counts[label] += 1
                if tensor_counts[label] == len(tensors[label]):
                    np.random.shuffle(tensors[label])
                    print('\nGenerator looped over:', tensor_counts[label],
                          'examples of label:', label, '\nShuffled them. Last tensor was:', tensor_path)
                    tensor_counts[label] = 0

                if with_positions:
                    positions.append(position_string_from_tensor_name(tensor_path))

                cur_example += 1
                if cur_example == args.batch_size:
                    break

        if with_positions:
            yield (batch, label_matrix, positions)
            positions = []
        else:
            yield (batch, label_matrix)
        label_matrix = np.zeros((args.batch_size, len(args.labels)))
        if with_positions and tm:
            tensor_shape = vqsr_cnn.tensor_shape_from_args(args)
            batch[args.tensor_name] = np.zeros(((args.batch_size,)+tensor_shape))

        if with_positions and vqsr_cnn.annotations_from_args(args):
            batch[args.annotation_set] = np.zeros((args.batch_size, len(args.annotations)))


def load_dna_annotations_positions_from_class_dirs(args, train_paths,
                                                   per_class_max=4000, include_dna=True, include_annotations=True):
    count = 0

    annotation_data = []
    reference_data = []
    labels_data = []
    positions = []

    for tp in train_paths:
        label_key = os.path.basename(tp)
        if label_key not in args.labels:
            print('Skipping label directory:', label_key, ' which is not in args label set:', args.labels.keys())
            continue
        label = args.labels[label_key]
        imgs = os.listdir(tp)
        count += 1
        print(count, " dir out of:", len(train_paths), tp, "has:", len(imgs))
        this_t = 0
        for t in imgs:
            this_t += 1
            if this_t > per_class_max:
                print('Per class max reached. bailing at', this_t)
                break

            fn, file_extension = os.path.splitext(t)
            if not file_extension.lower() == vqsr_cnn.TENSOR_SUFFIX:
                continue

            with h5py.File(tp+'/'+t, 'r') as hf:
                if include_annotations:
                    annotation_data.append(np.array(hf.get(args.annotation_set)))
                if include_dna:
                    reference_data.append(np.array(hf.get(args.tensor_name)))

            y_vector = np.zeros(len(args.labels)) # One hot Y vector of size labels, correct label is 1 others are 0
            y_vector[label] = 1.0
            labels_data.append(y_vector)
            positions.append(position_string_from_tensor_name(t))

    if include_dna and include_annotations:
        return np.asarray(reference_data), np.asarray(annotation_data), np.asarray(labels_data), np.asarray(positions)
    elif include_annotations:
        return np.asarray(annotation_data), np.asarray(labels_data), np.asarray(positions)
    elif include_dna:
        return np.asarray(reference_data), np.asarray(labels_data), np.asarray(positions)


def load_tensors_and_annotations_from_class_dirs(args, train_paths, per_class_max=2500, position_dict=None):
    annotations = []
    positions = []
    tensors = []
    labels = []
    count = 0

    for tp in train_paths:
        label_key = os.path.basename(tp)
        if label_key not in args.labels:
            print('Skipping label directory:', label_key, ' which is not in args label set:', args.labels.keys())
            continue

        label = args.labels[label_key]
        imgs = os.listdir(tp)
        count += 1
        this_t = 0
        for t in imgs:
            if this_t > per_class_max:
                print('Per class max reached. bailing at', this_t)
                break

            fn, file_extension = os.path.splitext(t)
            if not file_extension.lower() == vqsr_cnn.TENSOR_SUFFIX:
                continue

            with h5py.File(tp+'/'+t, 'r') as hf:
                tensors.append(np.array(hf.get(args.tensor_name)))
                annotations.append(np.array(hf.get(args.annotation_set)))

            y_vector = np.zeros(len(args.labels)) # One hot Y vector of size labels, correct label is 1 all others are 0
            y_vector[label] = 1.0
            labels.append(y_vector)
            positions.append(position_string_from_tensor_name(t))
            this_t += 1

        print(count, " dir out of:", len(train_paths), tp, "has:", len(imgs), 'Loaded:', this_t)

    return np.asarray(tensors), np.asarray(annotations), np.asarray(labels), np.asarray(positions)


def position_string_from_tensor_name(tensor_name):
    '''Genomic position as underscore delineated string from a filename.

    Includes an allele index if the filename includes _allele_
    This is ugly, we need file names ending with genomic position
    (e.g. my_tensor-12_1234.h5 returns 12_1234 and a_tensor_allele_1-8_128.hd5 returns 8_128_1)

    Arguments:
        tensor_name: the filename to parse
    Returns:
        Genomic position string Contig_Position or Contig_Position_AlleleIndex
    '''
    slash_split = tensor_name.split('/')
    dash_split = slash_split[-1].split('-')
    gsplit = dash_split[0].split('_')

    gpos = dash_split[-1]
    chrom = gpos.split('_')[0]
    pos = os.path.splitext(gpos.split('_')[1])[0]
    pos_str = chrom + '_' + pos

    for i,p in enumerate(gsplit):
        if p == 'allele':
            pos_str += '_'+str(gsplit[i+1])

    return pos_str


def get_path_to_train_valid_or_test(path, valid_ratio=0.1, test_ratio=0.2, valid_contig='-19_', test_contig='-20_'):
    dice = np.random.rand()
    if dice < valid_ratio or valid_contig in path:
        return os.path.join(path, 'valid/')
    elif dice < valid_ratio+test_ratio or test_contig in path:
        return os.path.join(path, 'test/')
    else:
        return os.path.join(path, 'train/')


def get_train_valid_test_paths(args):
    train_dir = args.data_dir + 'train/'
    valid_dir = args.data_dir + 'valid/'
    test_dir = args.data_dir + 'test/'
    train_paths = [train_dir + tp for tp in sorted(os.listdir(train_dir)) if os.path.isdir(train_dir + tp)]
    valid_paths = [valid_dir + vp for vp in sorted(os.listdir(valid_dir)) if os.path.isdir(valid_dir + vp)]
    test_paths = [test_dir + vp for vp in sorted(os.listdir(test_dir)) if os.path.isdir(test_dir + vp)]

    assert(len(train_paths) == len(valid_paths) == len(test_paths))

    return train_paths, valid_paths, test_paths


def plain_name(full_name):
    name = os.path.basename(full_name)
    return name.split('.')[0]


# Back to the top!
if "__main__" == __name__:
    run()