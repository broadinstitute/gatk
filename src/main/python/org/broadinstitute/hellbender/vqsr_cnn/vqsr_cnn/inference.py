# Imports
import os
import re
import sys
import h5py
import time
import json
import argparse
import numpy as np
from collections import Counter, defaultdict, namedtuple

# Keras Imports
import keras.backend as K
from keras.models import load_model


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Definitions ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dna_symbols = {'A':0, 'C':1, 'G':2, 'T':3}
inputs_indel = {'A':0, 'C':1, 'G':2, 'T':3, '*':4}

# Base calling ambiguities, See https://www.bioinformatics.org/sms/iupac.html
ambiguity_codes = {'K':[0,0,0.5,0.5], 'M':[0.5,0.5,0,0], 'R':[0.5,0,0,0.5], 'Y':[0,0.5,0.5,0], 'S':[0,0.5,0,0.5], 'W':[0.5,0,0.5,0], 
				  'B':[0,0.333,0.333,0.334], 'V':[0.333,0.333,0,0.334],'H':[0.333,0.333,0.334,0],'D':[0.333,0,0.333,0.334],
				  'X':[0.25,0.25,0.25,0.25], 'N':[0.25,0.25,0.25,0.25]}


tensor_maps_2d = ['read_tensor']
tensor_maps_1d = ['reference']

# Annotation sets
annotations = {
	'_' : [], # Allow command line to unset annotations
	'gatk_w_qual' : ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'QUAL', 'ReadPosRankSum'],
	'gatk' : ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum'],
	'annotations' : ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum'],
	'm2':['AF', 'AD_0', 'AD_1', 'MBQ', 'MFRL_0', 'MFRL_1', 'MMQ', 'MPOS'],
	'combine': ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum', 'AF', 'AD_0', 'AD_1', 'MBQ', 'MFRL_0', 'MFRL_1', 'MMQ', 'MPOS'],
	'gnomad': ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum', 'DP_MEDIAN', 'DREF_MEDIAN', 'GQ_MEDIAN', 'AB_MEDIAN'],
}

snp_indel_labels = {'NOT_SNP':0, 'NOT_INDEL':1, 'SNP':2, 'INDEL':3}

CODE2CIGAR = 'MIDNSHP=XB'
CIGAR2CODE = dict([y, x] for x, y in enumerate(CODE2CIGAR))
cigar_code = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4}
CIGAR_REGEX = re.compile("(\d+)([MIDNSHP=XB])")

eps = 1e-7
skip_char = '~'
indel_char = '*'
separator_char = '\t'
mapping_quality_max = 60.0

read_elements = 8
Read = namedtuple("Read", "seq qual cigar reverse mate_reverse first mapping_quality reference_start")
Variant = namedtuple("Variant", "contig pos ref alt type")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Inference ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def score_and_write_batch(args, model, file_out, fifo, batch_size, python_batch_size):
	'''Score and write a batch of variants with a 1D CNN.

	This function is tightly coupled with the NeuralNetInference
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
		fifo_data = fifo_line.split(separator_char)

		variant_data.append(fifo_data[0] + '\t' + fifo_data[1] + '\t' + fifo_data[2] + '\t' + fifo_data[3])
		reference_batch.append(reference_string_to_tensor(fifo_data[4]))
		annotation_batch.append(annotation_string_to_tensor(fifo_data[5]))
		variant_types.append(fifo_data[6])
		fidx = 7
		if args.tensor_map in tensor_maps_2d and len(fifo_data) > fidx:
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
				fidx += read_elements
				print('Got read, seq data fidx:', read_tuples[-1])
			_, ref_start, _ = get_variant_window(args, var)
			insert_dict = get_inserts(args, read_tuples, var)
			read_batch.append(read_tuples_to_read_tensor(args, read_tuples, ref_start, insert_dict))
			print('read_batch len:', len(read_batch))

	if args.tensor_map in tensor_maps_1d:
		predictions = model.predict([np.array(reference_batch), np.array(annotation_batch)], batch_size=python_batch_size)
	elif args.tensor_map in tensor_maps_2d:
		if len(read_batch) > 0:
			predictions = model.predict({'read_tensor':np.array(read_batch), 'annotations':np.array(annotation_batch)}, batch_size=python_batch_size)
	else:
		raise ValueError('Unknown tensor mapping.  Check architecture file.', args.tensor_map)

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
	dna_data = np.zeros( (len(reference), len(dna_symbols)) )
	for i,b in enumerate(reference):
		if b in dna_symbols:
			dna_data[i, dna_symbols[b]] = 1.0
		elif b in ambiguity_codes:
			dna_data[i] = ambiguity_codes[b]
		else:
			raise ValueError('Error! Unknown code:', b)

	return dna_data


def annotation_string_to_tensor(annotation_string):
	name_val_pairs = annotation_string.split(';')
	name_val_arrays = [p.split('=') for p in name_val_pairs]
	annotation_map = {str(p[0]).strip() : p[1] for p in name_val_arrays if len(p) > 1}
	annotation_data = np.zeros(( len(annotations['gatk']), ))
	
	for i,a in enumerate(annotations['gatk']):
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
				if t[0] == cigar_code['I']:
					insert_idx = cur_idx - index_dif
					if insert_idx not in insert_dict:
						insert_dict[insert_idx] = t[1]
					elif insert_dict[insert_idx] < t[1]:
						insert_dict[insert_idx] = t[1]

				if t[0] in [cigar_code['M'], cigar_code['I'], cigar_code['S'], cigar_code['D']]:
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
				if cur_op == cigar_code['M']:
					return read.seq[clamp(var_idx, 0, max_idx)]
				else:
					return CODE2CIGAR[cur_op]
		return 'Y'


def cigar_string_to_tuples(cigar):
	if not cigar or len(cigar) == 0:
		return []
	parts = CIGAR_REGEX.findall(cigar)
	# reverse order
	return [(CIGAR2CODE[y], int(x)) for x,y in parts]


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
	channel_map = get_tensor_channel_map_from_args(args)
	tensor = np.zeros( tensor_shape_from_args(args) )

	for j,read in enumerate(read_tuples):
		rseq, rqual = sequence_and_qualities_from_read(args, read, ref_start, insert_dict)
		flag_start = -1
		flag_end = 0

		for i,b in enumerate(rseq):

			if i == args.window_size:
				break
			if b == skip_char:
				continue
			elif flag_start == -1:
				flag_start = i
			else:
				flag_end = i

			if b in args.input_symbols:
				if b == indel_char:
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
			elif b in ambiguity_codes:
				if K.image_data_format() == 'channels_last':
					tensor[j, i, :4] = ambiguity_codes[b]
				else:
					tensor[:4, j, i] = ambiguity_codes[b]
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
			if K.image_data_format() == 'channels_last':
				tensor[j, flag_start:flag_end, channel_map['mapping_quality']] = float(read.mapping_quality)/mapping_quality_max
			else:
				tensor[channel_map['mapping_quality'], j, flag_start:flag_end] = float(read.mapping_quality)/mapping_quality_max
	return tensor


def sequence_and_qualities_from_read(args, read, ref_start, insert_dict):
	cur_idx = 0
	my_indel_dict = {}
	no_qual_filler = 0

	index_dif = ref_start - read.reference_start
	for t in cigar_string_to_tuples(read.cigar):
		my_ref_idx = cur_idx - index_dif
		if t[0] == cigar_code['I'] and my_ref_idx in insert_dict:
			my_indel_dict[my_ref_idx] = insert_dict[my_ref_idx] - t[1]
		elif t[0] == cigar_code['D']:
			my_indel_dict[my_ref_idx] = t[1]
		if t[0] in [cigar_code['M'], cigar_code['I'], cigar_code['S'], cigar_code['D']]:
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
		rseq = skip_char*(-index_dif) + rseq
		rqual = [no_qual_filler]*(-index_dif) + rqual

	for j in sorted(my_indel_dict.keys(), key=int, reverse=True):
		if j < 1:
			rseq = (indel_char*my_indel_dict[j]) + rseq
			rqual = ([no_qual_filler]*my_indel_dict[j]) + rqual
		else:
			rseq = rseq[:j] + (indel_char*my_indel_dict[j]) + rseq[j:]
			rqual = rqual[:j] + ([no_qual_filler]*my_indel_dict[j]) + rqual[j:]

	return rseq, rqual


def base_quality_to_phred_array(base_quality, base, base_dict):
	phred = np.zeros((4,))
	exponent = float(-base_quality) / 10.0
	p = 1.0-(10.0**exponent) # Convert to probability
	not_p = (1.0-p) / 3.0 # Error could be any of the other 3 bases
	not_base_quality = -10 * np.log10(not_p) # Back to Phred

	for b in base_dict.keys():
		if b == indel_char:
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
		elif b == indel_char:
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


def predictions_to_snp_scores(predictions):
	snp = predictions[:, snp_indel_labels['SNP']]
	not_snp = predictions[:, snp_indel_labels['NOT_SNP']]
	return np.log(eps + snp / (not_snp + eps))


def predictions_to_indel_scores(predictions):
	indel = predictions[:, snp_indel_labels['INDEL']]
	not_indel = predictions[:, snp_indel_labels['NOT_INDEL']]
	return np.log(eps + indel / (not_indel + eps))


def predictions_to_snp_indel_scores(predictions):
	snp_dict = predictions_to_snp_scores(predictions)
	indel_dict = predictions_to_indel_scores(predictions)
	return snp_dict, indel_dict


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Arguments ~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_args():
	parser = argparse.ArgumentParser()

	# Tensor defining dictionaries
	parser.add_argument('--tensor_map', default='read_tensor',
						help='Key which looks up the map from tensor channels to their meaning.')
	parser.add_argument('--input_symbols', default=dna_symbols,
						help='Dict mapping input symbols to their index within input tensors.')
	parser.add_argument('--labels', default=snp_indel_labels,
						help='Dict mapping label names to their index within label tensors.')

	# Tensor defining arguments
	parser.add_argument('--batch_size', default=32, type=int,
						help='Mini batch size for stochastic gradient descent algorithms.')
	parser.add_argument('--read_limit', default=128, type=int,
						help='Maximum number of reads to load.')
	parser.add_argument('--window_size', default=128, type=int,
						help='Size of sequence window to use as input, typically centered at a variant.')
	parser.add_argument('--channels_last', default=False, dest='channels_last', action='store_true',
						help='Store the channels in the last axis of tensors, tensorflow->true, theano->false')
	parser.add_argument('--base_quality_mode', default='phot', choices=['phot', 'phred', '1hot'],
						help='How to treat base qualities, must be in [phot, phred, 1hot]')

	# Annotation arguments
	parser.add_argument('--annotations', help='Array of annotation names, initialised via annotation_set argument')
	parser.add_argument('--annotation_set', default='gatk', choices=annotations.keys(),
						help='Key which maps to an annotations list (or None for architectures that do not take annotations).')

	# Genomic area of interest
	parser.add_argument('--start_pos', default=0, type=int,
						help='Genomic position to start from, helpful for parallelization.')
	parser.add_argument('--end_pos', default=0, type=int,
						help='Genomic position to end at, helpful for parallelization.')
	parser.add_argument('--chrom', help='Specific chromosome to load, helpful for parallelization.')

	# Input files and directories: vcfs, bams, beds, hd5, fasta
	parser.add_argument('--architecture', default='',
						help='A json file specifying semantics and architecture of a neural net.')
	parser.add_argument('--bam_file',
						help='Path to a BAM file to train from or generate tensors with.')
	parser.add_argument('--train_vcf',
						help='Path to a VCF that has verified true calls from NIST, platinum genomes, etc.')
	parser.add_argument('--negative_vcf',
						help='Haplotype Caller or VQSR generated VCF with raw annotation values [and quality scores].')
	parser.add_argument('--output_vcf', default=None,
						help='Optional VCF to write to.')
	parser.add_argument('--bed_file',
						help='Bed file specifying high confidence intervals associated with args.train_vcf.')
	parser.add_argument('--dataset',
						help='Directory of tensors, must be split into test/valid/train sets with directories for each label within.')
	parser.add_argument('--reference_fasta',
						help='The reference FASTA file (e.g. HG19 or HG38).')

	# Run specific arguments
	parser.add_argument('--id', default='no_id',
						help='Identifier for this run, user-defined string to keep experiments organized.')
	parser.add_argument('--random_seed', default=12878, type=int,
						help='Random seed to use throughout run.  Always use np.random.')

	# Parse, print, set annotations and seed
	args = parser.parse_args()
	args.annotations = annotations_from_args(args)
	np.random.seed(args.random_seed)
	print('Arguments are', args)

	return args


def annotations_from_args(args):
	if args.annotation_set and args.annotation_set in annotations:
		return annotations[args.annotation_set]
	return None


def tensor_shape_from_args(args):
	in_channels = len(get_tensor_channel_map_from_args(args))
	if K.image_data_format() == 'channels_last':
		tensor_shape = (args.read_limit, args.window_size, in_channels)
	else:
		tensor_shape = (in_channels, args.read_limit, args.window_size)
	return tensor_shape


def set_args_and_get_model_from_semantics(args, semantics_json):
	'''Recreate a model from a json file specifying model semantics.

	Update the args namespace from the semantics file values.
	Assert that the serialized tensor map and the recreated one are the same.

	Arguments:
		args.tensor_map: String which indicates tensor map to use or None
		args.window_size: sites included in the tensor map
		args.read_limit: Maximum reads included in the tensor map
		args.annotations: List of annotations or None
		semantics_json: Semantics json file (created with serialize_model_semantics())

	Returns:
		The Keras model
	'''
	with open(semantics_json, 'r') as infile:
		semantics = json.load(infile)

	if 'input_tensor_map' in semantics:
		args.tensor_map = semantics['input_tensor_map_name']
		args.window_size = semantics['window_size']
		args.read_limit = semantics['read_limit']
		tm = get_tensor_channel_map_from_args(args)
		assert(len(tm) == len(semantics['input_tensor_map']))
		#print('tm len:', len(tm), 'semantics map:', semantics['input_tensor_map'])
		for key in tm:
			assert(tm[key] == semantics['input_tensor_map'][key])

	if 'input_annotations' in semantics:
		args.annotations = semantics['input_annotations']

	args.input_symbols = semantics['input_symbols']
	args.labels = semantics['output_labels']

	if 'data_dir' in semantics:
		args.data_dir = semantics['data_dir']

	weight_path_hd5 = os.path.join(os.path.dirname(semantics_json), semantics['architecture'])
	print('Loading keras weight file from:', weight_path_hd5)
	model = load_model(weight_path_hd5, custom_objects=get_metric_dict(args.labels))
	model.summary()
	return model


def args_and_model_from_semantics(semantics_json):
	args = parse_args()
	model = set_args_and_get_model_from_semantics(args, semantics_json)
	return args, model


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Tensor Maps ~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_tensor_channel_map_from_args(args):
	'''Return tensor mapping dict given args.tensor_map'''
	if not args.tensor_map:
		return None

	if 'read_tensor' == args.tensor_map:
		return get_tensor_channel_map()
	elif 'reference' == args.tensor_map:
		return get_tensor_channel_map_1d_dna()
	else:
		raise ValueError('Unknown tensor mapping mode:', args.tensor_map)


def get_tensor_channel_map_1d_dna():
	'''1D Reference tensor with 4 channel DNA encoding.'''
	tensor_map = {}
	for k in dna_symbols.keys():
		tensor_map[k] = dna_symbols[k]

	return tensor_map


def get_tensor_channel_map_reference_reads():
	'''Read and reference tensor with 4 channel DNA encoding.
	Plus insertions and deletions.
	'''
	tensor_map = {}
	for k in inputs_indel.keys():
		tensor_map['read_'+k] = inputs_indel[k]
	for k in inputs_indel.keys():
		tensor_map['reference_'+k] = len(inputs_indel) + inputs_indel[k]

	return tensor_map


def get_tensor_channel_map():
	'''Read and reference tensor with 4 channel DNA encoding.
	Also includes read flags.
	'''
	tensor_map = {}
	for k in inputs_indel.keys():
		tensor_map['read_'+k] = inputs_indel[k]
	for k in inputs_indel.keys():
		tensor_map['reference_'+k] = len(inputs_indel) + inputs_indel[k]
	tensor_map['flag_bit_4'] = 10
	tensor_map['flag_bit_5'] = 11
	tensor_map['flag_bit_6'] = 12
	tensor_map['flag_bit_7'] = 13
	tensor_map['mapping_quality'] = 14
	return tensor_map


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~ Metrics ~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def precision(y_true, y_pred):
	'''Calculates the precision, a metric for multi-label classification of
	how many selected items are relevant.
	'''
	true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
	predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
	precision = true_positives / (predicted_positives + K.epsilon())
	return precision


def recall(y_true, y_pred):
	'''Calculates the recall, a metric for multi-label classification of
	how many relevant items are selected.
	'''
	true_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)))
	possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
	recall = true_positives / (possible_positives + K.epsilon())
	return recall


def per_class_recall(labels):
	recall_fxns = []

	for label_key in labels:
		label_idx = labels[label_key]
		string_fxn = 'def '+ label_key + '_recall(y_true, y_pred):\n'
		string_fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
		string_fxn += '\tpossible_positives = K.sum(K.round(K.clip(y_true, 0, 1)), axis=0)\n'
		string_fxn += '\treturn true_positives['+str(label_idx)+'] / (possible_positives['+str(label_idx)+'] + K.epsilon())\n'

		exec(string_fxn)
		recall_fxn = eval(label_key + '_recall')
		recall_fxns.append(recall_fxn)

	return recall_fxns


def per_class_precision(labels):
	precision_fxns = []

	for label_key in labels:
		label_idx = labels[label_key]
		string_fxn = 'def '+ label_key + '_precision(y_true, y_pred):\n'
		string_fxn += '\ttrue_positives = K.sum(K.round(K.clip(y_true*y_pred, 0, 1)), axis=0)\n'
		string_fxn += '\tpredicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)), axis=0)\n'
		string_fxn += '\treturn true_positives['+str(label_idx)+'] / (predicted_positives['+str(label_idx)+'] + K.epsilon())\n'

		exec(string_fxn)
		precision_fxn = eval(label_key + '_precision')
		precision_fxns.append(precision_fxn)

	return precision_fxns


def get_metric_dict(labels=snp_indel_labels):
	metrics = {'precision':precision, 'recall':recall}
	precision_fxns = per_class_precision(labels)
	recall_fxns = per_class_recall(labels)
	for i,label_key in enumerate(labels.keys()):
		metrics[label_key+'_precision'] = precision_fxns[i]
		metrics[label_key+'_recall'] = recall_fxns[i]

	return metrics