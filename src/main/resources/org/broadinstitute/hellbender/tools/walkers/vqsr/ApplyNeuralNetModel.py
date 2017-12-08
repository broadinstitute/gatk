#!/usr/bin/env python
# ApplyModel.py
#
# Load pre-trained model from HD5 files and apply to command line supplied tensors.
#
# September 2017
# Sam Friedman 
# sam@broadinstitute.org

# Python 2/3 friendly
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

# Imports
import os
import sys
import vcf
import h5py
import time
import pysam
import argparse
import numpy as np
from Bio import Seq, SeqIO
from collections import Counter, defaultdict

# Keras Imports
import keras.backend as K
from keras.models import load_model

# Defines
eps = 1e-7
tensor_exts = ['.h5', '.hd5']
snp_indel_labels = {'NOT_SNP':0, 'NOT_INDEL':1, 'SNP':2, 'INDEL':3}
total_read_flags = 12

skip_char = '~'
indel_char = '*'


data_path = '/dsde/data/deep/vqsr/'
reference_fasta = data_path + 'Homo_sapiens_assembly19.fasta'

dna_symbols = {'A':0, 'C':1, 'G':2, 'T':3}
inputs_indel = {'A':0, 'C':1, 'G':2, 'T':3, indel_char:4}
input_symbols = dna_symbols

annotations = ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum']

cigar_code = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4}

# Base calling ambiguities, See https://www.bioinformatics.org/sms/iupac.html
ambiguity_codes = {'K':[0,0,0.5,0.5], 'M':[0.5,0.5,0,0], 'R':[0.5,0,0,0.5], 'Y':[0,0.5,0.5,0], 'S':[0,0.5,0,0.5], 'W':[0.5,0,0.5,0], 
				  'B':[0,0.333,0.333,0.334], 'V':[0.333,0.333,0,0.334],'H':[0.333,0.333,0.334,0],'D':[0.333,0,0.333,0.334],
				  'X':[0.25,0.25,0.25,0.25], 'N':[0.25,0.25,0.25,0.25]}


def run():
	args = parse_args()

	# Load the model
	model = load_model(args.architecture, custom_objects=get_metric_dict(args.labels))
	model.summary()

	if '1d_as' == args.mode:
		annotations = ['MQ', 'DP', 'SOR', 'FS', 'QD', 'MQRankSum', 'ReadPosRankSum', 'DP_MEDIAN', 'DREF_MEDIAN', 'GQ_MEDIAN', 'AB_MEDIAN']
		annotate_vcf_1d_per_allele(args, model)
	elif '1d' == args.mode:
		annotate_vcf_1d(args, model)
	elif '2d' == args.mode:
		input_symbols = inputs_indel
		annotate_vcf_2d(args, model)
	else:
		print('Unknown mode:', args.mode)


def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('mode')

	parser.add_argument('--architecture', help='hd5 file with model architecture and weights.')
	parser.add_argument('--input_vcf', help='Path to a VCF to read and annotate with Neural Network Scores')
	parser.add_argument('--output_vcf', default='./out.vcf', help='Path to a VCF to write annotated with Neural Network Scores')
	parser.add_argument('--tensors', help='Directory containing tensors to evaluate with the model.')
	parser.add_argument('--labels', default=snp_indel_labels, help='Map label names to their index within label tensors.')
	parser.add_argument('--batch_size', default=32, type=int, help='Mini batch size for stochastic gradient descent algorithms.')
	parser.add_argument('--reference_fasta', default=reference_fasta, help='The reference FASTA file (e.g. HG19 or HG38).')
	parser.add_argument('--exclude_annotations', default=False, action='store_true', help='Exclude annotations in input tensor.')
	parser.add_argument('--exclude_dna', default=False, action='store_true', help='Exclude sequence data (reference and/or reads) in input tensor.')
	parser.add_argument('--read_limit', default=128, type=int, help='Maximum number of reads to load.')	
	parser.add_argument('--window_size', default=128, type=int, help='Size of sequence window to use as input, typically centered at a variant.')
	parser.add_argument('--bam_file', help='Path to a BAM file to train from or generate tensors with.')
	parser.add_argument('--interval_list', help='Path to a picard style interval list file.')
	parser.add_argument('--tensor_map', default='2d_mapping_quality', help='Key which looks up the map from tensor channels to their meaning.')
	parser.add_argument('--channels_last', default=False, dest='channels_last', action='store_true', help='Store the channels in the last axis of tensors, tensorflow->true, theano->false')	
	parser.add_argument('--base_quality_mode', default='phot', choices=['phot', 'phred', '1hot'], help='How to treat base qualities, must be in {phot, phred, 1hot}')
	

	# Parse, print
	args = parser.parse_args()
	print('Arguments are', args)

	return args


def predict_on_tensors(args, model):
	# Get the tensors
	tensors, annotations, positions = load_tensors_and_annotations_from_class_dirs(args)
	
	# Make predictions
	predictions = model.predict([tensors, annotations], batch_size=args.batch_size)	
	snp_dict = predictions_to_snp_scores(args, predictions, positions)
	indel_dict = predictions_to_indel_scores(args, predictions, positions)
	print('SNP Predictions:\n', snp_dict, '\n\nINDEL Predictions:\n', indel_dict)


def annotate_vcf_1d(args, model):
	stats = Counter()

	score_key = "CNN_1D"

	vcf_reader = pysam.VariantFile(args.input_vcf, "r")
	vcf_reader.header.info.add(score_key, "1", "Float", "Site-level Score from 1D Convolutional Neural Net.")
	vcf_writer = pysam.VariantFile(args.output_vcf, 'w', header=vcf_reader.header)
	print('got vcfs.')

	intervals = interval_file_to_dict(args.interval_list)
	print('got intervals.')

	reference = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, "fasta"))
	print('got ref.')

	tensor_channel_map = get_tensor_channel_map_from_args(args)

	positions = []
	variant_batch = []
	time_batch = time.time()
	annotation_batch = np.zeros((args.batch_size, len(annotations)))
	dna_batch = np.zeros((args.batch_size, args.window_size, len(input_symbols)))
	
	print('iterate over intervals...')
	for k in intervals:
		for start,stop in zip(intervals[k][0], intervals[k][1]):
			for variant in vcf_reader.fetch(k, start, stop):
				idx_offset, ref_start, ref_end = get_variant_window(args, variant)

				contig = reference[variant.contig]
				record = contig[variant.pos-idx_offset: variant.pos+idx_offset]

				if not args.exclude_annotations:
					if all(map(lambda x: x not in variant.info and x != "QUAL", annotations)):
						stats['Missing ALL annotations'] += 1
						continue # Require at least 1 annotation...

					annotation_data = np.zeros(( len(annotations), ))
					qual_and_dp_normalizer = 1000000.0

					for i,a in enumerate(annotations):
						if a == "QUAL":
							annotation_data[i] = variant.qual
						elif a in variant.info:
							annotation_data[i] = variant.info[a]

						if a == "DP" or a == "QUAL":
							 annotation_data[i] /= qual_and_dp_normalizer

					annotation_batch[stats['annotations_in_batch']] = annotation_data
					stats['annotations_in_batch'] += 1

				if not args.exclude_dna:
					dna_data = np.zeros( (args.window_size, len(input_symbols)) )
					assert(len(record.seq) == args.window_size)
					for i,b in enumerate(record.seq):
						if b in input_symbols:
							dna_data[i, input_symbols[b]] = 1.0
						elif b in ambiguity_codes:
							dna_data[i] = ambiguity_codes[b]
						else:
							print('Error! Unknown code:', b)
							stats['Error! Unknown input symbol:'+b] += 1
							continue

					dna_batch[stats['reference_tensors_in_batch']] = dna_data
					stats['reference_tensors_in_batch'] += 1

				positions.append(variant.contig + '_' + str(variant.pos))
				variant_batch.append(variant)

				if stats['reference_tensors_in_batch'] == args.batch_size or stats['annotations_in_batch'] == args.batch_size:
					input_data = [dna_batch, annotation_batch]
					apply_model_to_batch(args, model, input_data, variant_batch, vcf_writer, score_key, stats)

					# Reset the batch
					positions = []
					variant_batch = []
					dna_batch = np.zeros((args.batch_size, args.window_size, len(input_symbols)))
					annotation_batch = np.zeros((args.batch_size,len(annotations)))
					stats['reference_tensors_in_batch'] = stats['annotations_in_batch'] = 0

	if stats['reference_tensors_in_batch'] > 0 or stats['annotations_in_batch'] > 0: # get scores for the remaining variants
		input_data = [dna_batch, annotation_batch]
		apply_model_to_batch(args, model, input_data, variant_batch, vcf_writer, score_key, stats)
	
	for s in stats.keys():
		print(s, 'has:', stats[s])




def annotate_vcf_1d_per_allele(args, model):
	stats = Counter()

	score_key = 'AS_CNN'
	vcf_reader = pysam.VariantFile(args.input_vcf, "r")
	vcf_reader.header.info.add(score_key, '.', 'Float', 'Allele specific score from 1D Convolutional Neural Net.')
	vcf_writer = pysam.VariantFile(args.output_vcf, 'w', header=vcf_reader.header)

	print('got vcfs.')

	reference = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, "fasta"))
	print('got ref.')

	tensor_channel_map = get_tensor_channel_map_from_args(args)

	positions = []
	variant_batch = []
	time_batch = time.time()
	dna_batch = np.zeros((args.batch_size, args.window_size, len(input_symbols)))
	annotation_batch = np.zeros((args.batch_size, len(annotations)))
	
	print('iterate over vcf...')
	for variant in vcf_reader:
		for allele_index, allele in enumerate(variant.alts):
			idx_offset, ref_start, ref_end = get_variant_window(args, variant)

			contig = reference[variant.contig]
			record = contig[variant.pos-idx_offset: variant.pos+idx_offset]

			if not args.exclude_annotations:
				if all(map(lambda x: x not in variant.info and x != "QUAL", annotations)):
					stats['Missing ALL annotations'] += 1
					continue # Require at least 1 annotation...

				annotation_data = np.zeros(( len(annotations), ))

				qual_and_dp_normalizer = 1000000.0
				as_normalizer = 100.0

				for i,a in enumerate(annotations):
					if a in ['DP_MEDIAN', 'DREF_MEDIAN', 'GQ_MEDIAN', 'AB_MEDIAN']:
						annotation_data[i] = variant.info[a][allele_index]
					elif a == "QUAL":				
						annotation_data[i] = variant.qual
					elif a in variant.info:
						annotation_data[i] = variant.info[a]

					if a == "DP" or a == "QUAL":
						 annotation_data[i] /= qual_and_dp_normalizer	
		
					if a == "DP_MEDIAN" or a == "GQ_MEDIAN":
						 annotation_data[i] /= as_normalizer

				annotation_batch[stats['annotations_in_batch']] = annotation_data
				stats['annotations_in_batch'] += 1
					
			if not args.exclude_dna:
				dna_data = np.zeros( (args.window_size, len(input_symbols)) )
				assert(len(record.seq) == args.window_size)
				for i,b in enumerate(record.seq):
					if b in input_symbols:
						dna_data[i, input_symbols[b]] = 1.0
					elif b in ambiguity_codes:
						dna_data[i] = ambiguity_codes[b]
					else:
						print('Error! Unknown code:', b)
						stats['Error! Unknown input symbol:'+b] += 1
						continue
		
				dna_batch[stats['reference_tensors_in_batch']] = dna_data
				stats['reference_tensors_in_batch'] += 1
		
			positions.append(variant.contig + '_' + str(variant.pos) + '_' + str(allele_index))
			variant_batch.append((variant, allele_index))
			if stats['reference_tensors_in_batch'] == args.batch_size or stats['annotations_in_batch'] == args.batch_size:
				break

		if stats['reference_tensors_in_batch'] >= args.batch_size or stats['annotations_in_batch'] >= args.batch_size:
			t0 = time.time()
			predictions = model.predict([dna_batch, annotation_batch], batch_size=args.batch_size)	
			snp_dict = predictions_to_snp_scores(args, predictions, positions)
			indel_dict = predictions_to_indel_scores(args, predictions, positions)
			t1 = time.time()
			time_batch = time.time()
			if (stats['variants_written']) % (args.batch_size*50) == 0:
				print('CNN Batch predictions took:',(t1-t0), 'seconds. \nWhole batch time:', (t1-time_batch), 'Time per variant:', ((t1-time_batch)/args.batch_size))
								
			# loop over the batch of variants and add per allele scores
			scores = defaultdict(list)
			variants_to_write = {}
			for v_out, allele_index in variant_batch:
				g_position = v_out.contig + '_' + str(v_out.pos)
				a_position = g_position + '_' + str(allele_index)
				variants_to_write[g_position] = v_out
				if len(v_out.ref) == 1 and len(v_out.alts[allele_index]) == 1: # SNP means ref and alt both are length 1
					scores[g_position].append(float(snp_dict[a_position]))
				else:
					scores[g_position].append(float(indel_dict[a_position]))
				
			for k in sorted(variants_to_write.iterkeys()):
				variants_to_write[k].info[score_key] = tuple(scores[k])
				vcf_writer.write(variants_to_write[k])
				stats['variants_written'] += 1
				if (stats['variants_written']) % (args.batch_size*50) == 0:
					print('SNP Predictions:\n', snp_dict, '\n\nINDEL Predictions:\n', indel_dict, '\n\nLast variant:', variants_to_write[k])
					for s in stats.keys():
						print(s, 'has:', stats[s])
			
			# Reset the batch
			positions = []		
			variant_batch = []
			dna_batch = np.zeros((args.batch_size, args.window_size, len(input_symbols)))
			annotation_batch = np.zeros((args.batch_size,len(annotations)))		
			stats['reference_tensors_in_batch'] = stats['annotations_in_batch'] = 0

	for s in stats.keys():
		print(s, 'has:', stats[s])


def annotate_vcf_2d(args, model):
	stats = Counter()
	score_key = 'CNN_2D'

	samfile = pysam.AlignmentFile(args.bam_file, "rb")	
	print('got sam.')

	vcf_reader = pysam.VariantFile(args.input_vcf, "r")
	vcf_reader.header.info.add(score_key, '.', 'Float', 'Allele specific score from 2D Convolutional Neural Net.')
	vcf_writer = pysam.VariantFile(args.output_vcf, 'w', header=vcf_reader.header)
	print('got vcfs.')

	reference = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, "fasta"))
	print('got ref.')

	intervals = interval_file_to_dict(args.interval_list)
	print('got intervals.')

	tensor_channel_map = get_tensor_channel_map_from_args(args)
	if args.channels_last:
		tensor_shape = (args.read_limit, args.window_size, len(tensor_channel_map))
	else:
		tensor_shape = (len(tensor_channel_map), args.read_limit, args.window_size) 

	print('iterate over vcf')

	positions = []
	variant_batch = []
	time_batch = time.time()
	tensor_batch = np.zeros(((args.batch_size,) + tensor_shape))
	annotation_batch = np.zeros((args.batch_size, len(annotations)))
	
	print('iterate over intervals...')
	for k in intervals:
		for start,stop in zip(intervals[k][0], intervals[k][1]):
			for variant in vcf_reader.fetch(k, start, stop):
				#for allele_index, allele in enumerate(variant.alts):
				idx_offset, ref_start, ref_end = get_variant_window(args, variant)

				contig = reference[variant.contig]	
				record = contig[ ref_start : ref_end ]

				if not args.exclude_annotations:
					if all(map(lambda x: x not in variant.info and x != "QUAL", annotations)):
						stats['Missing ALL annotations'] += 1
						continue # Require at least 1 annotation...

					annotation_data = np.zeros(( len(annotations), ))
					for i,a in enumerate(annotations):
						if a == "QUAL":
							annotation_data[i] = variant.qual
						elif a in variant.info:
							annotation_data[i] = variant.info[a]
					annotation_batch[stats['annotations_in_batch']] = annotation_data
					stats['annotations_in_batch'] += 1

				if not args.exclude_dna:
					good_reads, insert_dict = get_good_reads_new(args, samfile, variant)
					reference_seq = record.seq
					for i in sorted(insert_dict.keys(), key=int, reverse=True):
						reference_seq = reference_seq[:i] + indel_char*insert_dict[i] + reference_seq[i:]

					sequences, qualities, mapping_qualities, flags = good_reads_to_arrays(args, good_reads, ref_start, insert_dict)

					if len(sequences) == 0:
						stats['No acceptable aligned reads'] += 1
						continue

					if args.channels_last:
						read_tensor = np.zeros( (args.read_limit, args.window_size, len(tensor_channel_map)) )
						read_tensor[:,:,:10] = reads_to_tensor(args, sequences, qualities, reference_seq)
					else:
						read_tensor = np.zeros( (len(tensor_channel_map), args.read_limit, args.window_size) )
						read_tensor[:10,:,:] = reads_to_tensor(args, sequences, qualities, reference_seq)
					
					add_flags_to_read_tensor(args, read_tensor, tensor_channel_map, flags)
					add_mq_to_read_tensor(args, read_tensor, tensor_channel_map, mapping_qualities)			
					tensor_batch[stats['read_tensors_in_batch']] = read_tensor
					stats['read_tensors_in_batch'] += 1

				positions.append(variant.contig + '_' + str(variant.pos))
				variant_batch.append(variant)

				if stats['read_tensors_in_batch'] == args.batch_size or stats['annotations_in_batch'] == args.batch_size:
					input_data = [tensor_batch, annotation_batch]
					apply_model_to_batch(args, model, input_data, positions, variant_batch, vcf_writer, score_key, stats)

					# Reset the batch
					positions = []		
					variant_batch = []
					tensor_batch = np.zeros(((args.batch_size,)+tensor_shape))
					annotation_batch = np.zeros((args.batch_size,len(annotations)))		
					stats['read_tensors_in_batch'] = stats['annotations_in_batch'] = 0

	if stats['read_tensors_in_batch'] > 0 or stats['annotations_in_batch'] > 0:
		apply_model_to_batch(args, model, [tensor_batch, annotation_batch], positions, variant_batch, vcf_writer, score_key, stats)

	for s in stats.keys():
		print(s, 'has:', stats[s])


def apply_model_to_batch(args, model, input_data, positions, variant_batch, vcf_writer, score_key, stats):
	predictions = model.predict(input_data, batch_size=args.batch_size)
	snp_dict = predictions_to_snp_scores(args, predictions, positions)
	indel_dict = predictions_to_indel_scores(args, predictions, positions)

	# loop over the batch of variants and write them out with a score
	for v_out in variant_batch:
		position = v_out.contig + '_' + str(v_out.pos)

		if len(v_out.ref) == 1 and len(v_out.alleles[1][0]) == 1: # SNP means ref and alt both are length 1
			v_out.info[score_key] = float(snp_dict[position])
		elif len(v_out.ref) > 1 or len(v_out.alleles[1][0]) > 1:
			v_out.info[score_key] = float(indel_dict[position])
		else:
			stats['Not SNP or INDEL'] += 1

		vcf_writer.write(v_out)
		stats['variants_written'] += 1


def get_good_reads_new(args, samfile, variant, sort_by='base'):
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

	for read in samfile.fetch(variant.contig, variant.pos-1, variant.pos+1):

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
			my_insert = 0
			for t in read.cigartuples:
				if t[0] == 1:
					insert_idx = cur_idx - index_dif
					if insert_idx not in insert_dict: 
						insert_dict[insert_idx] = t[1]
					elif insert_dict[insert_idx] < t[1]:
						insert_dict[insert_idx] = t[1]
				if t[0] in [cigar_code['M'], cigar_code['I'], cigar_code['S'], cigar_code['D']]:
					cur_idx += t[1]

		good_reads.append(read)

	if len(good_reads) > args.read_limit:
		good_reads = np.random.choice(good_reads, size=args.read_limit, replace=False).tolist()

	good_reads.sort(key=lambda x: x.reference_start + x.query_alignment_start)
	if sort_by == 'base':
		good_reads.sort(key=lambda x: x.seq[clamp(variant.pos-x.reference_start-1, 0, len(x.seq)-1)])

	return good_reads, insert_dict


def good_reads_to_arrays(args, good_reads, ref_start, insert_dict):
	qualities = []
	sequences = []
	mapping_qualities = []
	
	flags = np.zeros((total_read_flags, args.read_limit))
	no_qual_filler = 0
	
	for i,read in enumerate(good_reads):
		index_dif = ref_start - read.reference_start
		cur_idx = 0
		my_indel_dict = {}
		for t in read.cigartuples:
			my_ref_idx = cur_idx - index_dif
			if t[0] == 1:
				my_indel_dict[my_ref_idx] = insert_dict[my_ref_idx] - t[1]
			elif t[0] == 2:
				my_indel_dict[my_ref_idx] = t[1]
			if t[0] in [cigar_code['M'], cigar_code['I'], cigar_code['S'], cigar_code['D']]:
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
			rseq = skip_char*(-index_dif) + rseq
			rqual = [no_qual_filler]*(-index_dif) + rqual

		for j in sorted(my_indel_dict.keys(), key=int, reverse=True):
			if j < len(rseq):
				rseq = rseq[:j] + (indel_char*my_indel_dict[j]) + rseq[j:]
				rqual = rqual[:j] + ([no_qual_filler]*my_indel_dict[j]) + rqual[j:]

		flags[:, i] = flag_to_array(read.flag)
		mapping_qualities.append(read.mapping_quality)
		sequences.append(rseq)
		qualities.append(rqual)

	return sequences, qualities, mapping_qualities, flags

def clamp(n, minn, maxn):
	return max(min(maxn, n), minn)

def flag_to_array(flag):
	flags = []
	
	for i in range(total_read_flags):
		flags.append((flag>>i)&1)

	return np.array(flags)


def add_flags_to_read_tensor(args, tensor, tensor_channel_map, flags):
	for k in tensor_channel_map.keys():
		if 'flag' in k:
			flag_bit = int(k.split('_')[-1])
			for read_idx in range(flags.shape[1]):
				if args.channels_last:
					tensor[read_idx, :, tensor_channel_map[k]] = flags[flag_bit, read_idx]
				else:
					tensor[tensor_channel_map[k], read_idx, :] = flags[flag_bit, read_idx]


def add_mq_to_read_tensor(args, tensor, tensor_channel_map, mapping_qualities):
	if not 'mapping_quality' in tensor_channel_map:
		return
		
	mq_normalizer = 60.0
	for read_idx, mq in enumerate(mapping_qualities):
		if args.channels_last:
			tensor[read_idx, :, tensor_channel_map['mapping_quality']] = float(mq) / mq_normalizer
		else:
			tensor[tensor_channel_map['mapping_quality'], read_idx, :] = float(mq) / mq_normalizer


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
		print('Error! Unknown base quality mode:', args.base_quality_mode)


def reads_to_tensor(args, sequences, qualities=None, reference_seq=None):
	debug = False

	in_channels = get_reference_and_read_channels(args)
	if args.channels_last:
		tensor = np.zeros( (args.read_limit, args.window_size, in_channels) )
	else:
		tensor = np.zeros( (in_channels, args.read_limit, args.window_size) )

	for j,l in enumerate(sequences):
		for i,b in enumerate(l):
			if i == args.window_size:
				break
			if b in input_symbols:
				if b == indel_char or qualities is None:
					if args.channels_last:
						tensor[j, i, input_symbols[b]] = 1.0
					else:
						tensor[input_symbols[b], j, i] = 1.0
				else:
					hot_array = quality_from_mode(args, qualities[j][i], b, input_symbols)
					if args.channels_last:
						tensor[j, i, :4] = hot_array
					else:
						tensor[:4, j, i] = hot_array

			elif b in ambiguity_codes:
				if args.channels_last:
					tensor[j, i, :4] = ambiguity_codes[b]
				else:
					tensor[:4, j, i] = ambiguity_codes[b]
			elif b == skip_char:
				continue
			else:
				print('Error! Unknown symbol in seq block:', b)
				return

	if reference_seq:
		for i,b in enumerate(reference_seq):
			if i == args.window_size:
				break
			ref_offset = len(set(input_symbols.values()))
			if b in input_symbols:
				if args.channels_last:
					tensor[:, i, ref_offset+input_symbols[b]] = 1.0
				else:
					tensor[ref_offset+input_symbols[b], :, i] = 1.0
			elif b in ambiguity_codes:
				if args.channels_last:
					tensor[:, i, ref_offset:ref_offset+4] = np.transpose(np.tile(ambiguity_codes[b], (args.read_limit, 1)))
				else:
					tensor[ref_offset:ref_offset+4, :, i] = np.transpose(np.tile(ambiguity_codes[b], (args.read_limit, 1)))		

	if debug:
		np.set_printoptions(threshold=np.inf)
		print(reference_seq, '<- reference sequence')
		print(sequences, '\n\n\n ~~~~~~ Becomes Tensor: ~~~~~~~ \n\n')
		print(tensor)

	return tensor



def get_variant_window(args, variant):
	index_offset = (args.window_size//2)
	reference_start = variant.pos-(index_offset+1)
	reference_end = variant.pos+index_offset

	return index_offset, reference_start, reference_end


def predictions_to_snp_scores(args, predictions, positions):
	snp = predictions[:, args.labels['SNP']]
	not_snp = predictions[:, args.labels['NOT_SNP']]
	snp_scores = np.log(eps + snp / (not_snp + eps))
	return dict(zip(positions, snp_scores))


def predictions_to_indel_scores(args, predictions, positions):
	indel = predictions[:, args.labels['INDEL']]
	not_indel = predictions[:, args.labels['NOT_INDEL']]
	indel_scores = np.log(eps + indel / (not_indel + eps))
	return dict(zip(positions, indel_scores))


def predictions_to_snp_indel_scores(args, predictions, positions):
	snp_dict = predictions_to_snp_scores(args, predictions, positions)
	indel_dict = predictions_to_indel_scores(args, predictions, positions)
	return snp_dict, indel_dict


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


def get_metric_dict(labels):
	metrics = {'precision':precision, 'recall':recall}
	precision_fxns = per_class_precision(labels)
	recall_fxns = per_class_recall(labels)
	for i,label_key in enumerate(labels.keys()):
		metrics[label_key+'_precision'] = precision_fxns[i]
		metrics[label_key+'_recall'] = recall_fxns[i]

	return metrics



def total_input_channels_from_args(args):
	'''Get the number of channels in the tensor map'''		
	return len(get_tensor_channel_map_from_args(args))


def get_reference_and_read_channels(args):
	'''Get the number of read and reference channels in the tensor map'''		
	count = 0
	tm = get_tensor_channel_map_from_args(args)
	for k in tm.keys():
		if 'read' in k or 'reference' in k:
			count += 1
	return count


def get_tensor_channel_map_from_args(args):
	'''Return tensor mapping dict given args.tensor_map'''	
	if '2d_2bit' == args.tensor_map:
		return get_tensor_channel_map_2bit()
	elif '2d_mapping_quality' == args.tensor_map:
		return get_tensor_channel_map_mq()
	elif '2d' == args.tensor_map or '2d_annotations' == args.tensor_map:
		return get_tensor_channel_map()
	elif '1d_dna' == args.tensor_map:
		return get_tensor_channel_map_1d_dna()
	elif '1d' == args.tensor_map or '1d_annotations' == args.tensor_map:
		return get_tensor_channel_map_1d()
	elif args.tensor_map == 'bqsr':
		return bqsr_tensor_channel_map()
	else:
		print('Unknown tensor mapping mode:', args.tensor_map)


def get_tensor_channel_map_1d_dna():
	'''1D Reference tensor with 4 channel DNA encoding.'''
	tensor_map = {}
	for k in inputs.keys():
		tensor_map[k] = inputs[k]
	
	return tensor_map

def get_tensor_channel_map_1d():
	'''1D Reference tensor with 4 channel DNA encoding'''
	tensor_map = {}
	for k in inputs.keys():
		tensor_map[k] = inputs[k]
	
	return tensor_map


def get_tensor_channel_map_1d_plus_beds():
	'''1D Reference tensor with 4 channel DNA encoding.
	Also includes channels for binary labels from bed files.
	'''
	tensor_map = {}
	for k in inputs.keys():
		tensor_map[k] = inputs[k]
	
	ref_offset = len(inputs)
	for i,b in enumerate(reference_beds):
		tensor_map[b] = ref_offset + i
	
	return tensor_map


def get_tensor_channel_map():
	'''Read and reference tensor with 4 channel DNA encoding.
	Also includes read flags.
	'''
	tensor_map = {}
	for k in inputs_indel.keys():
		tensor_map['read_'+k] = inputs_indel[k]
	for k in inputs_indel.keys():
		tensor_map['reference_'+k] = 5 + inputs_indel[k]			
	tensor_map['flag_bit_4'] = 10
	tensor_map['flag_bit_5'] = 11	
	tensor_map['flag_bit_6'] = 12	
	tensor_map['flag_bit_7'] = 13
	return tensor_map


def get_tensor_channel_map_mq():
	'''Read and reference tensor with 4 channel DNA encoding.
	Also includes read flags.
	'''
	tensor_map = {}
	for k in inputs_indel.keys():
		tensor_map['read_'+k] = inputs_indel[k]
	for k in inputs_indel.keys():
		tensor_map['reference_'+k] = 5 + inputs_indel[k]			

	tensor_map['flag_bit_4'] = 10
	tensor_map['flag_bit_5'] = 11	
	tensor_map['flag_bit_6'] = 12	
	tensor_map['flag_bit_7'] = 13
	tensor_map['flag_bit_9'] = 14	
	tensor_map['flag_bit_10'] = 15

	tensor_map['mapping_quality'] = 16

	return tensor_map



def load_tensors_and_annotations_from_class_dirs(args):
	tensors = []
	positions = []
	annotations = []

	for tp in os.listdir(args.tensors):

		fn, file_extension = os.path.splitext(tp)
		if not file_extension.lower() in tensor_exts:
			continue

		gpos = tp.split('-')[2]
		chrom = gpos.split('_')[0]
		pos = os.path.splitext(gpos.split('_')[1])[0]

		positions.append(chrom + '_' + pos)

		try:
			with h5py.File(args.tensors+'/'+tp, 'r') as hf:
				tensors.append(np.array(hf.get('read_tensor')))
				annotations.append(np.array(hf.get('annotations')))
		except ValueError as e:
			print(str(e), '\nValue error at:', tp)

	return (np.asarray(tensors), np.asarray(annotations), np.asarray(positions))



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

	with open(interval_file)as f:
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


if __name__=='__main__':
	run()

