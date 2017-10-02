#!/usr/bin/env python
# ApplyNeuralNetModel.py
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
from collections import Counter

# Keras Imports
from keras.models import Model, save_model, load_model

# Defines
eps = 1e-7
tensor_exts = ['.h5', '.hd5']
snp_indel_labels = {'NOT_SNP':0, 'NOT_INDEL':1, 'SNP':2, 'INDEL':3}


def run():
	args = parse_args()
	
	# Load the model
	model = load_model(args.architecture)
	model.summary()

	# Get the tensors
	tensors, annotations, positions = load_tensors_and_annotations_from_class_dirs(args)
	
	# Make predictions
	predictions = model.predict([tensors, annotations], batch_size=args.batch_size)	
	snp_dict = predictions_to_snp_scores(args, predictions, positions)
	indel_dict = predictions_to_indel_scores(args, predictions, positions)
	print('SNP Predictions:\n', snp_dict, '\n\nINDEL Predictions:\n', indel_dict)

def parse_args():
	parser = argparse.ArgumentParser()

	parser.add_argument('--vcf', help='Path to a VCF to annotate with Neural Network Scores')
	parser.add_argument('--architecture', help='hd5 file of with model architecture and weights.')
	parser.add_argument('--tensors', help='Directory containing tensors to evaluate with the model.')
	parser.add_argument('--labels', default=snp_indel_labels, help='Map label names to their index within label tensors.')
	parser.add_argument('--batch_size', default=4, type=int, help='Mini batch size for stochastic gradient descent algorithms.')
	
	# Parse, print
	args = parser.parse_args()
	print('Arguments are', args)

	return args


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
		except ValueError, e:
			print(str(e), '\nValue error at:', tp)

	return (np.asarray(tensors), np.asarray(annotations), np.asarray(positions))


if __name__=='__main__':
	run()

