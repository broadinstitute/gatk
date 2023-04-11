# -*- coding: utf-8 -*-
import numpy as np
import unittest

from curate_input_array_files import curate_input_arrays

dir='curate_input_array_test_files/'
with open(dir + "output_sample_name_list_file_correct"):
    output_sample_name_list_correct = np.loadtxt(dir + "output_sample_name_list_file_correct", dtype=str).tolist()
with open(dir + "output_vcf_list_file_correct"):
    output_vcf_list_correct = np.loadtxt(dir + "output_vcf_list_file_correct", dtype=str).tolist()
with open(dir + "output_vcf_index_list_file_correct"):
    output_vcf_index_list_correct = np.loadtxt(dir + "output_vcf_index_list_file_correct", dtype=str).tolist()

class TestCurateInputArrays(unittest.TestCase):
    def test_curate_input_array_files_success(self):
        actual = curate_input_arrays(
            sample_map_to_be_loaded_file_name=dir + 'input_samples_to_be_loaded_map_file',
            sample_name_list_file_name=dir + 'input_sample_name_list_file',
            vcf_list_file_name=dir + 'input_vcf_list_file',
            vcf_index_list_file_name=dir + 'input_vcf_index_list_file',
            output_files='')
        self.maxDiff=None
        self.assertEqual(actual['sample_names_array'], output_sample_name_list_correct)
        self.assertEqual(actual['vcf_array'], output_vcf_list_correct)
        self.assertEqual(actual['vcf_indexes_array'], output_vcf_index_list_correct)
