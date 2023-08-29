# -*- coding: utf-8 -*-
import numpy as np
import unittest

from curate_input_array_files import curate_input_arrays, read_single_column_file

dir='curate_input_array_test_files/'

class TestCurateInputArrays(unittest.TestCase):
    def test_curate_input_array_files_success(self):
        with open(dir + "output_sample_name_list_file_correct") as samples, \
             open(dir + "output_vcf_list_file_correct") as vcfs, \
             open(dir + "output_vcf_index_list_file_correct") as vcf_indexes:

            output_sample_name_list_correct = np.loadtxt(samples, dtype=str).tolist()
            output_vcf_list_correct = np.loadtxt(vcfs, dtype=str).tolist()
            output_vcf_index_list_correct = np.loadtxt(vcf_indexes, dtype=str).tolist()

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

    def test_exome_data(self):
        output_sample_name_list_correct = read_single_column_file(dir + "expected_exome_sample_names.txt")
        output_vcf_list_correct = read_single_column_file(dir + "expected_exome_vcf_list.txt")
        output_vcf_index_list_correct = read_single_column_file(dir + "expected_exome_vcf_index_list.txt")
        actual = curate_input_arrays(
            sample_map_to_be_loaded_file_name=dir + 'exome_sample_map.csv',
            sample_name_list_file_name=dir + 'exome_sample_names.txt',
            vcf_list_file_name=dir + 'exome_vcf_file_names.txt',
            vcf_index_list_file_name=dir + 'exome_vcf_index_file_names.txt',
            output_files='')

        self.maxDiff=None
        self.assertEqual(actual['sample_names_array'], output_sample_name_list_correct)
        self.assertEqual(actual['vcf_array'], output_vcf_list_correct)
        self.assertEqual(actual['vcf_indexes_array'], output_vcf_index_list_correct)
