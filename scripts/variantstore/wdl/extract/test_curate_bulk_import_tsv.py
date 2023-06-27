# -*- coding: utf-8 -*-
import numpy as np
import unittest

from curate_bulk_import_tsv import curate_bulk_import_data

test_dir = 'curate_bulk_import_tsv_test_files/'


class TestCurateInputArrays(unittest.TestCase):
    def test_curate_input_array_files_success(self):
        with open(test_dir + "bulk_import_output_expected.tsv") as expected_file:

            expected = np.loadtxt(expected_file, dtype=str)

            actual = curate_bulk_import_data(
                sample_map_to_be_loaded_file_name=test_dir + 'input_samples_to_be_loaded_map.csv',
                bulk_import_input_file_name=test_dir + 'bulk_import_input.tsv')
            self.maxDiff = None

            self.assertEqual(actual, expected)
