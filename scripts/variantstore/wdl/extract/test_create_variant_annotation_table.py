# -*- coding: utf-8 -*-
import uuid
import time
from datetime import datetime

import csv
import json
import gzip
import argparse
import unittest

from create_variant_annotation_table import make_annotated_json_row

class TestMakeAnnotatedJsonRow(unittest.TestCase):
    def test_make_annotated_json_row_success(self):
        dir='variant_annotation_table_test_files/'
        with open(dir + 'vat_expected_pathogenic.json') as vat_expected_pathogenic, open(dir + 'vat_test_pathogenic.json') as vat_test_pathogenic:
            cysticFibrosisExpectedVAT = json.load(vat_expected_pathogenic)
            cysticFibrosisNirvanaOutput = json.load(vat_test_pathogenic)
            actual = make_annotated_json_row(
                row_position=117559590,
                row_ref="ATCT",
                row_alt="A",
                variant_line=cysticFibrosisNirvanaOutput,
                transcript_line=None)
            expected = cysticFibrosisExpectedVAT
            self.maxDiff=None
            self.assertEqual(actual, expected)

    def test_clinvar_success(self):
        dir='variant_annotation_table_test_files/'
        with open(dir + 'vat_expected_likely_benign.json') as vat_expected_likely_benign, open(dir + 'vat_test_likely_benign.json') as vat_test_likely_benign:
            likelyBenignExpectedVAT = json.load(vat_expected_likely_benign)
            likelyBenignNirvanaOutput = json.load(vat_test_likely_benign)
            actual = make_annotated_json_row(
                row_position=5226567,
                row_ref="A",
                row_alt="C",
                variant_line=likelyBenignNirvanaOutput,
                transcript_line=None)
            expected = likelyBenignExpectedVAT
            self.maxDiff=None
            self.assertIsNone(actual.get('clinvar_id')) # This should not exist---there are no Clinvar values with matching alleles
            self.assertIsNone(actual.get('clinvar_phenotype')) # No Clinvar values with matching alleles
            self.assertIsNone(actual.get('clinvar_classification')) # No Clinvar values with matching alleles
            self.assertEqual(actual, expected)
