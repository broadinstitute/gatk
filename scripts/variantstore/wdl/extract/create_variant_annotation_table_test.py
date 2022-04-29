# -*- coding: utf-8 -*-
import uuid
import time
from datetime import datetime

import csv
import json
#import ijson
import gzip
import argparse
import unittest

from create_variant_annotation_table import make_annotated_json_row

class TestMakeAnnotatedJsonRow(unittest.TestCase):
    def test_make_annotated_json_row_success(self):
        vat_expected_pathogenic = open('variant_annotation_table_test_files/vat_expected_pathogenic.json')
        vat_test_pathogenic = open('variant_annotation_table_test_files/vat_test_pathogenic.json')
        cysticFibrosisExpectedVAT = json.load(vat_expected_pathogenic)
        cysticFibrosisNirvanaOutput = json.load(vat_test_pathogenic)
        actual = make_annotated_json_row(
          row_position=117559590,
          row_ref="ATCT",
          row_alt="A",
          variant_line=cysticFibrosisNirvanaOutput,
          transcript_line=None)
        ## This is dumb....I do this because json and python dicts aren't playing nice
        cysticFibrosisExpectedVAT["gnomad_nfr_ac"] = None
        cysticFibrosisExpectedVAT["gnomad_nfr_an"] = None
        cysticFibrosisExpectedVAT["gnomad_nfr_af"] = None
        cysticFibrosisExpectedVAT["gnomad_failed_filter"] = None
        expected = cysticFibrosisExpectedVAT
        vat_expected_pathogenic.close()
        vat_test_pathogenic.close()
        self.maxDiff=None
        self.assertEqual(actual, expected)


    def test_clinvar_success(self):
        vat_expected_likely_benign = open('variant_annotation_table_test_files/vat_expected_likely_benign.json')
        vat_test_likely_benign = open('variant_annotation_table_test_files/vat_test_likely_benign.json')
        likelyBenignExpectedVAT = json.load(vat_expected_likely_benign)
        likelyBenignNirvanaOutput = json.load(vat_test_likely_benign)
        actual = make_annotated_json_row(
          row_position=5226567,
          row_ref="A",
          row_alt="C",
          variant_line=likelyBenignNirvanaOutput,
          transcript_line=None)
        likelyBenignExpectedVAT["dbsnp_rsid"] = None
        expected = likelyBenignExpectedVAT
        vat_expected_likely_benign.close()
        vat_test_likely_benign.close()
        self.maxDiff=None
        self.assertEqual(actual, expected)

