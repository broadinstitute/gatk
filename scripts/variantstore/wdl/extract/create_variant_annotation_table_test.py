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

vat_expected_pathogenic = open('variant_annotation_table_test_files/vat_expected_pathogenic.json')
vat_test_pathogenic = open('variant_annotation_table_test_files/vat_test_pathogenic.json')
cysticFibrosisExpectedVAT = json.load(vat_expected_pathogenic)
cysticFibrosisNirvanaOutput = json.load(vat_test_pathogenic)

class TestMakeAnnotatedJsonRow(unittest.TestCase):
    def test_make_annotated_json_row_success(self):
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
        self.maxDiff=None
        self.assertEqual(actual, expected)

#likelyBenignPositionObj = "variant_annotation_table_test_files/vat_test_likely_benign.json"

