# -*- coding: utf-8 -*-
import uuid
import time
from datetime import datetime

import csv
import json
import gzip
import argparse
import unittest
import copy

from create_variant_annotation_table import make_annotated_json_row

dir='variant_annotation_table_test_files/'
with open(dir + 'vat_expected_pathogenic.json') as vat_expected_pathogenic, open(dir + 'vat_test_pathogenic.json') as vat_test_pathogenic:
    cysticFibrosisExpectedVAT = json.load(vat_expected_pathogenic)
    cysticFibrosisNirvanaOutput = json.load(vat_test_pathogenic)
with open(dir + 'vat_expected_likely_benign.json') as vat_expected_likely_benign, open(dir + 'vat_test_likely_benign.json') as vat_test_likely_benign:
    likelyBenignExpectedVAT = json.load(vat_expected_likely_benign)
    likelyBenignNirvanaOutput = json.load(vat_test_likely_benign)


class TestMakeAnnotatedJsonRow(unittest.TestCase):

    def test_make_annotated_json_row_success(self):
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


    def test_clinvar_ordering(self):
        clinvar_swap = [{ 'id': 'RCV01', \
                      'significance': [ 'sleater', 'kinney', 'guitar', 'solo', 'pathogenic'], \
                      'refAllele': 'TCT', \
                      'altAllele': '-', \
                      'phenotypes': ['brunette'],
                      'lastUpdatedDate': '2020-03-01'}, \
                      { 'id': 'RCV02', \
                      'significance': [ 'pathogenic', 'LikELy paTHoGenIc', 'conflicting data from submitters'], \
                      'refAllele': 'TCT', \
                      'altAllele': '-', \
                      'phenotypes': ['blonde'],
                      'lastUpdatedDate': '2020-03-02'}]
        clinvarOrderingNirvanaOutput = copy.deepcopy(cysticFibrosisNirvanaOutput)
        clinvarOrderingNirvanaOutput["clinvar"] = clinvar_swap
        actual = make_annotated_json_row(
            row_position=117559590,
            row_ref="ATCT",
            row_alt="A",
            variant_line=clinvarOrderingNirvanaOutput,
            transcript_line=None)
        self.maxDiff=None
        self.assertEqual(actual.get('clinvar_classification'), ['likely pathogenic', 'pathogenic', 'conflicting data from submitters', 'guitar', 'kinney', 'sleater', 'solo'])
        self.assertEqual(actual.get('clinvar_id'), ['RCV01', 'RCV02'])
        self.assertEqual(actual.get('clinvar_phenotype'), ['blonde', 'brunette'])
        self.assertEqual(actual.get('clinvar_last_updated'), '2020-03-02')


    def test_clinvar_inclusion(self):
        clinvar_swap = [{ 'id': 'RCV01', \
                      'significance': [], \
                      'refAllele': 'TCT', \
                      'altAllele': '-', \
                      'phenotypes': [], \
                      'lastUpdatedDate': '2020-03-01'}, \
                      { 'id': 'nope', \
                       'significance': [ 'carrie'], \
                      'refAllele': 'TCT', \
                      'altAllele': '-', \
                      'phenotypes': ['did this go through?'], \
                      'lastUpdatedDate': '2020-03-02'}, \
                      { 'id': 'RCV02', \
                      'refAllele': 'T', \
                      'altAllele': '-'}, \
                      { 'id': 'RCV03', \
                      'refAllele': 'TCT', \
                      'altAllele': 'G'}]
        clinvarInclusionNirvanaOutput = copy.deepcopy(cysticFibrosisNirvanaOutput)
        clinvarInclusionNirvanaOutput["clinvar"] = clinvar_swap
        actual = make_annotated_json_row(
            row_position=117559590,
            row_ref="ATCT",
            row_alt="A",
            variant_line=clinvarInclusionNirvanaOutput,
            transcript_line=None)
        self.maxDiff=None
        self.assertEqual(len(actual.get('clinvar_id')), 1)
        self.assertFalse(actual.get('clinvar_phenotype')) # Empty lists are false in python
        self.assertFalse(actual.get('clinvar_classification')) # Empty lists are false in python
        self.assertEqual(actual.get('clinvar_last_updated'), '2020-03-01')