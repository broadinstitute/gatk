import json
import unittest
import copy
import gzip
import tempfile
from pathlib import Path
import io
import contextlib
import logging
from importlib import reload

from create_vt_bqloadjson_from_annotations import make_annotation_json, make_annotated_json_row

dir='variant_annotation_table_test_files/'

with open(dir + 'vat_expected_pathogenic.json') as vat_expected_pathogenic, open(dir + 'vat_test_pathogenic.json') as vat_test_pathogenic:
    cysticFibrosisExpectedVAT = json.load(vat_expected_pathogenic)
    cysticFibrosisNirvanaOutput = json.load(vat_test_pathogenic)
with open(dir + 'vat_expected_likely_benign.json') as vat_expected_likely_benign, open(dir + 'vat_test_likely_benign.json') as vat_test_likely_benign:
    likelyBenignExpectedVAT = json.load(vat_expected_likely_benign)
    likelyBenignNirvanaOutput = json.load(vat_test_likely_benign)

class TestCreateVtBqloadjsonFromAnnotations(unittest.TestCase):
    def test_normal(self):
        normal_annotated_positions_json=dir + "test_annotated_positions.json"
        expected_output_annotated_positions_bq_json=dir + "expected_test_output_annotated_positions_bq.json"

        output_json = tempfile.NamedTemporaryFile(prefix="test_output_annotated_positions_bq", suffix=".json.gz")
        logging.basicConfig(
            format='%(asctime)s %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')
        make_annotation_json(normal_annotated_positions_json, output_json.name, logging)
        contents = gzip.open(output_json, 'rt').read()
        expected_contents = Path(expected_output_annotated_positions_bq_json).read_text()

        # contents = Path(output_json).read_text()
        # expected_contents = Path(expected_output_annotated_positions_bq_json).read_text()
        self.assertEqual(contents, expected_contents)

    # This test tries to read an old-style json file and verifies that the program will fail if fed such a file
    def test_old(self):
        # An old style (before we separated positions and genes jsons into separate files) - expected to fail
        old_annotated_positions_json=dir + "test_old_annotated_positions.json"

        output_json = tempfile.NamedTemporaryFile(prefix="test_old_output_annotated_positions_bq", suffix=".json.gz")

        f = io.StringIO()
        with self.assertRaises(SystemExit) as cm, contextlib.redirect_stderr(f):
            # logging runs within its own thread, we need to configure/instantiateit here so the test can get at its stderr
            # but, if a previous test has set up logging already, then its already running in its own thread
            # So we shut down and reload it for the test to work.
            logging.shutdown()
            reload(logging)
            logging.basicConfig(
                format='%(asctime)s %(levelname)-8s %(message)s',
                level=logging.INFO,
                datefmt='%Y-%m-%d %H:%M:%S')
            make_annotation_json(old_annotated_positions_json, output_json.name, logging)
        self.assertEqual(cm.exception.code, 0)  # We expect it to pass with a warning
        # Print the output from the method for understanding the (expected) failure.
        print(f.getvalue())

        self.assertTrue("WARNING: Found no items in annotated json file" in f.getvalue())

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
        clinvar_swap = [{ 'id': 'RCV01',
                          'significance': [ 'sleater', 'kinney', 'guitar', 'solo', 'pathogenic'],
                          'refAllele': 'TCT',
                          'altAllele': '-',
                          'phenotypes': ['brunette'],
                      'lastUpdatedDate': '2020-03-01'},
                        { 'id': 'RCV02',
                          'significance': [ 'pathogenic', 'LikELy paTHoGenIc', 'conflicting data from submitters'],
                          'refAllele': 'TCT',
                          'altAllele': '-',
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
        clinvar_swap = [{ 'id': 'RCV01',
                          'significance': [],
                          'refAllele': 'TCT',
                          'altAllele': '-',
                          'phenotypes': [],
                          'lastUpdatedDate': '2020-03-01'},
                        { 'id': 'nope',
                          'significance': [ 'carrie'],
                          'refAllele': 'TCT',
                          'altAllele': '-',
                          'phenotypes': ['did this go through?'],
                          'lastUpdatedDate': '2020-03-02'},
                        { 'id': 'RCV02',
                          'refAllele': 'T',
                          'altAllele': '-'},
                        { 'id': 'RCV03',
                          'refAllele': 'TCT',
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