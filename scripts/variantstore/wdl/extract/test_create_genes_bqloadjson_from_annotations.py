import unittest
import gzip
import tempfile
from pathlib import Path
import io
import contextlib
import logging
from importlib import reload

from create_genes_bqloadjson_from_annotations import make_annotation_json

dir='variant_annotation_table_test_files/'

class TestCreateGenesBqloadjsonFromAnnotations(unittest.TestCase):
    def test_normal(self):
        normal_annotated_genes_json=dir + "test_annotated_genes.json"
        expected_output_annotated_genes_bq_json=dir + "expected_test_output_annotated_genes_bq.json"

        output_json = tempfile.NamedTemporaryFile(prefix="test_output_annotated_genes_bq", suffix=".json.gz")
        logging.basicConfig(
            format='%(asctime)s %(levelname)-8s %(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S')
        make_annotation_json(normal_annotated_genes_json, output_json.name, logging)
        contents = gzip.open(output_json, 'rt').read()
        expected_contents = Path(expected_output_annotated_genes_bq_json).read_text()
        self.assertEqual(contents, expected_contents)

    # This test reads an empty genes json file
    def test_empty(self):
        empty_annotated_genes_json=dir + "test_empty_annotated_genes.json"

        output_json = tempfile.NamedTemporaryFile(prefix="test_empty_output_annotated_genes_bq", suffix=".json.gz")

        f = io.StringIO()
        with self.assertRaises(SystemExit) as cm, contextlib.redirect_stderr(f):
            # logging runs within its own thread, we need to configure/instantiate it here so the test can get at its stderr
            # but, if a previous test has set up logging already, then its already running in its own thread
            # So we shut down and reload it for the test to work.
            logging.shutdown()
            reload(logging)
            logging.basicConfig(
                format='%(asctime)s %(levelname)-8s %(message)s',
                level=logging.INFO,
                datefmt='%Y-%m-%d %H:%M:%S')
            make_annotation_json(empty_annotated_genes_json, output_json.name, logging)
        self.assertEqual(cm.exception.code, 0)  # We expect it to pass, with a warning
        # Print the output from the method for understanding the (expected) failure.
        print(f.getvalue())

        self.assertTrue("WARNING: Found no items in annotated json file" in f.getvalue())

    # This test tries to read an old-style json file and verifies that the program will fail if fed such a file
    def test_old(self):
        # An old style (before we separated positions and genes jsons into separate files) - expected to fail
        old_annotated_genes_json=dir + "test_old_annotated_genes.json"

        output_json = tempfile.NamedTemporaryFile(prefix="test_old_output_annotated_genes_bq", suffix=".json.gz")

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
            make_annotation_json(old_annotated_genes_json, output_json.name, logging)
        self.assertEqual(cm.exception.code, 0)  # We expect it to pass, with a warning
        # Print the output from the method for understanding the (expected) failure.
        print(f.getvalue())

        self.assertTrue("WARNING: Found no items in annotated json file" in f.getvalue())
