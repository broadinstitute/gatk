import unittest
import gzip
import tempfile

from create_genes_bqloadjson_from_annotations import make_annotation_json

dir='variant_annotation_table_test_files/'
normal_annotated_genes_json=dir + "test_annotated_genes.json.gz"
expected_output_annotated_genes_bq_json=dir + "expected_test_output_annotated_genes_bq.json.gz"

# An old style (before we separated positions and genes jsons into separate files) - expected to fail
old_annotated_genes_json=dir + "test_old_annotated_genes.json.gz"

class TestCreateGenesBqloadjsonFromAnnotations(unittest.TestCase):
    def test_normal(self):
        output_json = tempfile.NamedTemporaryFile(prefix="test_output_annotated_genes_bq", suffix=".json.gz")
        make_annotation_json(normal_annotated_genes_json, output_json.name)
        contents = gzip.open(output_json, 'r').read()
        expected_contents = gzip.open(expected_output_annotated_genes_bq_json, 'r').read()
        self.assertEqual(contents, expected_contents)

    # This test tries to read an old-style json file and verifies that the program will fail if fed such a file
    @unittest.expectedFailure
    def test_old(self):
        output_json = tempfile.NamedTemporaryFile(prefix="test_old_output_annotated_genes_bq", suffix=".json.gz")
        make_annotation_json(old_annotated_genes_json, output_json.name)
