import unittest

from create_vat_inputs import parse_ancestry_file


class TestHailCreateVatInputs(unittest.TestCase):

    def test_parse_ancestry_file(self):
        # TODO more interesting test data: this is 10 samples all with 'eas' ancestry.

        with open('create_vat_inputs_test_files/quickstart_ancestry.tsv') as ancestry:
            expected = {
                'ERS4367795': 'eas',
                'ERS4367796': 'eas',
                'ERS4367797': 'eas',
                'ERS4367798': 'eas',
                'ERS4367799': 'eas',
                'ERS4367800': 'eas',
                'ERS4367801': 'eas',
                'ERS4367803': 'eas',
                'ERS4367804': 'eas',
                'ERS4367805': 'eas',
            }
            actual = parse_ancestry_file(ancestry)
            self.assertEqual(actual, expected)

    # TODO more tests
