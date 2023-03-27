import unittest

from strip_custom_annotations_from_sites_only_vcf import parse_annot_header

class TestStripCustomAnnotationsFromSitesOnlyVcf(unittest.TestCase):

    def test_parse_annot_header(self):
        input_custom_annotations_tsv = 'strip_custom_annotations_from_sites_only_vcf_test_files/custom_annotations_template.tsv'

        with open(input_custom_annotations_tsv) as x: expected_header = x.read()
        with open(input_custom_annotations_tsv) as tsv:
            actual_header, actual_column_names, actual_names_to_column_map = parse_annot_header(tsv, input_custom_annotations_tsv)
            self.assertEqual(actual_header, expected_header)
            self.assertEqual(actual_column_names, ['AC', 'AN', 'AF', 'Hom', 'AC_eas', 'AN_eas', 'AF_eas', 'Hom_eas'])
            self.assertEqual(actual_names_to_column_map, {'#CHROM': 0, 'POS': 1, 'REF': 2, 'ALT': 3, 'AC':4, 'AN':5,
                                                          'AF':6, 'Hom':7, 'AC_eas': 8, 'AN_eas': 9, 'AF_eas': 10, 'Hom_eas': 11})





