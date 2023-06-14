import unittest

from filter_interval_list_chromosomes import filter_chromosomes


class TestFilterIntervalListChromosomes(unittest.TestCase):
    def test_filter_interval_list(self):
        import filecmp
        import tempfile
        test_dir = 'filter_interval_list_chromosomes'
        full_interval_list_path = 'wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list'

        with tempfile.NamedTemporaryFile() as actual_output_file:
            chromosomes = ['chr20', 'chrX', 'chrY']
            filter_chromosomes(actual_output_file.name, f'{test_dir}/{full_interval_list_path}', *chromosomes)

            expected_filtered_chromosomes = f'{test_dir}/chrX_chrY_chr20_filtered.interval_list'
            self.assertTrue(filecmp.cmp(actual_output_file.name,
                                        expected_filtered_chromosomes,
                                        shallow=False), f'fail on filtering {chromosomes}')
