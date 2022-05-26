import unittest

from scale_xy_bed_values import scale_xy_bed_values


class TestScaleXYBedValues(unittest.TestCase):
    def test_downsampled_scaling(self):
        import filecmp
        import tempfile

        for x, y in [(10, 10), (1, 10), (10, 1)]:
            with tempfile.NamedTemporaryFile() as actual_output_bed:
                scale_xy_bed_values('scale_xy_bed_values_test_files/intervals_downsampled_5.bed',
                                    actual_output_bed.name,
                                    x,
                                    y)

                expected_output_bed = f'scale_xy_bed_values_test_files/intervals_downsampled_5_scaled_{x}_{y}.bed'
                self.assertTrue(filecmp.cmp(expected_output_bed,
                                            actual_output_bed.name,
                                            shallow=False), f'fail on X scaling {x} and Y scaling {y}')
