from os.path import abspath
import unittest

from scale_xy_bed_values import scale_xy_bed_values


class TestScaleXYBedValues(unittest.TestCase):
    def test_downsampled_scaling(self):
        import filecmp
        import tempfile

        for x, y in [(10, 10), (1, 10), (10, 1)]:
            with tempfile.NamedTemporaryFile() as tmp:
                scale_xy_bed_values('scale_xy_bed_values_test_files/intervals_downsampled_5.bed',
                                    tmp.name,
                                    x,
                                    y)
                self.assertTrue(filecmp.cmp(f'scale_xy_bed_values_test_files/intervals_downsampled_5_scaled_{x}_{y}.bed',
                                            tmp.name,
                                            shallow=False), f'fail on X scaling {x} and Y scaling {y}')
