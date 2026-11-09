#!/usr/bin/env python3
"""Unit tests for the sample-set check and the left-alignment classifier in
hail_create_participant_mapping_inputs.

Hail is not installed in the Variants image or in CI, so a test that skipped without it would never run. The module
under test is therefore imported against a stub hail. The stub's hadoop_open is the builtin open, and its len, if_else
and case evaluate eagerly on Python strings, so key_class runs exactly the expression it builds, on plain values. That
covers the classification logic; Hail's typing of the expression is exercised only by a real run.
"""

import os
import sys
import tempfile
import types
import unittest


class _Case:
    """hl.case() evaluated eagerly: the first true condition's value, else the default."""

    def __init__(self):
        self._matched = False
        self._value = None

    def when(self, condition, value):
        if not self._matched and condition:
            self._matched, self._value = True, value
        return self

    def default(self, value):
        return self._value if self._matched else value


_hail = types.ModuleType('hail')
_hail.hadoop_open = open
_hail.len = len
_hail.if_else = lambda condition, if_true, if_false: if_true if condition else if_false
_hail.case = _Case
# Anything else, such as hl.MatrixTable or hl.expr.CallExpression, appears only in annotations, which
# hail_create_vat_inputs evaluates at import.
_hail.expr = types.ModuleType('hail.expr')
_hail.expr.__getattr__ = lambda name: object
_hail.__getattr__ = lambda name: object
sys.modules['hail'] = _hail

import hail_create_participant_mapping_inputs as mapping  # noqa: E402


class TestKeyClass(unittest.TestCase):

    def assertClass(self, ref, alt, expected):
        self.assertEqual(mapping.key_class(ref, alt), expected, f'{ref}-{alt}')

    def test_snv(self):
        self.assertClass('A', 'T', 'snv')

    def test_insertions(self):
        # The module docstring's examples: GG can move left onto the preceding G, AG cannot move past A.
        self.assertClass('G', 'GG', 'insertion_shiftable')
        self.assertClass('A', 'AG', 'insertion')
        # Tandem repeats: an inserted AC after C is also an inserted CA one base left, and that one step is enough.
        self.assertClass('C', 'CAC', 'insertion_shiftable')
        self.assertClass('A', 'ACAC', 'insertion')

    def test_deletions(self):
        self.assertClass('GG', 'G', 'deletion_shiftable')
        self.assertClass('GA', 'G', 'deletion')
        self.assertClass('CAC', 'C', 'deletion_shiftable')
        self.assertClass('CACA', 'C', 'deletion')

    def test_keys_the_anchor_test_does_not_cover(self):
        # Split keys are minimized, so neither should occur; write_keys_to_normalize fails the job if one does.
        self.assertClass('AT', 'GC', 'mnp_or_complex')
        self.assertClass('A', 'TT', 'unanchored')
        self.assertClass('TT', 'A', 'unanchored')
        self.assertNotIn('mnp_or_complex', mapping.SUPPORTED_CLASSES)
        self.assertNotIn('unanchored', mapping.SUPPORTED_CLASSES)

    def test_shiftable_classes_are_supported(self):
        self.assertTrue(set(mapping.SHIFTABLE_CLASSES) <= set(mapping.SUPPORTED_CLASSES))


class TestCheckSampleSet(unittest.TestCase):

    # name, status, is_loaded, as GetSampleInfo exports sample_info.
    SAMPLE_INFO = [
        ('1001', 'active', 'true'),
        ('1002', 'active', 'true'),
        ('1003', 'withdrawn', 'true'),
        ('1004', 'control', 'true'),
    ]

    def setUp(self):
        self.summary = {}

    def check(self, vds_samples, rows=None):
        with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as f:
            for row in self.SAMPLE_INFO if rows is None else rows:
                f.write('\t'.join(row) + '\n')
            # Blank lines are skipped.
            f.write('\n')
        self.addCleanup(os.remove, f.name)
        mapping.check_sample_set(set(vds_samples), f.name, self.summary)

    def assertProblem(self, vds_samples, problem, rows=None):
        with self.assertRaises(ValueError) as raised:
            self.check(vds_samples, rows)
        self.assertIn(problem, str(raised.exception))

    def test_exactly_the_active_samples_passes(self):
        self.check({'1001', '1002'})
        self.assertEqual(self.summary['sample_info'], {'active': 2, 'withdrawn': 1, 'control': 1})

    def test_withdrawn_sample_in_the_vds(self):
        self.assertProblem({'1001', '1002', '1003'}, '1 withdrawn in the VDS')

    def test_control_in_the_vds(self):
        self.assertProblem({'1001', '1002', '1004'}, '1 controls in the VDS')

    def test_sample_not_in_sample_info(self):
        self.assertProblem({'1001', '1002', '9999'}, '1 in the VDS but not in sample_info')

    def test_loaded_active_sample_missing(self):
        self.assertProblem({'1001'}, '1 active, is_loaded, missing from the VDS')

    def test_unloaded_active_sample_missing(self):
        # is_loaded only explains the absence; the sample is still expected.
        rows = self.SAMPLE_INFO + [('1005', 'active', 'false')]
        self.assertProblem({'1001', '1002'}, '1 active, not is_loaded, missing from the VDS', rows)

    def test_every_problem_is_reported(self):
        with self.assertRaises(ValueError) as raised:
            self.check({'1001', '1003', '1004', '9999'})
        message = str(raised.exception)
        for problem in ['withdrawn in the VDS', 'controls in the VDS', 'in the VDS but not in sample_info',
                        'active, is_loaded, missing from the VDS']:
            self.assertIn(problem, message)


if __name__ == '__main__':
    unittest.main()
