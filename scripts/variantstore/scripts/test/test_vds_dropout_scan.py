#!/usr/bin/env python3
"""Unit tests for the Hail-free parts of vds_dropout_scan.

The module guards its ``import hail`` so everything below runs without Hail installed.
What is covered here is the logic that would be silently wrong rather than loudly
broken: superpartition arithmetic that has to line up with the `vet_NNN` table names,
sample selection that has to be reproducible across runs and VDSes, bin boundaries, and
the argument validation that stops a misconfigured run before it books a cluster.

The Hail aggregation itself is not covered; it is exercised by the probe action against
real data.
"""

import argparse
import contextlib
import io
import unittest

import vds_dropout_scan as vds


class TestSuperpartitionArithmetic(unittest.TestCase):
    """Must match CAST(CEIL(sample_id / 4000.0) AS INT64) in the Avro export WDL."""

    def test_boundaries(self):
        self.assertEqual(1, vds.superpartition_for(1))
        self.assertEqual(1, vds.superpartition_for(4000))
        self.assertEqual(2, vds.superpartition_for(4001))
        self.assertEqual(2, vds.superpartition_for(8000))

    def test_known_affected_superpartitions(self):
        """The two real Foxtrot dropouts, from the VS-1946 sample ID ranges."""
        self.assertEqual(83, vds.superpartition_for(328_001))
        self.assertEqual(83, vds.superpartition_for(332_000))
        self.assertEqual(64, vds.superpartition_for(252_001))
        self.assertEqual(64, vds.superpartition_for(256_000))
        # One past each range lands in the neighbouring superpartition.
        self.assertEqual(84, vds.superpartition_for(332_001))
        self.assertEqual(65, vds.superpartition_for(256_001))

    def test_foxtrot_scale(self):
        self.assertEqual(134, vds.superpartition_for(535_000))

    def test_custom_size(self):
        self.assertEqual(1, vds.superpartition_for(100, superpartition_size=100))
        self.assertEqual(2, vds.superpartition_for(101, superpartition_size=100))

    def test_non_positive_sample_id_raises(self):
        for bad in (0, -1):
            with self.assertRaises(ValueError):
                vds.superpartition_for(bad)


class TestBinArithmetic(unittest.TestCase):

    def test_bin_start_is_one_based_inclusive(self):
        self.assertEqual(1, vds.bin_start_for(1))
        self.assertEqual(1, vds.bin_start_for(50_000))
        self.assertEqual(50_001, vds.bin_start_for(50_001))

    def test_known_dropout_boundaries(self):
        """The real chr4 window's edges, to pin the quantization the report will show."""
        self.assertEqual(56_550_001, vds.bin_start_for(56_585_368))
        self.assertEqual(57_000_001, vds.bin_start_for(57_035_833))

    def test_bin_index_and_start_agree(self):
        for position in (1, 49_999, 50_000, 50_001, 56_585_368):
            index = vds.bin_index_for(position)
            self.assertEqual(vds.bin_start_for(position), index * 50_000 + 1)


class TestStride(unittest.TestCase):

    def test_stride_one_keeps_everything(self):
        self.assertTrue(all(vds.stride_keeps_bin(i, 1) for i in range(20)))

    def test_stride_five_keeps_one_in_five(self):
        kept = [i for i in range(20) if vds.stride_keeps_bin(i, 5)]
        self.assertEqual([0, 5, 10, 15], kept)

    def test_dropout_wider_than_stride_always_contains_a_scanned_bin(self):
        """The pigeonhole guarantee striding rests on."""
        stride, bin_size = 5, 50_000
        span_bins = (550_000 // bin_size) + 1  # the real chr19 window
        for offset in range(stride * 2):
            window = range(offset, offset + span_bins)
            self.assertTrue(any(vds.stride_keeps_bin(i, stride) for i in window))

    def test_invalid_stride_raises(self):
        with self.assertRaises(ValueError):
            vds.stride_keeps_bin(0, 0)


class TestSampleSelection(unittest.TestCase):

    def make_assignment(self, n_superpartitions=10, per_superpartition=4000):
        assignment = {}
        for sp in range(1, n_superpartitions + 1):
            for i in range(per_superpartition):
                assignment[f's{(sp - 1) * per_superpartition + i:07d}'] = sp
        return assignment

    def test_takes_requested_depth_from_each_superpartition(self):
        chosen = vds.stratified_sample(self.make_assignment(), 100)
        self.assertEqual(10, len(chosen))
        for sp, names in chosen.items():
            self.assertEqual(100, len(names), msg=f'superpartition {sp}')

    def test_selection_is_deterministic(self):
        assignment = self.make_assignment()
        self.assertEqual(
            vds.stratified_sample(assignment, 100),
            vds.stratified_sample(assignment, 100),
        )

    def test_selection_is_independent_of_input_ordering(self):
        """Reproducibility must not depend on dict iteration order."""
        assignment = self.make_assignment()
        reversed_assignment = dict(reversed(list(assignment.items())))
        self.assertEqual(
            vds.stratified_sample(assignment, 100),
            vds.stratified_sample(reversed_assignment, 100),
        )

    def test_different_seeds_select_differently(self):
        assignment = self.make_assignment()
        a = vds.stratified_sample(assignment, 100, seed='seed-a')
        b = vds.stratified_sample(assignment, 100, seed='seed-b')
        self.assertNotEqual(a, b)

    def test_selection_is_stable_when_unrelated_samples_are_added(self):
        """Adding a superpartition must not perturb the others' selections.

        This is what lets a later callset be screened on the same samples as an earlier
        one for the superpartitions they share.
        """
        base = self.make_assignment(n_superpartitions=5)
        extended = dict(base)
        for i in range(4000):
            extended[f'new{i:07d}'] = 6

        chosen_base = vds.stratified_sample(base, 100)
        chosen_extended = vds.stratified_sample(extended, 100)
        for sp in range(1, 6):
            self.assertEqual(chosen_base[sp], chosen_extended[sp], msg=f'superpartition {sp}')

    def test_short_superpartition_contributes_everything(self):
        """The final superpartition of a callset holds fewer than 4000 samples."""
        assignment = self.make_assignment(n_superpartitions=2)
        assignment.update({f'tail{i}': 3 for i in range(20)})
        chosen = vds.stratified_sample(assignment, 100)
        self.assertEqual(20, len(chosen[3]))

    def test_depth_exceeding_every_superpartition_keeps_all(self):
        assignment = self.make_assignment(n_superpartitions=3, per_superpartition=10)
        chosen = vds.stratified_sample(assignment, 100)
        self.assertEqual({1: 10, 2: 10, 3: 10}, {k: len(v) for k, v in chosen.items()})

    def test_invalid_depth_raises(self):
        with self.assertRaises(ValueError):
            vds.stratified_sample({'a': 1}, 0)

    def test_sort_key_is_stable_across_processes(self):
        """Hard-coded because a salted hash would break cross-run reproducibility."""
        self.assertEqual(
            vds.sample_sort_key('1234567', seed='vs-1998'),
            vds.sample_sort_key('1234567', seed='vs-1998'),
        )
        self.assertNotEqual(
            vds.sample_sort_key('1234567', seed='vs-1998'),
            vds.sample_sort_key('1234567', seed='other'),
        )


class TestSampleMapParsing(unittest.TestCase):

    def test_parses_tab_separated(self):
        text = 'sample_name\tsample_id\n5545879\t328001\n5648802\t328002\n'
        self.assertEqual({'5545879': 328001, '5648802': 328002},
                         vds.parse_sample_map(io.StringIO(text)))

    def test_parses_comma_separated(self):
        """bq query --format=csv output, if the caller forgets the tr."""
        text = 'sample_name,sample_id\n5545879,328001\n'
        self.assertEqual({'5545879': 328001}, vds.parse_sample_map(io.StringIO(text)))

    def test_blank_lines_are_skipped(self):
        text = 'sample_name\tsample_id\n5545879\t328001\n\n5648802\t328002\n'
        self.assertEqual(2, len(vds.parse_sample_map(io.StringIO(text))))

    def test_missing_header_raises(self):
        with self.assertRaises(ValueError):
            vds.parse_sample_map(io.StringIO('5545879\t328001\n'))

    def test_empty_input_raises(self):
        with self.assertRaises(ValueError):
            vds.parse_sample_map(io.StringIO('sample_name\tsample_id\n'))

    def test_non_integer_sample_id_raises(self):
        with self.assertRaises(ValueError):
            vds.parse_sample_map(io.StringIO('sample_name\tsample_id\nx\tnotanumber\n'))


class TestSuperpartitionAssignment(unittest.TestCase):

    def test_restricts_to_present_samples(self):
        sample_map = {'a': 1, 'b': 4001, 'c': 8001}
        assigned = vds.assign_superpartitions(sample_map, ['a', 'b'])
        self.assertEqual({'a': 1, 'b': 2}, assigned)

    def test_sample_missing_from_map_raises(self):
        """Silently dropping it would bias every superpartition total."""
        with self.assertRaises(ValueError) as ctx:
            vds.assign_superpartitions({'a': 1}, ['a', 'ghost'])
        self.assertIn('ghost', str(ctx.exception))


class TestContigParsing(unittest.TestCase):

    def test_default_is_all_primary_contigs(self):
        self.assertEqual(24, len(vds.parse_contig_list(None)))
        self.assertEqual(24, len(vds.parse_contig_list('')))

    def test_explicit_list(self):
        self.assertEqual(['chr4', 'chr19'], vds.parse_contig_list('chr4,chr19'))

    def test_whitespace_tolerated(self):
        self.assertEqual(['chr4', 'chr19'], vds.parse_contig_list(' chr4 , chr19 '))

    def test_non_gvs_contig_raises(self):
        """Alt and decoy contigs have no GVS location encoding."""
        with self.assertRaises(ValueError):
            vds.parse_contig_list('chr4,chr1_KI270706v1_random')

    def test_contig_list_includes_sex_chromosomes(self):
        contigs = vds.parse_contig_list(None)
        self.assertIn('chrX', contigs)
        self.assertIn('chrY', contigs)


class TestClusterRunnerCompatibility(unittest.TestCase):
    """run_in_hail_cluster.py renders every arguments-JSON key as `--key value`.

    That rules out `store_true` flags, which would receive an argument they cannot take,
    and `action='append'` used alone, since the JSON holds one value per key.
    """

    def test_boolean_flags_accept_explicit_values(self):
        for text in ('true', 'True', 'TRUE', 't', 'yes', '1'):
            self.assertTrue(vds.parse_bool(text), msg=text)
        for text in ('false', 'False', 'f', 'no', '0'):
            self.assertFalse(vds.parse_bool(text), msg=text)

    def test_boolean_passthrough(self):
        self.assertTrue(vds.parse_bool(True))
        self.assertFalse(vds.parse_bool(False))

    def test_bad_boolean_raises(self):
        with self.assertRaises(ValueError):
            vds.parse_bool('maybe')

    def test_overwrite_parses_from_a_value(self):
        parser = vds.build_parser()
        args = parser.parse_args([
            '--action', 'materialize', '--vds-path', 'gs://x', '--output-path', 'gs://o',
            '--superpartitions-path', 'sp.tsv', '--sample-list-path', 's.tsv',
            '--sample-map-path', 'm.tsv', '--overwrite', 'true',
        ])
        self.assertIs(True, args.overwrite)

    def test_no_store_true_flags_remain(self):
        parser = vds.build_parser()
        for action in parser._actions:
            self.assertNotEqual(
                'store_true', type(action).__name__.replace('_Store', 'store').lower(),
                msg=str(action.option_strings),
            )
            self.assertFalse(
                isinstance(action, argparse._StoreTrueAction),
                msg=f'{action.option_strings} is store_true and cannot survive the '
                    'arguments-JSON round trip',
            )

    def test_intervals_accept_a_single_comma_separated_string(self):
        self.assertEqual(
            ['chr4:56000000-58000000', 'chr19:39000000-41000000'],
            vds.parse_interval_list(['chr4:56000000-58000000,chr19:39000000-41000000']),
        )

    def test_intervals_accept_repetition(self):
        self.assertEqual(
            ['chr4:1-2', 'chr19:3-4'],
            vds.parse_interval_list(['chr4:1-2', 'chr19:3-4']),
        )

    def test_intervals_accept_a_mixture(self):
        self.assertEqual(
            ['a', 'b', 'c'], vds.parse_interval_list(['a,b', 'c']),
        )

    def test_empty_intervals(self):
        self.assertEqual([], vds.parse_interval_list(None))
        self.assertEqual([], vds.parse_interval_list([]))
        self.assertEqual([], vds.parse_interval_list(['  ,  ']))


class TestOutputFormatting(unittest.TestCase):

    def test_summary_rows_skip_zero_cells(self):
        totals = {('chr4', 1): {83: 0.0, 84: 7500.0}}
        rows = vds.format_summary_rows(totals, 50_000)
        self.assertEqual(['chr4\t1\t50001\t84\t7500'], rows)

    def test_summary_rows_are_ordered_by_contig_then_position(self):
        totals = {
            ('chr19', 50_001): {1: 10.0},
            ('chr4', 100_001): {1: 10.0},
            ('chr4', 1): {1: 10.0},
        }
        rows = vds.format_summary_rows(totals, 50_000)
        self.assertEqual(
            ['chr4\t1\t50001\t1\t10', 'chr4\t100001\t150001\t1\t10',
             'chr19\t50001\t100001\t1\t10'],
            rows,
        )

    def test_summary_row_end_is_exclusive(self):
        rows = vds.format_summary_rows({('chr1', 1): {1: 5.0}}, 50_000)
        self.assertTrue(rows[0].endswith('\t1\t50001\t1\t5'))

    def test_summary_header_matches_detector_expectation(self):
        """The two scripts are only useful together, so the contract is asserted here."""
        import vds_dropout_detect as vdd
        self.assertEqual('\t'.join(vdd.SUMMARY_COLUMNS), vds.SUMMARY_HEADER)
        self.assertEqual('\t'.join(vdd.SUPERPARTITION_COLUMNS), vds.SUPERPARTITION_HEADER)

    def test_superpartition_rows_count_samples(self):
        chosen = {83: ['a', 'b'], 1: ['c']}
        self.assertEqual(['1\t1', '83\t2'], vds.format_superpartition_rows(chosen))

    def test_sample_list_round_trips(self):
        chosen = {83: ['5545879'], 64: ['3074794']}
        sample_map = {'5545879': 328_001, '3074794': 252_001}
        rows = vds.format_sample_list_rows(chosen, sample_map)
        self.assertEqual(['3074794\t252001\t64', '5545879\t328001\t83'], rows)


class TestProbeExtrapolation(unittest.TestCase):

    def test_scales_linearly(self):
        self.assertAlmostEqual(
            100.0, vds.extrapolate_runtime(1.0, vds.GENOME_LENGTH // 100, vds.GENOME_LENGTH),
            places=0,
        )

    def test_zero_span_raises(self):
        with self.assertRaises(ValueError):
            vds.extrapolate_runtime(1.0, 0)

    def test_report_mentions_shuffle_check(self):
        report = vds.format_probe_report(60.0, 10_000_000, 1_340, 200)
        self.assertIn('shuffle', report)
        self.assertIn('genome-wide estimate', report)


class TestArgumentValidation(unittest.TestCase):
    """Catch misconfiguration before it books a cluster."""

    def parse(self, *argv):
        parser = vds.build_parser()
        args = parser.parse_args(argv)
        vds.validate_args(args, parser)
        return args

    def assert_rejected(self, *argv):
        """Assert these arguments are rejected, quietly.

        argparse prints its whole usage block to stderr on every rejection, so asserting
        SystemExit directly makes a passing run look like a wall of failures during the
        Docker image build. The message is still available here if a test needs it.
        """
        captured = io.StringIO()
        with contextlib.redirect_stderr(captured):
            with self.assertRaises(SystemExit):
                self.parse(*argv)
        return captured.getvalue()

    def test_materialize_requires_outputs(self):
        self.assert_rejected('--action', 'materialize', '--vds-path', 'gs://x',
                             '--sample-map-path', 'm.tsv')

    def test_materialize_accepts_full_argument_set(self):
        args = self.parse(
            '--action', 'materialize', '--vds-path', 'gs://in.vds',
            '--sample-map-path', 'm.tsv', '--output-path', 'gs://out.vds',
            '--superpartitions-path', 'sp.tsv', '--sample-list-path', 'samples.tsv',
        )
        self.assertEqual('materialize', args.action)
        self.assertEqual(vds.DEFAULT_SAMPLES_PER_SUPERPARTITION,
                         args.samples_per_superpartition)

    def test_materialize_requires_sample_map_not_just_list(self):
        self.assert_rejected('--action', 'materialize', '--vds-path', 'gs://x',
                             '--output-path', 'gs://o', '--superpartitions-path', 'sp.tsv',
                             '--sample-list-path', 'samples.tsv')

    def test_scan_requires_a_sample_source(self):
        self.assert_rejected('--action', 'scan', '--vds-path', 'gs://x',
                             '--summary-path', 's.tsv', '--superpartitions-path', 'sp.tsv')

    def test_scan_accepts_recorded_sample_list(self):
        args = self.parse('--action', 'scan', '--vds-path', 'gs://x',
                          '--summary-path', 's.tsv', '--superpartitions-path', 'sp.tsv',
                          '--sample-list-path', 'samples.tsv')
        self.assertEqual('samples.tsv', args.sample_list_path)

    def test_full_depth_requires_output(self):
        self.assert_rejected('--action', 'full-depth', '--vds-path', 'gs://x',
                             '--sample-list-path', 'samples.tsv')

    def test_bad_stride_rejected(self):
        message = self.assert_rejected('--action', 'probe', '--vds-path', 'gs://x',
                                       '--sample-list-path', 's.tsv', '--stride', '0')
        self.assertIn('--stride must be at least 1', message)

    def test_bad_bin_size_rejected(self):
        message = self.assert_rejected('--action', 'probe', '--vds-path', 'gs://x',
                                       '--sample-list-path', 's.tsv', '--bin-size', '0')
        self.assertIn('--bin-size must be at least 1', message)

    def test_intervals_are_repeatable(self):
        args = self.parse('--action', 'scan', '--vds-path', 'gs://x',
                          '--summary-path', 's.tsv', '--superpartitions-path', 'sp.tsv',
                          '--sample-list-path', 'samples.tsv',
                          '--intervals', 'chr4:56000000-58000000',
                          '--intervals', 'chr19:39000000-41000000')
        self.assertEqual(2, len(args.intervals))

    def test_unknown_action_rejected(self):
        self.assert_rejected('--action', 'nonsense', '--vds-path', 'gs://x')

    def test_flags_are_kebab_case(self):
        """run_in_hail_cluster.py turns arguments-JSON keys into --key value pairs."""
        parser = vds.build_parser()
        for action in parser._actions:
            for option in action.option_strings:
                self.assertNotIn('_', option, msg=option)


class TestHailGuard(unittest.TestCase):

    def test_hail_backed_helpers_fail_clearly_without_hail(self):
        """The module must import without Hail so these tests can run at all."""
        if vds.hl is not None:  # pragma: no cover - only on a cluster
            self.skipTest('hail is installed')
        with self.assertRaises(RuntimeError) as ctx:
            vds._require_hail()
        self.assertIn('Dataproc', str(ctx.exception))


if __name__ == '__main__':
    unittest.main()
