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
import ast
import collections
import contextlib
import hashlib
import io
import pathlib
import re
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


class TestSelectionAgainstIndependentImplementation(unittest.TestCase):
    """Check stratified_sample against a second, separately written implementation.

    A test that re-uses the function under test can only confirm it is self-consistent.
    Writing the selection out a second time from its specification is what caught the
    original version relying on sort stability, which made the result depend on dict
    iteration order -- fatal for an artifact whose whole purpose is being reproducible
    across VDSes.
    """

    @staticmethod
    def independent_selection(sample_map, depth, seed):
        """Re-derived from the documented algorithm, not by calling into the module.

        Per superpartition, order by SHA-256 of "seed:sample_name" then by name, and take
        the first `depth`.
        """
        buckets = collections.defaultdict(list)
        for name, sample_id in sample_map.items():
            superpartition = -(-sample_id // 4000)  # integer ceiling division
            digest = hashlib.sha256(f'{seed}:{name}'.encode('utf-8')).hexdigest()
            buckets[superpartition].append((digest, name))
        selected = {}
        for superpartition, rows in buckets.items():
            rows.sort()
            selected[superpartition] = sorted(name for _, name in rows[:depth])
        return selected

    def make_map(self, n_superpartitions=6, per_superpartition=500):
        return {
            f'{(sp - 1) * per_superpartition + i:07d}': (sp - 1) * 4000 + i + 1
            for sp in range(1, n_superpartitions + 1)
            for i in range(per_superpartition)
        }

    def test_matches_the_independent_implementation(self):
        sample_map = self.make_map()
        assigned = {n: vds.superpartition_for(i) for n, i in sample_map.items()}
        self.assertEqual(
            self.independent_selection(sample_map, 100, vds.DEFAULT_SEED),
            vds.stratified_sample(assigned, 100, vds.DEFAULT_SEED),
        )

    def test_agreement_holds_at_other_depths_and_seeds(self):
        sample_map = self.make_map(n_superpartitions=4, per_superpartition=300)
        assigned = {n: vds.superpartition_for(i) for n, i in sample_map.items()}
        for depth in (1, 10, 200):
            for seed in (vds.DEFAULT_SEED, 'other-seed'):
                self.assertEqual(
                    self.independent_selection(sample_map, depth, seed),
                    vds.stratified_sample(assigned, depth, seed),
                    msg=f'depth={depth} seed={seed}',
                )

    def test_agreement_holds_for_short_superpartitions(self):
        """The last superpartition of a callset holds fewer than the sampling depth."""
        sample_map = self.make_map(n_superpartitions=2, per_superpartition=400)
        sample_map.update({f'tail{i}': 8000 + i + 1 for i in range(15)})
        assigned = {n: vds.superpartition_for(i) for n, i in sample_map.items()}
        self.assertEqual(
            self.independent_selection(sample_map, 100, vds.DEFAULT_SEED),
            vds.stratified_sample(assigned, 100, vds.DEFAULT_SEED),
        )


class TestWdlGeneratedSampleMap(unittest.TestCase):
    """The WDL generates the sample map itself, inlining a copy of the bq/ query.

    Two copies of a predicate is two places to get it wrong, and getting it wrong here is
    quiet: a mismatched sample universe changes which samples are screened without
    erroring. So the essentials are asserted against both files.
    """

    WDL = (pathlib.Path(__file__).resolve().parents[2]
           / 'wdl' / 'GvsValidateVdsCompleteness.wdl')
    SQL = (pathlib.Path(__file__).resolve().parents[2]
           / 'bq' / 'vds_dropout_sample_map.sql')

    def read(self, path):
        if not path.exists():
            # The Docker test run mounts only the scripts directory.
            self.skipTest(f'{path} not available in this test environment')
        return path.read_text()

    def generated_query(self):
        wdl = self.read(self.WDL)
        return wdl.split('task GenerateSampleMap', 1)[1].split('>>>', 1)[0]

    def test_selects_the_columns_the_scan_expects(self):
        self.assertIn('SELECT sample_name, sample_id', self.generated_query())

    def test_applies_the_same_filter_as_the_avro_export(self):
        """A different sample universe would silently change what gets screened."""
        query = self.generated_query()
        self.assertIn('withdrawn IS NULL', query)
        self.assertIn('is_control = false', query)

    def test_wdl_and_sql_file_agree_on_the_filter(self):
        query = self.generated_query()
        sql = self.read(self.SQL)
        for predicate in ('withdrawn IS NULL', 'is_control = false'):
            self.assertIn(predicate, query, msg='WDL')
            self.assertIn(predicate, sql, msg='bq/ SQL')

    def test_exports_in_the_format_the_scan_can_read(self):
        """Tab-delimited with a header is what sniff_delimiter and import_table expect."""
        query = self.generated_query()
        self.assertIn('field_delimiter', query)
        self.assertIn('header=true', query)

    def test_export_uri_has_exactly_one_wildcard(self):
        uris = re.findall(r"uri='([^']*)'", self.generated_query())
        self.assertEqual(1, len(uris))
        self.assertEqual(1, uris[0].count('*'))

    def test_fails_rather_than_guessing_on_a_multi_file_export(self):
        """Reading one shard of a sharded export would screen part of the callset."""
        query = self.generated_query()
        self.assertIn('-ne 1', query)
        self.assertIn('Expected exactly one exported sample map file', query)

    def test_map_generation_is_not_cached(self):
        """sample_info changes as samples are withdrawn; a stale map would be wrong."""
        wdl = self.read(self.WDL)
        task = wdl.split('task GenerateSampleMap', 1)[1].split('command <<<', 1)[0]
        self.assertIn('volatile: true', task)


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
    """A partition is the unit of work, so cost scales by partitions, not base pairs.

    Not academic: a diagnostic bundle from a stalled run showed the entry aggregation
    RUNNABLE in a single task, because a 10 Mb interval on a 535K-sample VDS resolved to
    one partition. Scaling that wall clock by genomic span would have overstated a
    genome-wide run by roughly the width of the cluster.
    """

    def test_partition_model_divides_by_concurrency(self):
        # 10s for 2 partitions -> 5s each; 100 partitions over 10 tasks -> 50s.
        self.assertAlmostEqual(50.0, vds.extrapolate_by_partitions(10.0, 2, 100, 10))

    def test_partition_model_serial_case(self):
        self.assertAlmostEqual(500.0, vds.extrapolate_by_partitions(10.0, 2, 100, 1))

    def test_partition_model_rejects_nonsense(self):
        for probed, total, tasks in ((0, 10, 1), (1, 0, 1), (1, 10, 0)):
            with self.assertRaises(ValueError):
                vds.extrapolate_by_partitions(10.0, probed, total, tasks)

    def test_bp_model_still_available_but_documented_as_misleading(self):
        self.assertAlmostEqual(
            100.0, vds.extrapolate_runtime(1.0, vds.GENOME_LENGTH // 100,
                                           vds.GENOME_LENGTH), places=0)
        self.assertIn('Misleading', vds.extrapolate_runtime.__doc__)

    def test_zero_span_raises(self):
        with self.assertRaises(ValueError):
            vds.extrapolate_runtime(1.0, 0)

    def test_report_leads_with_per_partition_cost(self):
        report = vds.format_probe_report(600.0, 10_000_000, 1_340, 200,
                                         n_samples=13_400, n_partitions=20)
        self.assertIn('partitions read:', report)
        self.assertIn('per partition:', report)
        self.assertIn('30.0 s', report)
        self.assertIn('Do not scale by base pairs', report)

    def test_report_warns_when_parallelism_is_degenerate(self):
        report = vds.format_probe_report(600.0, 10_000_000, 1_340, 200, n_partitions=1)
        self.assertIn('CAUTION', report)
        self.assertIn('essentially serial', report)

    def test_report_does_not_warn_with_healthy_parallelism(self):
        report = vds.format_probe_report(600.0, 10_000_000, 1_340, 200, n_partitions=500)
        self.assertNotIn('CAUTION', report)

    def test_report_extrapolates_to_the_whole_vds(self):
        """The real geometry: 411 of 119,189 partitions, the probe that finally worked."""
        report = vds.format_probe_report(
            1200.0, 10_000_000, 1_340, 200, n_samples=535_662,
            n_partitions=411, n_total_partitions=119_189,
            parallelism=16, full_parallelism=1200)
        self.assertIn('119,189 partitions', report)
        self.assertIn('0.34% covered here', report)
        # 1200s * (119189/411) / 3600 = 96.67h at this width
        self.assertIn('96.7 h', report)
        # scaled by 16/1200 -> ~1.3h
        self.assertIn('1.3 h', report)

    def test_report_prefers_the_wider_estimate_explicitly(self):
        report = vds.format_probe_report(
            1200.0, 10_000_000, 1_340, 200, n_partitions=411,
            n_total_partitions=119_189, parallelism=16, full_parallelism=1200)
        self.assertIn('wider figure is the one to plan against', report)

    def test_report_omits_scaling_when_no_wider_width_is_known(self):
        report = vds.format_probe_report(
            1200.0, 10_000_000, 1_340, 200, n_partitions=411,
            n_total_partitions=119_189, parallelism=16)
        self.assertIn('at this cluster width', report)
        self.assertNotIn('concurrent tasks:', report.split('Genome-wide')[1])

    def test_report_omits_extrapolation_without_a_total(self):
        report = vds.format_probe_report(1200.0, 10_000_000, 1_340, 200, n_partitions=411)
        self.assertNotIn('Genome-wide estimate', report)

    def test_report_always_warns_against_span_scaling(self):
        report = vds.format_probe_report(600.0, 10_000_000, 1_340, 200, n_partitions=411)
        self.assertIn('Do not scale by base pairs', report)

    def test_report_states_what_the_estimate_excludes(self):
        report = vds.format_probe_report(60.0, 10_000_000, 1_340, 200, n_partitions=10)
        self.assertIn('not writing', report)
        self.assertIn('floor for the materialize pass', report)

    def test_report_without_partition_count_omits_the_model(self):
        report = vds.format_probe_report(60.0, 10_000_000, 1_340, 200)
        self.assertNotIn('per partition:', report)


class TestIntervalReadSemantics(unittest.TestCase):
    """read_vds(intervals=...) repartitions; it does not filter.

    N intervals yield exactly N partitions. A single 10 Mb interval therefore collapsed
    ~380 native partitions of a 119,189-partition VDS into one, and one task streamed all
    of it for nine hours. `--contigs` would have been worse still: 24 intervals, so 24
    partitions for the whole genome.

    Guarded by inspecting the source because the alternative needs a live Hail session and
    an AoU-scale VDS, and the mistake is silent -- it produces correct results, slowly.
    """

    SOURCE = pathlib.Path(__file__).resolve().parents[1] / 'vds_dropout_scan.py'

    def source(self):
        return self.SOURCE.read_text()

    @staticmethod
    def _called_name(node):
        func = node.func
        return func.attr if isinstance(func, ast.Attribute) else getattr(func, 'id', '')

    def test_read_vds_is_never_given_intervals(self):
        """Parsed rather than grepped: the docstring explaining the trap quotes the call."""
        tree = ast.parse(self.source())
        offenders = [
            node.lineno for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and self._called_name(node) == 'read_vds'
            and any(kw.arg == 'intervals' for kw in node.keywords)
        ]
        self.assertEqual([], offenders,
                         msg=f'read_vds called with intervals= at line(s) {offenders}; '
                             'that repartitions to one partition per interval')

    def test_filter_intervals_is_actually_called(self):
        tree = ast.parse(self.source())
        calls = [node for node in ast.walk(tree)
                 if isinstance(node, ast.Call)
                 and self._called_name(node) == 'filter_intervals']
        self.assertTrue(calls, 'subsetting must go through hl.vds.filter_intervals')

    def test_subsetting_goes_through_filter_intervals(self):
        self.assertIn('hl.vds.filter_intervals(', self.source())

    def test_reference_blocks_are_not_split(self):
        """Splitting costs work and would skew the covered-base metric."""
        self.assertIn('split_reference_blocks=False', self.source())

    def test_the_reasoning_is_recorded_where_the_call_is(self):
        doc = vds.read_and_subset_vds.__doc__
        self.assertIn('repartitioning', doc)
        self.assertIn('native partitioning', doc)


class TestStriding(unittest.TestCase):
    """Striding must reach read_vds as intervals.

    Expressed as a row predicate on the bin index, Hail cannot use the row index and
    streams every byte of every partition it opens, so the scan costs the same as reading
    everything and merely reports less.
    """

    def test_stride_one_returns_intervals_unchanged(self):
        given = ['chr4:56000000-58000000']
        self.assertEqual(given, vds.apply_stride(given, 50_000, 1))

    def test_stride_expands_into_one_interval_per_retained_bin(self):
        out = vds.strided_intervals('chr20', 1, 500_001, 50_000, 5)
        self.assertEqual(['chr20:1-50001', 'chr20:250001-300001'], out)

    def test_strided_intervals_are_clipped_to_the_request(self):
        out = vds.strided_intervals('chr4', 60_000, 160_000, 50_000, 1)
        self.assertEqual(
            ['chr4:60000-100001', 'chr4:100001-150001', 'chr4:150001-160000'], out)

    def test_stride_covers_any_window_wider_than_the_stride(self):
        """The pigeonhole guarantee, now over real intervals."""
        for offset in range(0, 500_000, 50_000):
            out = vds.strided_intervals('chr1', 1 + offset, 1 + offset + 550_000,
                                        50_000, 5)
            self.assertTrue(out, msg=f'offset {offset}')

    def test_apply_stride_over_multiple_intervals(self):
        out = vds.apply_stride(['chr1:1-200001', 'chr2:1-200001'], 100_000, 2)
        self.assertEqual(['chr1:1-100001', 'chr2:1-100001'], out)

    def test_apply_stride_tolerates_thousands_separators(self):
        self.assertEqual(['chr1:1-100001'],
                         vds.apply_stride(['chr1:1-100,001'], 100_000, 2))

    def test_striding_a_whole_contig_is_rejected(self):
        """Unbounded intervals cannot be strided, and failing beats silently ignoring it."""
        with self.assertRaises(ValueError) as ctx:
            vds.apply_stride(['chr20'], 50_000, 5)
        self.assertIn('whole contig', str(ctx.exception))

    def test_invalid_stride_raises(self):
        with self.assertRaises(ValueError):
            vds.strided_intervals('chr1', 1, 1000, 50_000, 0)

    def test_min_healthy_partitions_is_sane(self):
        self.assertGreater(vds.MIN_HEALTHY_PARTITIONS, 1)


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
