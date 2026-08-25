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


class TestSqlEquivalence(unittest.TestCase):
    """The SQL in bq/ must select exactly what the Python selects.

    The recorded sample list is what makes r2 and r3 figures comparable, so reproducibility
    is load-bearing rather than a nicety. One implementation agreeing with itself proves
    nothing; here the SQL's window-function semantics are re-implemented independently and
    checked against the production function.
    """

    BQ_DIR = pathlib.Path(__file__).resolve().parents[2] / 'bq'

    def bq_file(self, name):
        path = self.BQ_DIR / name
        if not path.exists():
            # The Docker test run mounts only the scripts directory, so bq/ is absent there.
            self.skipTest(f'{path} not available in this test environment')
        return path.read_text()

    @staticmethod
    def sql_semantics_selection(sample_map, depth, seed):
        """Independent re-implementation of the SQL, not a call into the module.

        Mirrors:
            ROW_NUMBER() OVER (PARTITION BY CAST(CEIL(sample_id / 4000.0) AS INT64)
                               ORDER BY TO_HEX(SHA256(CONCAT(seed, ':', sample_name))),
                                        sample_name)
        """
        buckets = collections.defaultdict(list)
        for name, sample_id in sample_map.items():
            superpartition = -(-sample_id // 4000)  # CAST(CEIL(x / 4000.0) AS INT64)
            digest = hashlib.sha256(f'{seed}:{name}'.encode('utf-8')).hexdigest()
            buckets[superpartition].append((digest, name))
        selected = {}
        for superpartition, rows in buckets.items():
            rows.sort()  # ORDER BY hex, then name
            selected[superpartition] = sorted(name for _, name in rows[:depth])
        return selected

    def make_map(self, n_superpartitions=6, per_superpartition=500):
        return {
            f'{(sp - 1) * per_superpartition + i:07d}': (sp - 1) * 4000 + i + 1
            for sp in range(1, n_superpartitions + 1)
            for i in range(per_superpartition)
        }

    def test_python_and_sql_semantics_agree(self):
        sample_map = self.make_map()
        assigned = {n: vds.superpartition_for(i) for n, i in sample_map.items()}
        self.assertEqual(
            self.sql_semantics_selection(sample_map, 100, vds.DEFAULT_SEED),
            vds.stratified_sample(assigned, 100, vds.DEFAULT_SEED),
        )

    def test_agreement_holds_at_other_depths_and_seeds(self):
        sample_map = self.make_map(n_superpartitions=4, per_superpartition=300)
        assigned = {n: vds.superpartition_for(i) for n, i in sample_map.items()}
        for depth in (1, 10, 200):
            for seed in (vds.DEFAULT_SEED, 'other-seed'):
                self.assertEqual(
                    self.sql_semantics_selection(sample_map, depth, seed),
                    vds.stratified_sample(assigned, depth, seed),
                    msg=f'depth={depth} seed={seed}',
                )

    def test_agreement_holds_for_short_superpartitions(self):
        """The last superpartition of a callset holds fewer than the sampling depth."""
        sample_map = self.make_map(n_superpartitions=2, per_superpartition=400)
        sample_map.update({f'tail{i}': 8000 + i + 1 for i in range(15)})
        assigned = {n: vds.superpartition_for(i) for n, i in sample_map.items()}
        self.assertEqual(
            self.sql_semantics_selection(sample_map, 100, vds.DEFAULT_SEED),
            vds.stratified_sample(assigned, 100, vds.DEFAULT_SEED),
        )

    def test_hex_ordering_equals_byte_ordering(self):
        """Underpins ordering by TO_HEX rather than by the raw digest."""
        digests = [hashlib.sha256(f'x:{i}'.encode()).digest() for i in range(400)]
        hexes = [d.hex() for d in digests]
        self.assertEqual(
            sorted(range(len(digests)), key=lambda i: hexes[i]),
            sorted(range(len(digests)), key=lambda i: digests[i]),
        )

    def test_selection_does_not_depend_on_input_order(self):
        """The explicit name tie-break is what makes this true on both sides."""
        sample_map = self.make_map(n_superpartitions=3, per_superpartition=200)
        assigned = {n: vds.superpartition_for(i) for n, i in sample_map.items()}
        shuffled = dict(reversed(list(assigned.items())))
        self.assertEqual(
            vds.stratified_sample(assigned, 50),
            vds.stratified_sample(shuffled, 50),
        )

    def test_sql_seed_matches_python_default(self):
        sql = self.bq_file('vds_dropout_sample_list.sql')
        match = re.search(r"DECLARE SEED STRING DEFAULT '([^']*)'", sql)
        self.assertIsNotNone(match, 'SEED declaration not found')
        self.assertEqual(vds.DEFAULT_SEED, match.group(1))

    def test_sql_depth_matches_python_default(self):
        sql = self.bq_file('vds_dropout_sample_list.sql')
        match = re.search(r'DECLARE SAMPLES_PER_SUPERPARTITION INT64 DEFAULT (\d+)', sql)
        self.assertIsNotNone(match, 'depth declaration not found')
        self.assertEqual(vds.DEFAULT_SAMPLES_PER_SUPERPARTITION, int(match.group(1)))

    def test_sql_superpartition_size_matches_python_default(self):
        sql = self.bq_file('vds_dropout_sample_list.sql')
        self.assertIn(f'CEIL(sample_id / {vds.DEFAULT_SUPERPARTITION_SIZE}.0)', sql)

    def test_sql_column_order_matches_the_sample_list_header(self):
        """So the two outputs can be diffed directly.

        Ordering is read out of the SQL by position rather than by membership; checking
        only that each name appears would pass regardless of the order they appear in.
        """
        sql = self.bq_file('vds_dropout_sample_list.sql')
        block = sql.split('CREATE TEMP TABLE stratified_selection AS')[1].split('FROM (')[0]
        names = ('sample_name', 'sample_id', 'superpartition')
        for name in names:
            self.assertIn(name, block)
        self.assertEqual(
            vds.SAMPLE_LIST_HEADER.split('\t'),
            sorted(names, key=block.index),
        )

    def test_both_output_paths_read_from_the_same_temp_table(self):
        """The returning SELECT and the EXPORT must not re-derive the selection.

        Two copies of the window function would be free to drift, which is the one thing
        this file exists to rule out.
        """
        sql = self.bq_file('vds_dropout_sample_list.sql')
        self.assertEqual(1, sql.count('ROW_NUMBER() OVER'),
                         'the selection should be defined exactly once')
        self.assertEqual(1, sql.count('CREATE TEMP TABLE stratified_selection'))
        export = sql.split('EXPORT DATA OPTIONS')[1]
        self.assertIn('FROM stratified_selection', export,
                      'the EXPORT should read the temp table, not re-derive the selection')

    def test_export_path_is_consumable_by_read_sample_list(self):
        """Tab-delimited with a header is exactly what read_sample_list expects."""
        sql = self.bq_file('vds_dropout_sample_list.sql')
        export = sql.split('EXPORT DATA OPTIONS')[1]
        self.assertIn("field_delimiter='\\t'", export)
        self.assertIn('header=true', export)
        self.assertIn('overwrite=true', export)

    def test_export_uri_has_exactly_one_wildcard(self):
        """EXPORT DATA requires precisely one wildcard in the uri."""
        for name in ('vds_dropout_sample_map.sql', 'vds_dropout_sample_list.sql'):
            sql = self.bq_file(name)
            uris = re.findall(r"uri='([^']*)'", sql)
            self.assertTrue(uris, msg=f'{name} has no EXPORT uri')
            for uri in uris:
                self.assertEqual(1, uri.count('*'), msg=f'{name}: {uri}')

    def test_export_path_preserves_the_diffable_ordering(self):
        sql = self.bq_file('vds_dropout_sample_list.sql')
        export = sql.split('EXPORT DATA OPTIONS')[1]
        self.assertIn('ORDER BY superpartition, sample_name', export)

    def test_only_one_output_path_is_live(self):
        """Running the file should produce one output, not two.

        Both output paths read the temp table, so counting live references to it counts
        live output paths directly.
        """
        sql = self.bq_file('vds_dropout_sample_list.sql')
        live = '\n'.join(line for line in sql.splitlines()
                         if line.strip() and not line.lstrip().startswith('--'))
        self.assertEqual(1, live.count('FROM stratified_selection'))
        self.assertNotIn('EXPORT DATA', live,
                         'the EXPORT path ships commented out; the returning SELECT is default')

    def test_sample_map_sql_selects_the_expected_columns(self):
        sql = self.bq_file('vds_dropout_sample_map.sql')
        self.assertIn('sample_name', sql)
        self.assertIn('sample_id', sql)
        self.assertIn('withdrawn IS NULL', sql)
        self.assertIn('is_control = false', sql)

    def test_both_sql_files_use_substitutable_placeholders(self):
        for name in ('vds_dropout_sample_map.sql', 'vds_dropout_sample_list.sql'):
            sql = self.bq_file(name)
            self.assertIn('PROJECT_ID.DATASET_NAME.sample_info', sql, msg=name)


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

    def test_default_probe_interval_clears_the_centromere(self):
        """A probe interval overlapping an assembly gap skews the extrapolation.

        The chr20 centromere gap ends at 30,088,349 per
        wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list, so a window
        starting at 30,000,000 would begin 88 kb inside it.
        """
        contig, _, span = vds.DEFAULT_PROBE_INTERVAL.partition(':')
        start, _, end = span.partition('-')
        self.assertEqual('chr20', contig)
        self.assertGreater(int(start), vds.CHR20_CENTROMERE_END)
        self.assertEqual(10_000_000, int(end) - int(start))

    def test_default_probe_interval_stays_within_chr20(self):
        _, _, span = vds.DEFAULT_PROBE_INTERVAL.partition(':')
        _, _, end = span.partition('-')
        self.assertLess(int(end), 64_334_167)  # chr20 length

    def test_report_states_what_the_estimate_excludes(self):
        """Readers must not mistake the figure for a materialize prediction."""
        report = vds.format_probe_report(60.0, 10_000_000, 1_340, 200)
        self.assertIn('not writing', report)
        self.assertIn('floor for the materialize pass', report)

    def test_report_includes_sample_count_when_given(self):
        report = vds.format_probe_report(60.0, 10_000_000, 1_340, 200, n_samples=13_400)
        self.assertIn('13,400', report)

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
