#!/usr/bin/env python3
"""Unit tests for the Hail-free parts of vds_dropout_scan.

The module guards its ``import hail`` so everything below runs without Hail installed.
What is covered here is the logic that would be silently wrong rather than loudly
broken: superpartition arithmetic that has to line up with the `vet_NNN` table names,
sample selection that has to be reproducible across runs and VDSes, bin boundaries, and
the argument validation that stops a misconfigured run before it books a cluster.

The Hail aggregation itself is not covered; it is exercised by running a scan over a
bounded interval against real data.
"""

import argparse
import ast
import contextlib
import io
import os
import pathlib
import re
import tempfile
import types
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


class TestWdlGeneratedSampleMap(unittest.TestCase):
    """The WDL generates the sample map itself, and holds the only copy of the query.

    Getting the sample universe wrong here is quiet: it changes which samples are screened
    without erroring, so the essentials are asserted rather than assumed. The predicates
    must match what GvsExtractAvroFilesForHail.wdl applies when exporting Avro, or the peer
    comparison is drawn against a different cohort than the VDS holds.
    """

    WDL = (pathlib.Path(__file__).resolve().parents[2]
           / 'wdl' / 'GvsValidateVdsCompleteness.wdl')
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


class TestMapCoverageIsEnforced(unittest.TestCase):
    """A VDS sample absent from the map must abort the scan, not be skipped.

    Screening part of a superpartition biases every peer comparison the detector makes, and
    silently: the numbers still look plausible. The check itself now runs inside Hail, over
    the annotated column table, so it is guarded here at the source level rather than by
    calling it.
    """

    SOURCE = pathlib.Path(__file__).resolve().parents[1] / 'vds_dropout_scan.py'

    def test_unmatched_columns_raise(self):
        body = self.SOURCE.read_text().split('def _screening_matrix', 1)[1]
        body = body.split('\ndef ', 1)[0]
        self.assertIn('if n_unmatched:', body)
        self.assertIn('raise ValueError', body)

    def test_the_reason_is_stated_where_it_is_enforced(self):
        body = self.SOURCE.read_text().split('def _screening_matrix', 1)[1]
        body = body.split('\ndef ', 1)[0]
        self.assertIn('bias', body)
        self.assertIn('silently', body)


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


class TestExecutorSummary(unittest.TestCase):
    """Cluster width is logged for the cost model, so it must degrade one field at a time.

    The first version built both figures inside one try block and iterated
    getExecutorMemoryStatus().keySet(), which py4j exposes as a Java object rather than a
    Python iterable. The TypeError discarded the task-slot count too, which had been
    working -- so a diagnostic meant to explain a slow run instead reported nothing.
    """

    class _MemStatus:
        def __init__(self, size): self._size = size
        def size(self): return self._size

    def _context(self, size=5, parallelism=1856, break_status=False):
        outer = self

        class Sc:
            def getExecutorMemoryStatus(inner):
                if break_status:
                    raise TypeError("'JavaObject' object is not iterable")
                return outer._MemStatus(size)

        class Jsc:
            def sc(inner): return Sc()

        return types.SimpleNamespace(_jsc=Jsc(), defaultParallelism=parallelism)

    def summary_with(self, context_or_raiser):
        original = vds.hl
        vds.hl = types.SimpleNamespace(spark_context=context_or_raiser)
        try:
            return vds.executor_summary()
        finally:
            vds.hl = original

    def test_reports_both_figures(self):
        text = self.summary_with(lambda: self._context(size=5, parallelism=1856))
        # The driver is in the map but runs no tasks, hence 4 rather than 5.
        self.assertIn('4 executor(s)', text)
        self.assertIn('~1856 task slot(s)', text)

    def test_task_slots_survive_a_broken_executor_count(self):
        """The regression: one unavailable figure must not take the other with it."""
        text = self.summary_with(lambda: self._context(break_status=True))
        self.assertIn('executor count unavailable', text)
        self.assertIn('~1856 task slot(s)', text)

    def test_missing_context_is_reported_not_raised(self):
        def raiser():
            raise RuntimeError('no context')
        text = self.summary_with(raiser)
        self.assertIn('unavailable', text)
        self.assertIn('no context', text)

    def test_never_raises(self):
        """This only annotates a log line; it must not be able to fail a run."""
        def raiser():
            raise RuntimeError('boom')
        for ctx in (lambda: self._context(), lambda: self._context(break_status=True), raiser):
            self.assertIsInstance(self.summary_with(ctx), str)

    def test_executor_count_never_goes_negative(self):
        """An empty map would otherwise report -1 executors."""
        text = self.summary_with(lambda: self._context(size=0))
        self.assertIn('0 executor(s)', text)


class TestContigCheckpointing(unittest.TestCase):
    """A genome-wide scan runs for hours, so it checkpoints per contig and resumes.

    The verification on merge is the load-bearing part: a resume that silently skipped a
    contig would emit a summary that looks complete and reports the rest as clean. Building
    a silent omission into a tool designed to detect silent omissions would be a poor
    trade for restartability.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.summary = os.path.join(self.tmp.name, 'summary.tsv')

    def write_shard(self, contig, rows, mark_done=True):
        shard, marker = vds.shard_paths(self.summary, contig)
        vds.write_lines(shard, vds.SUMMARY_HEADER, [
            f'{contig}\t{i * 50000 + 1}\t{i * 50000 + 50001}\t83\t7500'
            for i in range(rows)])
        if mark_done:
            vds.write_lines(marker, 'contig\trows', [f'{contig}\t{rows}'])
        return shard

    def args(self, vds_path='gs://bucket/r2.vds', mode='variants', bin_size=50_000):
        return types.SimpleNamespace(
            vds_path=vds_path, mode=mode, bin_size=bin_size, summary_path=self.summary)

    def test_marker_records_provenance(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args())
        with open(marker) as handle:
            header, row = handle.read().strip().split('\n')
        self.assertEqual(vds.MARKER_HEADER, header)
        self.assertIn('gs://bucket/r2.vds', row)
        self.assertIn('variants', row)
        self.assertIn('50000', row)

    def test_matching_marker_verifies(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args())
        vds.verify_marker(marker, 'chr1', self.args())   # must not raise

    def test_different_vds_aborts(self):
        """The hazard: same output_prefix, different VDS, silently wrong summary."""
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args(vds_path='gs://bucket/r2.vds'))
        with self.assertRaises(RuntimeError) as ctx:
            vds.verify_marker(marker, 'chr1', self.args(vds_path='gs://bucket/r3.vds'))
        self.assertIn('different run', str(ctx.exception))
        self.assertIn('r3.vds', str(ctx.exception))

    def test_different_bin_size_aborts(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args(bin_size=50_000))
        with self.assertRaises(RuntimeError):
            vds.verify_marker(marker, 'chr1', self.args(bin_size=10_000))

    def test_different_mode_aborts(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args(mode='variants'))
        with self.assertRaises(RuntimeError):
            vds.verify_marker(marker, 'chr1', self.args(mode='references'))

    def test_mismatch_error_says_how_to_recover(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args(vds_path='a'))
        with self.assertRaises(RuntimeError) as ctx:
            vds.verify_marker(marker, 'chr1', self.args(vds_path='b'))
        self.assertIn('Delete', str(ctx.exception))

    def test_legacy_marker_is_accepted_with_a_warning(self):
        """Deliberate: an in-flight scan's checkpoints stay usable across the upgrade."""
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_lines(marker, 'contig\trows', ['chr1\t42'])
        with contextlib.redirect_stdout(io.StringIO()) as out:
            vds.verify_marker(marker, 'chr1', self.args())
        self.assertIn('predates provenance recording', out.getvalue())
        self.assertIn('delete', out.getvalue().lower())

    def test_malformed_marker_aborts(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_lines(marker, vds.MARKER_HEADER, ['chr1\t42'])
        with self.assertRaises(RuntimeError) as ctx:
            vds.verify_marker(marker, 'chr1', self.args())
        self.assertIn('malformed', str(ctx.exception))

    def test_shard_and_marker_are_distinct_paths(self):
        shard, marker = vds.shard_paths('/x/summary.tsv', 'chr4')
        self.assertEqual('/x/summary.tsv.chr4', shard)
        self.assertEqual('/x/summary.tsv.chr4.done', marker)
        self.assertNotEqual(shard, marker)

    def test_marker_governs_resume_not_shard_existence(self):
        """A shard truncated mid-write must not be mistaken for a finished one."""
        shard = self.write_shard('chr1', 2, mark_done=False)
        _, marker = vds.shard_paths(self.summary, 'chr1')
        self.assertTrue(vds._path_exists(shard))
        self.assertFalse(vds._path_exists(marker))

    def test_merge_concatenates_all_shards(self):
        shards = [self.write_shard('chr1', 3), self.write_shard('chr2', 2)]
        total = vds.concatenate_shards(shards, self.summary, ['chr1', 'chr2'])
        self.assertEqual(5, total)

    def test_merged_file_has_one_header(self):
        shards = [self.write_shard('chr1', 2), self.write_shard('chr2', 2)]
        vds.concatenate_shards(shards, self.summary, ['chr1', 'chr2'])
        with open(self.summary) as handle:
            lines = [l for l in handle.read().split('\n') if l.strip()]
        self.assertEqual(vds.SUMMARY_HEADER, lines[0])
        self.assertEqual(1, sum(1 for l in lines if l == vds.SUMMARY_HEADER))
        self.assertEqual(5, len(lines))

    def test_missing_contig_aborts_the_merge(self):
        """The safety net: an incomplete screen must never look complete."""
        shards = [self.write_shard('chr1', 2)]
        with self.assertRaises(RuntimeError) as ctx:
            vds.concatenate_shards(shards, self.summary, ['chr1', 'chr2'])
        self.assertIn('chr2', str(ctx.exception))
        self.assertIn('incomplete', str(ctx.exception))

    def test_error_says_how_to_recover(self):
        shards = [self.write_shard('chr1', 2)]
        with self.assertRaises(RuntimeError) as ctx:
            vds.concatenate_shards(shards, self.summary, ['chr1', 'chr2'])
        self.assertIn('.done', str(ctx.exception))

    def test_shard_with_a_foreign_header_is_rejected(self):
        shard = os.path.join(self.tmp.name, 'summary.tsv.chr9')
        vds.write_lines(shard, 'something\telse', ['chr9\t1'])
        with self.assertRaises(ValueError) as ctx:
            vds.concatenate_shards([shard], self.summary, ['chr9'])
        self.assertIn('corrupt or from another run', str(ctx.exception))

    def test_empty_shard_is_tolerated_when_the_contig_is_not_expected(self):
        """A contig with no data at all yields an empty shard; only expectation matters."""
        shard = os.path.join(self.tmp.name, 'summary.tsv.chrY')
        vds.write_lines(shard, vds.SUMMARY_HEADER, [])
        self.assertEqual(0, vds.concatenate_shards([shard], self.summary, []))


class TestWriteRetry(unittest.TestCase):
    """The summary write is the last step of a multi-hour job."""

    def test_rows_are_materialized_so_a_retry_can_repeat_them(self):
        """A generator half-drained by a failed attempt would silently truncate the retry."""
        import inspect
        source = inspect.getsource(vds.write_lines)
        self.assertIn('list(rows)', source)

    def test_retry_budget_is_bounded_and_nonzero(self):
        self.assertGreater(vds.WRITE_ATTEMPTS, 1)
        self.assertLess(vds.WRITE_ATTEMPTS, 10)
        self.assertGreater(vds.WRITE_RETRY_DELAY_SECONDS, 0)

    def test_write_returns_the_row_count(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'out.tsv')
            self.assertEqual(3, vds.write_lines(path, 'h', ['a', 'b', 'c']))


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


class TestIntervalReadSemantics(unittest.TestCase):
    """read_vds(intervals=...) repartitions; it does not filter.

    Also the reason locus sampling was dropped entirely: with the reader pruning to native
    partitions, a genome-wide variant scan is under an hour and a reference scan a couple
    of hours, so there is nothing worth buying by reading only part of the genome -- and
    the first question this tool exists to answer is exhaustive by nature.

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

    def test_scan_requires_a_sample_map(self):
        """Superpartition membership is a function of sample_id, which a VDS lacks."""
        self.assert_rejected('--action', 'scan', '--vds-path', 'gs://x',
                             '--summary-path', 's.tsv', '--superpartitions-path', 'sp.tsv')

    def test_full_depth_requires_output(self):
        self.assert_rejected('--action', 'full-depth', '--vds-path', 'gs://x',
                             '--sample-map-path', 'm.tsv')

    def test_bad_bin_size_rejected(self):
        message = self.assert_rejected('--action', 'scan', '--vds-path', 'gs://x',
                                       '--summary-path', 's.tsv',
                                       '--superpartitions-path', 'sp.tsv',
                                       '--sample-map-path', 'm.tsv', '--bin-size', '0')
        self.assertIn('--bin-size must be at least 1', message)

    def test_intervals_are_repeatable(self):
        args = self.parse('--action', 'scan', '--vds-path', 'gs://x',
                          '--summary-path', 's.tsv', '--superpartitions-path', 'sp.tsv',
                          '--sample-map-path', 'm.tsv',
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
