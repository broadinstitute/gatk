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
import sys
import tempfile
import threading
import time
import types
import unittest
import unittest.mock

import vds_dropout_scan as vds


def variantstore_dir():
    """The `scripts/variantstore` tree, wherever this test run can see it.

    The docker-mode run in `run_python_unit_tests.sh` bind-mounts it and exports this
    variable to say where, because the working copy is not otherwise visible inside the
    container. Falling back to the path relative to this file keeps a bare
    `PYTHONPATH=. python3 test/test_vds_dropout_scan.py` working.
    """
    mounted = os.environ.get('GVS_VARIANTSTORE_DIR')
    return pathlib.Path(mounted) if mounted else pathlib.Path(__file__).resolve().parents[2]


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


class TestReferenceLengthField(unittest.TestCase):
    """Reference-block length is spelled LEN or END, and read_vds supplies both today.

    So this is insurance against a Hail version that stops synthesizing one, not a fix for a
    VDS we cannot read. The Hail expression building around the choice cannot be exercised
    off-cluster; the choice itself can, which is why it is a separate function.
    """

    def test_len_schema(self):
        self.assertEqual('LEN', vds.reference_length_field(['GQ', 'GT', 'LEN']))

    def test_end_schema(self):
        self.assertEqual('END', vds.reference_length_field(['END', 'GQ', 'GT']))

    def test_len_wins_when_both_present(self):
        self.assertEqual('LEN', vds.reference_length_field(['END', 'GQ', 'LEN']))

    def test_neither_is_an_error(self):
        with self.assertRaises(ValueError) as caught:
            vds.reference_length_field(['GQ', 'GT'])
        self.assertIn('neither LEN nor END', str(caught.exception))


class TestBinArithmetic(unittest.TestCase):

    def test_default_bin_size(self):
        """Pinned because it sets the detection floor and is baked into the Hail pass."""
        self.assertEqual(10_000, vds.DEFAULT_BIN_SIZE)

    def test_bin_start_is_one_based_inclusive(self):
        self.assertEqual(1, vds.bin_start_for(1, 50_000))
        self.assertEqual(1, vds.bin_start_for(50_000, 50_000))
        self.assertEqual(50_001, vds.bin_start_for(50_001, 50_000))

    def test_known_dropout_boundaries(self):
        """The real chr4 window's edges, to pin the quantization the report will show."""
        self.assertEqual(56_550_001, vds.bin_start_for(56_585_368, 50_000))
        self.assertEqual(57_000_001, vds.bin_start_for(57_035_833, 50_000))
        self.assertEqual(56_580_001, vds.bin_start_for(56_585_368, 10_000))
        self.assertEqual(57_030_001, vds.bin_start_for(57_035_833, 10_000))

    def test_bin_index_and_start_agree(self):
        for bin_size in (10_000, 50_000):
            for position in (1, 49_999, 50_000, 50_001, 56_585_368):
                index = vds.bin_index_for(position, bin_size)
                self.assertEqual(vds.bin_start_for(position, bin_size),
                                 index * bin_size + 1)


class TestWdlGeneratedSampleMap(unittest.TestCase):
    """The WDL generates the sample map itself, and holds the only copy of the query.

    Getting the sample universe wrong here is quiet: it changes which samples are screened
    without erroring, so the essentials are asserted rather than assumed. The predicates
    must match what GvsExtractAvroFilesForHail.wdl applies when exporting Avro, or the peer
    comparison is drawn against a different cohort than the VDS holds.
    """

    WDL = variantstore_dir() / 'wdl' / 'GvsValidateVdsCompleteness.wdl'
    def read(self, path):
        if not path.exists():
            # Left as a skip rather than a failure so the file can be run from a
            # checkout that does not have the WDL tree, but the docker-mode run does mount
            # it -- see `variantstore_dir`.
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


class TestPlaceholderOutputs(unittest.TestCase):
    """Both task outputs must exist and describe themselves accurately.

    Cromwell resolves task outputs regardless of which branch ran, so both files are created
    up front. Each has to describe itself rather than its sibling: an adjudicate_<mode>.sql
    whose text talks about reports is a plausibly named output containing something
    irrelevant, which reads as a real result.
    """

    WDL = variantstore_dir() / 'wdl' / 'GvsValidateVdsCompleteness.wdl'

    def test_scan_log_is_uploaded_however_the_task_exits(self):
        """The log has to survive the failures it is there to explain.

        Copied on the last line instead, after every step that can fail, it would be
        uploaded only by runs that succeed -- and under errexit a failing task would leave
        the log from an hours-long run in Cromwell's stdout alone. That is exactly
        backwards: a successful run barely needs its log.
        """
        body = self.scan_task()
        end = "fi' EXIT"
        self.assertIn(end, body, 'the log upload must be registered as an EXIT trap')
        trap = body[body.index('trap '):body.index(end) + len(end)]
        self.assertIn('scan_~{action}_~{mode}.log', trap)
        # `|| true` so a failed upload cannot replace the exit code that explains the failure.
        self.assertIn('|| true', trap)
        # And the log may not exist yet, since the trap is installed before it is written.
        self.assertIn('[[ -f scan.log ]]', trap)

    def test_scan_log_is_not_also_copied_at_the_end(self):
        """A second copy would be dead code implying the upload depends on reaching the end."""
        body = self.scan_task()
        after_trap = body[body.index("fi' EXIT"):]
        self.assertNotIn('gsutil cp scan.log', after_trap)

    def test_placeholder_does_not_deny_that_scan_produces_output(self):
        """The scan case is the only one in which a placeholder is ever read.

        The placeholders are written up front and overwritten by detect, so a placeholder
        surviving a scan means the task died in between. Interpolating the action into a
        sentence ending "only scan does" reads, for action=scan, as "action 'scan' does not
        produce one; only scan does" -- self-contradictory exactly when someone is trying to
        work out why their run failed.
        """
        text = self.wdl()
        self.assertNotIn("does not produce one; only scan does", text)
        # One sentence cannot serve both cases, so the wording is chosen at run time. The
        # branch has to sit immediately above the assignment for that to be what picks it.
        block = text[text.index('if [[ "~{action}" == "scan" ]]'):text.index('> report.tsv')]
        scan_case, _, other_case = block.partition('else')
        self.assertIn('placeholder=', scan_case)
        self.assertIn('did not get as far as the detect step', scan_case)
        # `~{action}` may appear here -- the log filename carries it -- but nothing in the
        # scan case may claim the action produces no such file, which is the whole bug.
        self.assertNotIn('does not produce', scan_case)
        # The other actions genuinely do not produce these files, so naming the action there
        # is the useful thing to say -- it is only the scan case that reads as a denial.
        self.assertIn("action '~{action}' does not produce this file", other_case)

    def test_placeholder_points_at_the_log_rather_than_itself(self):
        text = self.wdl()
        block = text[text.index('placeholder='):text.index('> report.tsv')]
        self.assertIn('scan_~{action}_~{mode}.log', block)
        self.assertIn('symptom, not the cause', block)

    def test_missing_detect_inputs_are_diagnosed_not_just_fatal(self):
        """A resumed scan need not rewrite the superpartition table, so this is reachable.

        Without the check the task died on a bare `gsutil cp` failure that named neither the
        object nor the reason, while the Hail log looked entirely healthy.
        """
        text = self.wdl()
        self.assertIn('gsutil -q stat', text)
        block = text[text.index('missing=()'):text.index('gsutil cp "~{summary_path}"')]
        self.assertIn('every contig already', block)
        self.assertIn('vds_dropout_detect.py', block)
        self.assertIn('output_prefix', block)

    def wdl(self):
        if not self.WDL.exists():
            self.skipTest(f'{self.WDL} not available in this test environment')
        return self.WDL.read_text()

    def scan_task(self):
        return self.wdl().split('task ScanVdsForDropouts', 1)[1]

    def test_placeholders_are_not_copies_of_each_other(self):
        body = self.scan_task()
        self.assertNotIn('cp report.tsv adjudicate.sql', body)

    def test_each_placeholder_names_its_own_purpose(self):
        body = self.scan_task().split('command <<<', 1)[1]
        head = body.split('bq query', 1)[0]
        self.assertIn('No findings report', head)
        self.assertIn('No adjudication SQL', head)

    def test_skipped_adjudication_is_explained_in_the_file(self):
        body = self.scan_task()
        self.assertIn('NO_SQL_GENERATED', body)
        self.assertIn('unproven', body)

    def test_skipped_adjudication_warns_on_stderr(self):
        """A silently missing adjudication leaves candidates looking settled."""
        body = self.scan_task()
        self.assertIn('WARNING: bq_project_id/bq_dataset_name were not supplied', body)

    def test_the_file_says_how_to_recover_without_re_running_the_scan(self):
        body = self.scan_task()
        self.assertIn('vds_dropout_detect.py', body)
        self.assertIn('No need to re-run the scan', body)

    def test_no_dead_existence_guard_on_the_sql_upload(self):
        """The placeholder guarantees the path, so `-f` could never be false."""
        body = self.scan_task()
        self.assertNotIn('if [[ -f ./adjudicate.sql ]]', body)


class TestOutputPrefixPaths(unittest.TestCase):
    """How the task turns `output_prefix` into GCS object names."""

    WDL = TestPlaceholderOutputs.WDL
    wdl = TestPlaceholderOutputs.wdl

    def test_output_prefix_slashes_are_normalized(self):
        """A prefix carrying a doubled slash produces objects nothing can read back.

        Every path is built by concatenation, so `gs://b/p/` gives `gs://b/p//summary_x.tsv`.
        The scan writes through `hl.hadoop_open`, and `org.apache.hadoop.fs.Path` collapses
        consecutive slashes, so the object lands at the single-slash name. gsutil does not
        collapse them -- a GCS object name is a literal string -- so the `gsutil cp` that
        feeds detect reported `No URLs matched` for a file that was sitting right there.
        Nothing about the run looked wrong until then, which is why this needs a test.
        """
        text = self.wdl()
        self.assertIn('sub(sub(output_prefix, "/+", "/"), "/$", "")', text)
        self.assertIn('sub(collapsed_output_prefix, "^gs:/", "gs://")', text)
        # No lookbehind in any pattern: Cromwell's `sub` is specified against POSIX ERE, and
        # an engine rejecting one would fail the workflow outright -- worse than the bug
        # being fixed. Checked against the `sub` calls, not the file, since the comment
        # explaining the choice necessarily quotes the construct it rejects.
        for line in text.splitlines():
            expression = line.split('#', 1)[0]
            if 'sub(' in expression:
                self.assertNotIn('(?<', expression)
        # Every task must receive the sanitized value; one raw passthrough reopens the hole.
        passthroughs = [line.strip() for line in text.splitlines()
                        if line.strip().startswith('output_prefix =')]
        self.assertTrue(passthroughs, 'expected at least one output_prefix passthrough')
        for passthrough in passthroughs:
            self.assertEqual('output_prefix = clean_output_prefix,', passthrough)

    @staticmethod
    def _clean(prefix: str) -> str:
        """The WDL's normalization, in Python. `^`, `+` and `$` mean the same in both."""
        collapsed = re.sub('/$', '', re.sub('/+', '/', prefix))
        return re.sub('^gs:/', 'gs://', collapsed)

    def test_normalization_leaves_exactly_one_slash_between_components(self):
        for given, want in [
            ('gs://bucket/p/', 'gs://bucket/p'),
            ('gs://bucket/p//', 'gs://bucket/p'),
            ('gs://bucket/p', 'gs://bucket/p'),
            # Interior, which a trailing-slash strip misses and which fails identically.
            ('gs://bucket/a//b', 'gs://bucket/a/b'),
            ('gs://bucket/a//b/', 'gs://bucket/a/b'),
            ('gs://bucket/', 'gs://bucket'),
            ('/tmp//out/', '/tmp/out'),
        ]:
            with self.subTest(prefix=given):
                self.assertEqual(want, self._clean(given))

    def test_normalization_preserves_the_scheme(self):
        """The collapse flattens `gs://` too, so the repair has to put it back."""
        self.assertTrue(self._clean('gs://bucket/p/').startswith('gs://'))

    def test_normalized_paths_survive_concatenation(self):
        """The property that actually matters: what the task builds must be readable."""
        for prefix in ['gs://b/p', 'gs://b/p/', 'gs://b/p//', 'gs://b/a//c/']:
            with self.subTest(prefix=prefix):
                built = self._clean(prefix) + '/summary_references.tsv'
                self.assertNotIn('//', built.removeprefix('gs://'))


class TestCommandBlockIndentation(unittest.TestCase):
    """This WDL indents its command blocks, so no line in them may sit at column zero.

    Cromwell dedents a `command <<< >>>` block by its common leading whitespace, so one
    line starting at column zero -- a wrapped string, a pasted comment -- drops that common
    prefix to nothing. Nothing is then stripped, the heredoc terminators keep indentation
    that bash requires them not to have, and the script dies with `unexpected end of file`
    reported at its last line, nowhere near the cause. `womtool validate` passes
    throughout, and the break only shows up once the task runs in the cloud, which is an
    expensive place to learn it.

    Only the cause is checked here, because it is a choice this file makes rather than a
    rule: writing the whole block at column zero is equally valid, and three other WDLs in
    the repo do exactly that. The consequences that are rules for every WDL -- terminators
    and whitespace-sensitive bodies landing at column zero after the dedent, and the inline
    Python compiling -- belong to `check-wdl-heredocs` and are asserted over the whole
    variantstore tree in test_check_wdl_heredocs.py.
    """

    WDL = TestPlaceholderOutputs.WDL
    wdl = TestPlaceholderOutputs.wdl

    def command_blocks(self):
        blocks, current = [], None
        for line in self.wdl().splitlines():
            if current is None:
                if re.search(r'command\s*<<<', line):
                    current = []
                continue
            if line.strip() == '>>>':
                blocks.append(current)
                current = None
                continue
            current.append(line)
        self.assertIsNone(current, 'unterminated command block')
        self.assertTrue(blocks, 'expected at least one command block')
        return blocks

    def test_no_command_line_starts_at_column_zero(self):
        """The cause, checked directly, because the symptom appears far from it."""
        for block in self.command_blocks():
            for line in block:
                if line.strip():
                    self.assertTrue(
                        line.startswith(' '),
                        f'line at column zero defeats the dedent: {line!r}')


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

    Building both figures inside one try block invites the opposite: py4j exposes
    getExecutorMemoryStatus().keySet() as a Java object rather than a Python iterable, so
    iterating it raises TypeError and takes the working task-slot count down with it -- a
    diagnostic meant to explain a slow run then reports nothing at all.
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


class TestWidthHeartbeat(unittest.TestCase):
    """Cluster width is a profile, not a number, so it is sampled rather than snapshotted.

    A reading taken before an aggregation starts is taken before autoscaling has seen the
    work and is systematically low; one at the end misses a slow ramp. Reading either as
    characteristic is an easy way to get a cost estimate wrong by a large factor, so the log
    records the shape instead.
    """

    def setUp(self):
        self.original = vds.executor_summary
        self.addCleanup(lambda: setattr(vds, 'executor_summary', self.original))

    def capture(self, sampler, duration=0.5, interval=0.1):
        vds.executor_summary = sampler
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            with vds.width_heartbeat('label', interval_seconds=interval):
                time.sleep(duration)
        return buffer.getvalue()

    def test_samples_more_than_once(self):
        output = self.capture(lambda: '590 executor(s), ~2360 task slot(s)')
        self.assertGreater(output.count('label'), 1)

    def test_reports_elapsed_and_width(self):
        output = self.capture(lambda: '590 executor(s), ~2360 task slot(s)')
        self.assertIn('min elapsed', output)
        self.assertIn('2360 task slot(s)', output)

    def test_thread_stops_when_the_block_exits(self):
        self.capture(lambda: 'steady')
        remaining = [t for t in threading.enumerate() if t.name == 'width-heartbeat']
        self.assertEqual([], remaining)

    def test_thread_is_a_daemon_so_it_cannot_hold_the_process_open(self):
        names = []
        vds.executor_summary = lambda: 'steady'
        with contextlib.redirect_stdout(io.StringIO()):
            with vds.width_heartbeat('label', interval_seconds=5):
                names = [(t.name, t.daemon) for t in threading.enumerate()
                         if t.name == 'width-heartbeat']
        self.assertEqual([('width-heartbeat', True)], names)

    def test_a_failing_sampler_neither_raises_nor_stops_the_heartbeat(self):
        """This annotates a log line; it has no business failing a ten-hour run."""
        def boom():
            raise RuntimeError('py4j unavailable')
        self.capture(boom)   # must not raise
        remaining = [t for t in threading.enumerate() if t.name == 'width-heartbeat']
        self.assertEqual([], remaining)

    def test_interval_default_is_sane_for_a_long_aggregation(self):
        self.assertGreaterEqual(vds.WIDTH_HEARTBEAT_SECONDS, 60)
        self.assertLessEqual(vds.WIDTH_HEARTBEAT_SECONDS, 1800)


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

    def args(self, vds_path='gs://bucket/r2.vds', mode='variants', bin_size=50_000,
             superpartition_size=vds.DEFAULT_SUPERPARTITION_SIZE, injections=()):
        return types.SimpleNamespace(
            vds_path=vds_path, mode=mode, bin_size=bin_size,
            superpartition_size=superpartition_size, summary_path=self.summary,
            injections=injections)

    def merge(self, shards, expected):
        """`concatenate_shards` with its progress output captured.

        Those lines matter -- the merge is otherwise silent for tens of minutes -- but
        interleaved with unittest's own output they make a real failure harder to read. One
        test below asserts on them directly.
        """
        with contextlib.redirect_stdout(io.StringIO()):
            return vds.concatenate_shards(shards, self.summary, expected)

    @contextlib.contextmanager
    def tiny_chunks(self, size=7):
        """Shrink the copy's chunk so the shards here span many of them.

        Every shard a test writes is a few hundred bytes, orders of magnitude under one
        real chunk, so without this the boundary-crossing paths -- the newline tally and
        the trailing-newline guard -- are never reached by the suite at all.

        Both constants, because the merge has two implementations: the byte copy and the
        Hail-handle fallback behind it. Shrinking only the one the default path reads would
        leave the other's boundary handling untested while the suite still passed.
        """
        originals = (vds.MERGE_CHUNK_BYTES, vds.MERGE_CHUNK_CHARACTERS)
        vds.MERGE_CHUNK_BYTES = size
        vds.MERGE_CHUNK_CHARACTERS = size
        try:
            yield
        finally:
            vds.MERGE_CHUNK_BYTES, vds.MERGE_CHUNK_CHARACTERS = originals

    def test_marker_records_provenance(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args())
        with open(marker) as handle:
            header, row = handle.read().strip().split('\n')
        self.assertEqual(vds.MARKER_HEADER, header)
        self.assertIn('gs://bucket/r2.vds', row)
        self.assertIn('variants', row)
        self.assertIn('50000', row)
        self.assertIn(str(vds.DEFAULT_SUPERPARTITION_SIZE), row)

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

    def test_different_superpartition_size_aborts(self):
        """Lowering it is the silent case, so it is the one the marker has to catch.

        Old superpartition ids are a subset of the ids a smaller size produces, so nothing
        downstream errors: the resumed contigs are simply grouped at the previous size and
        compared against sample counts taken at the new one.
        """
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args(superpartition_size=4000))
        with self.assertRaises(RuntimeError) as ctx:
            vds.verify_marker(marker, 'chr1', self.args(superpartition_size=2000))
        self.assertIn('superpartition_size', str(ctx.exception))

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

    def test_unrecognized_header_aborts(self):
        """A marker whose provenance cannot be read is not evidence that a shard is good."""
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_lines(marker, 'contig\trows', ['chr1\t42'])
        with self.assertRaises(RuntimeError) as ctx:
            vds.verify_marker(marker, 'chr1', self.args())
        self.assertIn('unrecognized marker header', str(ctx.exception))
        self.assertIn('Delete', str(ctx.exception))

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

    def test_injected_shard_is_refused_by_a_clean_resume(self):
        """The hazard this guards: an injected shard holds a hole no VDS has.

        A clean run resuming over it would merge that hole into its summary, flag it, and
        report a dropout in a VDS that does not have one -- the exact silent wrong answer
        the other provenance fields exist to prevent.
        """
        _, marker = vds.shard_paths(self.summary, 'chr1')
        injected = (vds.Injection('chr20', 1_000_000, 1_100_000, 83),)
        vds.write_marker(marker, 'chr1', 42, self.args(injections=injected))
        with self.assertRaises(RuntimeError) as caught:
            vds.verify_marker(marker, 'chr1', self.args())
        self.assertIn('injections', str(caught.exception))

    def test_clean_shard_is_refused_by_an_injected_resume(self):
        """Symmetric, and the one that would otherwise understate an injection's effect."""
        _, marker = vds.shard_paths(self.summary, 'chr1')
        vds.write_marker(marker, 'chr1', 42, self.args())
        injected = (vds.Injection('chr20', 1_000_000, 1_100_000, 83),)
        with self.assertRaises(RuntimeError):
            vds.verify_marker(marker, 'chr1', self.args(injections=injected))

    def test_matching_injection_verifies(self):
        _, marker = vds.shard_paths(self.summary, 'chr1')
        injected = (vds.Injection('chr20', 1_000_000, 1_100_000, 83),)
        vds.write_marker(marker, 'chr1', 42, self.args(injections=injected))
        vds.verify_marker(marker, 'chr1', self.args(injections=injected))  # must not raise

    def test_marker_governs_resume_not_shard_existence(self):
        """A shard truncated mid-write must not be mistaken for a finished one."""
        shard = self.write_shard('chr1', 2, mark_done=False)
        _, marker = vds.shard_paths(self.summary, 'chr1')
        self.assertTrue(vds._path_exists(shard))
        self.assertFalse(vds._path_exists(marker))

    def test_merge_concatenates_all_shards(self):
        shards = [self.write_shard('chr1', 3), self.write_shard('chr2', 2)]
        total = self.merge(shards, ['chr1', 'chr2'])
        self.assertEqual(5, total)

    def test_merged_file_has_one_header(self):
        shards = [self.write_shard('chr1', 2), self.write_shard('chr2', 2)]
        self.merge(shards, ['chr1', 'chr2'])
        with open(self.summary) as handle:
            lines = [l for l in handle.read().split('\n') if l.strip()]
        self.assertEqual(vds.SUMMARY_HEADER, lines[0])
        self.assertEqual(1, sum(1 for l in lines if l == vds.SUMMARY_HEADER))
        self.assertEqual(5, len(lines))

    def test_missing_contig_aborts_the_merge(self):
        """The safety net: an incomplete screen must never look complete."""
        shards = [self.write_shard('chr1', 2)]
        with self.assertRaises(RuntimeError) as ctx:
            self.merge(shards, ['chr1', 'chr2'])
        self.assertIn('chr2', str(ctx.exception))
        self.assertIn('incomplete', str(ctx.exception))

    def test_error_says_how_to_recover(self):
        shards = [self.write_shard('chr1', 2)]
        with self.assertRaises(RuntimeError) as ctx:
            self.merge(shards, ['chr1', 'chr2'])
        self.assertIn('.done', str(ctx.exception))

    def test_shard_with_a_foreign_header_is_rejected(self):
        shard = os.path.join(self.tmp.name, 'summary.tsv.chr9')
        vds.write_lines(shard, 'something\telse', ['chr9\t1'])
        with self.assertRaises(ValueError) as ctx:
            self.merge([shard], ['chr9'])
        self.assertIn('corrupt or from another run', str(ctx.exception))

    def test_empty_shard_is_tolerated_when_the_contig_is_not_expected(self):
        """A contig with no data at all yields an empty shard; only expectation matters."""
        shard = os.path.join(self.tmp.name, 'summary.tsv.chrY')
        vds.write_lines(shard, vds.SUMMARY_HEADER, [])
        self.assertEqual(0, self.merge([shard], []))

    def test_a_shard_whose_rows_are_for_another_contig_is_reported_missing(self):
        """The chunked copy reads one data row per shard, and that is the row it trusts.

        Cheaper than tallying every row's contig, and it fails in the safe direction: a
        shard holding the wrong contig raises here rather than passing quietly.
        """
        shard, _ = vds.shard_paths(self.summary, 'chr1')
        vds.write_lines(shard, vds.SUMMARY_HEADER, ['chr2\t1\t50001\t83\t7500'])
        with self.assertRaises(RuntimeError) as ctx:
            self.merge([shard], ['chr1'])
        self.assertIn('chr1', str(ctx.exception))

    def test_a_shard_without_a_trailing_newline_does_not_glue_rows_together(self):
        """Both chunk sizes, because a truncated shard can end anywhere in a chunk."""
        for label, chunks in (('one chunk', contextlib.nullcontext()),
                              ('many chunks', self.tiny_chunks())):
            with self.subTest(label):
                first = os.path.join(self.tmp.name, 'summary.tsv.chr1')
                with open(first, 'w') as handle:
                    handle.write(f'{vds.SUMMARY_HEADER}\n'
                                 'chr1\t1\t50001\t83\t1\n'
                                 'chr1\t50001\t100001\t83\t2')
                second = self.write_shard('chr2', 1)
                with chunks:
                    total = self.merge([first, second], ['chr1', 'chr2'])
                self.assertEqual(3, total)
                with open(self.summary) as handle:
                    rows = handle.read().splitlines()
                self.assertEqual([vds.SUMMARY_HEADER], rows[:1])
                self.assertEqual(4, len(rows))
                self.assertTrue(rows[3].startswith('chr2\t'), rows[3])

    def test_chunk_boundaries_do_not_change_the_result(self):
        """A boundary may fall mid-row, so the copy must not interpret what it copies."""
        shards = [self.write_shard('chr1', 5), self.write_shard('chr2', 4)]
        with self.tiny_chunks():
            total = self.merge(shards, ['chr1', 'chr2'])
        self.assertEqual(9, total)
        with open(self.summary) as handle:
            rows = handle.read().splitlines()
        self.assertEqual(10, len(rows))
        self.assertEqual(vds.SUMMARY_HEADER, rows[0])
        self.assertEqual(5, sum(1 for row in rows if row.startswith('chr1\t')))
        self.assertEqual(4, sum(1 for row in rows if row.startswith('chr2\t')))

    def test_the_merge_reports_progress_as_it_goes(self):
        """It is the longest non-Hail step, so silence there is the most expensive kind.

        A step with no Hail UI behind it is the one whose progress can only come from the
        script, and the question silence prompts -- whether the job has stalled -- is one no
        output could answer.
        """
        shards = [self.write_shard('chr1', 2), self.write_shard('chr2', 2)]
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            vds.concatenate_shards(shards, self.summary, ['chr1', 'chr2'])
        output = printed.getvalue()
        for shard in shards:
            self.assertIn(os.path.basename(shard), output)
        self.assertIn('4 rows', output)


class TestMergeAvoidsHail(unittest.TestCase):
    """The merge copies bytes directly rather than through `hl.hadoop_open`.

    Reading the summary back through Hail's text handles ran at roughly 180 kB/s, because
    every read and write crosses the Python/JVM boundary 8 KB at a time and the buffer
    cannot be enlarged -- see `TestHadoopOpenIsNotBufferTuned`. Copying 1.3 GB that way
    took just under two hours, which is what pushed the Foxtrot references scan past its
    cluster TTL and discarded ten hours of completed work. So the fast path has to stay
    the default, and the fallback has to stay correct for when it is not available.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.summary = os.path.join(self.tmp.name, 'summary.tsv')

    def write_shard(self, contig, rows):
        shard, _ = vds.shard_paths(self.summary, contig)
        vds.write_lines(shard, vds.SUMMARY_HEADER, [
            f'{contig}\t{i * 50000 + 1}\t{i * 50000 + 50001}\t83\t7500'
            for i in range(rows)])
        return shard

    def test_local_paths_never_reach_google_cloud_storage(self):
        """A local run must not import the package, let alone construct a client.

        `storage.Client()` goes looking for credentials, so building one for a merge of
        two files in a temp directory would make the suite fail on any machine without
        application default credentials -- and would make it talk to GCP on machines with
        them.
        """
        shard = self.write_shard('chr1', 2)
        with unittest.mock.patch.dict(sys.modules, {'google.cloud': None}):
            openers = vds._binary_openers([self.summary, shard])
        self.assertIsNotNone(openers, 'a local merge must not need the GCS client')
        read, write = openers
        with read(shard) as handle:
            self.assertIsInstance(handle.read(1), bytes)
        target = os.path.join(self.tmp.name, 'out')
        with write(target) as handle:
            handle.write(b'x')
        self.assertEqual(1, os.path.getsize(target))

    def test_a_gcs_path_anywhere_requires_the_gcs_client(self):
        """One remote path is enough to need it: shards and summary need not agree on
        scheme, and a merge that read remote shards into a local summary would still be
        moving bulk bytes over the network. The converse -- that the client is actually
        used when it imports -- cannot be asserted here, since constructing one goes
        looking for credentials."""
        with unittest.mock.patch.dict(sys.modules, {'google.cloud': None}):
            self.assertIsNone(vds._binary_openers(['gs://bucket/summary.tsv']))
            self.assertIsNone(vds._binary_openers([self.summary, 'gs://bucket/s.chr1']))

    def test_bucket_and_object_are_split(self):
        self.assertEqual(('b', 'a/summary.tsv'), vds._gcs_bucket_and_name(
            'gs://b/a/summary.tsv'))

    def test_a_bucket_without_an_object_is_an_error(self):
        """Rather than composing a request for an object named the empty string."""
        for path in ('gs://bucket', 'gs://bucket/', 'gs:///object'):
            with self.assertRaises(ValueError):
                vds._gcs_bucket_and_name(path)

    def test_the_fallback_produces_the_same_bytes(self):
        """It is reached only on a cluster missing the package, so equivalence is asserted
        here or nowhere. A merge that silently differed between the two would make the
        summary depend on which packages the cluster happened to have."""
        shards = [self.write_shard('chr1', 5), self.write_shard('chr2', 4)]
        with contextlib.redirect_stdout(io.StringIO()):
            fast = vds.concatenate_shards(shards, self.summary, ['chr1', 'chr2'])
        with open(self.summary, 'rb') as handle:
            fast_bytes = handle.read()

        with contextlib.redirect_stdout(io.StringIO()):
            slow = vds._concatenate_shards_via_hail(shards, self.summary, ['chr1', 'chr2'])
        with open(self.summary, 'rb') as handle:
            slow_bytes = handle.read()

        self.assertEqual(9, fast)
        self.assertEqual(fast, slow)
        self.assertEqual(fast_bytes, slow_bytes)

    def test_the_fallback_still_enforces_coverage(self):
        shards = [self.write_shard('chr1', 3)]
        with self.assertRaises(RuntimeError) as caught:
            with contextlib.redirect_stdout(io.StringIO()):
                vds._concatenate_shards_via_hail(shards, self.summary, ['chr1', 'chr2'])
        self.assertIn('chr2', str(caught.exception))

    def test_a_missing_package_warns_rather_than_failing(self):
        """The merge is the last step of a job that has already run for hours, so losing
        it to an ImportError would discard all of that. Degrade loudly instead."""
        shards = [self.write_shard('chr1', 3)]
        printed = io.StringIO()
        with unittest.mock.patch.object(vds, '_binary_openers', return_value=None):
            with contextlib.redirect_stdout(printed):
                total = vds.concatenate_shards(shards, self.summary, ['chr1'])
        self.assertEqual(3, total)
        self.assertIn('google-cloud-storage', printed.getvalue())

    def test_a_corrupt_header_is_reported_as_text(self):
        """The byte path must not surface a `b'...'` repr in an operator-facing message."""
        shard, _ = vds.shard_paths(self.summary, 'chr1')
        vds.write_lines(shard, 'contig\tsomething\telse', ['chr1\t1\t2'])
        with self.assertRaises(ValueError) as caught:
            with contextlib.redirect_stdout(io.StringIO()):
                vds.concatenate_shards([shard], self.summary, ['chr1'])
        message = str(caught.exception)
        self.assertIn('contig\\tsomething\\telse', message)
        self.assertNotIn("b'", message)


class TestHadoopOpenIsNotBufferTuned(unittest.TestCase):
    """A guard at the source level, because no test here can reach the code it guards.

    `_open_read` and `_open_write` only call `hl.hadoop_open` for a `gs://` path, so every
    test in this file takes the local `open` branch. Raising `buffer_size` looked like free
    throughput for the merge and instead made Hail's reader overrun its destination --
    `ValueError: memoryview assignment: lvalue and rvalue have different structures` on the
    first read -- which surfaced only on a cluster, after a 24-contig scan had finished.
    """

    SOURCE = pathlib.Path(__file__).resolve().parents[1] / 'vds_dropout_scan.py'

    def test_no_buffer_size_is_passed_to_hadoop_open(self):
        for line in self.SOURCE.read_text().splitlines():
            if 'hadoop_open(' in line:
                self.assertNotIn('buffer_size', line, line.strip())

    def test_the_reason_is_recorded_where_it_would_be_reintroduced(self):
        body = self.SOURCE.read_text().split('def _open_read', 1)[1].split('\ndef ', 1)[0]
        self.assertIn('buffer_size', body)
        self.assertIn('memoryview', body)


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

    Also the reason locus sampling was dropped entirely: a genome-wide scan is an overnight
    job at any sampling rate -- 5 h 01 m of Hail aggregation for variants and about 10 hours
    for references on Foxtrot r2, the merge and the surrounding tasks adding well under an
    hour -- so reading only part of the genome does not change how the tool is run,
    and the first question it exists to answer is exhaustive by nature.

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


class TestInjectionParsing(unittest.TestCase):
    """`--inject-dropout` is the only way to watch a hole become a depleted count.

    r2's two dropouts are variant-only, so nothing in the real data exercises references
    mode end to end. The spec is parsed rather than trusted because a silently misread
    window would produce a run that looks like a clean scan and proves nothing.
    """

    def test_parses_a_plain_spec(self):
        self.assertEqual(
            vds.Injection('chr20', 1000000, 1100000, 83),
            vds.parse_injection('chr20:1000000-1100000:83'))

    def test_accepts_grouped_digits(self):
        """A ten-digit coordinate typed by hand is the likeliest thing to get wrong."""
        self.assertEqual(
            vds.parse_injection('chr20:1000000-1100000:83'),
            vds.parse_injection('chr20:1_000_000-1,100,000:83'))

    def test_round_trips_through_str(self):
        """`__str__` is what lands in the marker, so it has to be re-parseable."""
        injection = vds.Injection('chr19', 40_000_001, 40_650_000, 64)
        self.assertEqual(injection, vds.parse_injection(str(injection)))

    def test_rejects_a_malformed_spec(self):
        for spec in ('chr20:1000000-1100000', 'chr20-1000000:83', '', 'chr20:a-b:83'):
            with self.subTest(spec=spec), self.assertRaises(ValueError):
                vds.parse_injection(spec)

    def test_rejects_an_inverted_window(self):
        with self.assertRaises(ValueError) as caught:
            vds.parse_injection('chr20:1100000-1000000:83')
        self.assertIn('past end', str(caught.exception))

    def test_no_flag_means_no_injection(self):
        self.assertEqual((), vds.parse_injection_list(None))
        self.assertEqual((), vds.parse_injection_list([]))

    def test_uninjected_matrix_is_untouched(self):
        """The no-op path must not reach Hail, so a clean scan is unaffected off-cluster."""
        sentinel = object()
        self.assertIs(sentinel, vds.apply_injections(sentinel, ()))

    def test_provenance_of_no_injection_is_not_empty(self):
        """An empty field would be ambiguous with a truncated marker row."""
        self.assertEqual('-', vds.injection_provenance(()))

    def test_bad_spec_is_a_parser_error_not_a_traceback(self):
        parser = vds.build_parser()
        args = parser.parse_args([
            '--action', 'scan', '--vds-path', 'gs://b/x.vds',
            '--sample-map-path', 'gs://b/map.tsv', '--summary-path', 'gs://b/s.tsv',
            '--superpartitions-path', 'gs://b/sp.tsv',
            '--inject-dropout', 'nonsense'])
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            vds.validate_args(args, parser)

    def test_validate_args_populates_injections(self):
        parser = vds.build_parser()
        args = parser.parse_args([
            '--action', 'scan', '--vds-path', 'gs://b/x.vds',
            '--sample-map-path', 'gs://b/map.tsv', '--summary-path', 'gs://b/s.tsv',
            '--superpartitions-path', 'gs://b/sp.tsv',
            '--inject-dropout', 'chr20:1000000-1100000:83'])
        vds.validate_args(args, parser)
        self.assertEqual((vds.Injection('chr20', 1000000, 1100000, 83),), args.injections)

    def test_a_clean_run_has_an_empty_injection_tuple(self):
        parser = vds.build_parser()
        args = parser.parse_args([
            '--action', 'scan', '--vds-path', 'gs://b/x.vds',
            '--sample-map-path', 'gs://b/map.tsv', '--summary-path', 'gs://b/s.tsv',
            '--superpartitions-path', 'gs://b/sp.tsv'])
        vds.validate_args(args, parser)
        self.assertEqual((), args.injections)
