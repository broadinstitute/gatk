#!/usr/bin/env python3
"""
Scan a GVS Hail VDS for "rectangle" data dropouts and emit a small summary table.

A rectangle dropout is a contiguous genomic window in which one GVS superpartition has
little or no data while every other superpartition has the usual amount.  This is the
shape produced when a single Avro export shard is lost, truncated, or never read: the
`EXPORT DATA ... ORDER BY location` in `GvsExtractAvroFilesForHail.wdl` means each
numbered output file holds a contiguous location range for exactly one superpartition.

This script does the expensive Hail work and nothing else.  All judging lives in
`vds_dropout_detect.py`, which reads the summary this produces.  The split means
thresholds can be re-tuned offline for free, and the judging logic can be unit tested
without Hail or a real VDS.

Actions
-------
``--action scan``
    Aggregate a VDS into per-(bin, superpartition) totals and write the summary and
    superpartition TSVs that `vds_dropout_detect.py` consumes.  Reads every partition
    overlapping the requested contigs or intervals; there is no sampling of loci, so the
    result supports an exhaustive claim about what is and is not missing.

``--action full-depth``
    For named intervals and superpartitions, count per-sample data and report how many
    samples carry nothing.  Meant for the handful of windows the screen flags, not for the
    genome: the screen aggregates to (bin, superpartition), so this is what turns a flagged
    rectangle into a per-sample tally.

Metrics
-------
``--mode variants``
    Counts defined ``variant_data`` entries.  One entry per (sample, locus), which maps
    one-to-one onto rows of the corresponding `vet_%` table, so the summary is directly
    comparable to a BigQuery row count.  Note that a variant call with ``call_GQ = 0``
    has a missing ``LGT`` but a present entry, so definedness is tested on ``LA``.

``--mode references``
    Sums reference-block coverage.  Hail spells the block length two ways -- ``LEN``
    directly, or ``END`` as the inclusive last position -- and ``read_vds`` supplies both,
    so either serves; ``reference_length_field`` says which one this uses and why.  Each block is attributed entirely to the bin holding its
    start locus; with 10 kb bins and a 1000 bp maximum block length the edge error is
    under 10%, and it is a uniform bias that cancels in the superpartition-versus-peers
    comparison the detector performs.

    Reference dropouts leave blocks *absent* rather than GQ 0.  AoU ingests with
    ``drop_state = "ZERO"`` so GQ0 blocks never reach BigQuery, and nothing in VDS
    construction back-fills gaps -- `import_gvs.py` builds `reference_data` purely from
    the Avro rows, and `update_reference_data_ploidy` only annotates `GT` on entries that
    already exist.

Sizing a run
------------
A partition is the unit of work: one task streams one partition end to end, so cost is
per-partition time times partitions divided by concurrent tasks, and is **not**
proportional to genomic span.

Both terms are cheap to obtain without a trial run.  Partition geometry is metadata::

    vds = hl.vds.read_vds(path)
    vds.variant_data.n_partitions()      # 119,189 for Foxtrot r2
    vds.reference_data.n_partitions()

Per-partition time comes from the scan itself: contigs are checkpointed, so the first one
to finish reports its own partition count and duration, and that work counts toward the
result rather than being discarded.  There is deliberately no separate timing probe, because
a measurement of this cannot be both cheap and representative -- a narrow interval measures
near-serial throughput, and one wide enough to exercise the cluster is already a substantial
fraction of the real run.

Every sample is screened
------------------------
There is no sampling, of samples or of loci.  A full-width pass is not cheap: measured on
Foxtrot r2, the Hail aggregation over a 535K-sample VDS took 5 h 01 m for variants and about
10 hours for references, with the shard merge and the tasks either side of it adding under
half an hour on top of that.  But it is an
overnight job at any sampling rate, so a subset would not change how the tool is run, and the
first question it exists to settle is exhaustive by nature.  Screening everything is also
simpler than screening a stratified subset, strictly more sensitive -- a sampled screen cannot
see a handful of individually lost samples -- and needs none of the machinery that makes a
sampled screen safe to compare across VDSes.

Sample map
----------
Superpartition membership comes from a TSV of sample name to sample ID, since a VDS knows
only names.  ``GvsValidateVdsCompleteness.wdl`` generates it from ``sample_info`` when given
a BigQuery project and dataset, and its ``GenerateSampleMap`` task is the one place the
query lives.

To supply one by hand instead, produce two columns with a header::

    sample_name<TAB>sample_id

Either tab- or comma-separated input is accepted, so ``bq query --format=csv`` output works
as-is.  The sample universe should match what the Avro export used -- non-withdrawn,
non-control -- or the peer comparison is drawn against a different cohort than the VDS
holds.

The map must cover every sample in the VDS.  A VDS sample missing from it is a hard error
rather than a skip, because screening part of a superpartition biases the peer comparison
the detector relies on, and would do so silently.  Samples in the map but absent from the
VDS are ignored.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime
import os
import re
import sys
import threading
import time
from collections import defaultdict
from typing import IO, Iterable, Mapping, NamedTuple, Sequence

try:
    import hail as hl
except ModuleNotFoundError:  # pragma: no cover - exercised only off-cluster
    hl = None

DEFAULT_SUPERPARTITION_SIZE = 4000
DEFAULT_BIN_SIZE = 10_000
# Below this many partitions the aggregation cannot use the cluster, which on an AoU-scale
# VDS turns a nominally small interval into hours of single-threaded streaming.
MIN_HEALTHY_PARTITIONS = 8

# How often to sample cluster width during a long aggregation. A contig takes tens of
# minutes, so this is a handful of lines each.
WIDTH_HEARTBEAT_SECONDS = 300
MODES = ('variants', 'references')
ACTIONS = ('scan', 'full-depth')

# GRCh38 primary contigs GVS encodes. Alt and decoy contigs carry no GVS location
# encoding and are not screened.
GVS_CONTIGS = tuple([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])


SUMMARY_HEADER = 'contig\tbin_start\tbin_end\tsuperpartition\tobserved'
SUPERPARTITION_HEADER = 'superpartition\tn_samples'
FULL_DEPTH_HEADER = (
    'contig\tstart\tend\tsuperpartition\tn_samples\tn_zero\tfraction_zero\ttotal_observed'
)


# ---------------------------------------------------------------------------
# Pure helpers -- no Hail, unit tested directly
# ---------------------------------------------------------------------------


def superpartition_for(sample_id: int, superpartition_size: int = DEFAULT_SUPERPARTITION_SIZE) -> int:
    """Return the 1-based GVS superpartition holding ``sample_id``.

    Mirrors ``CAST(CEIL(sample_id / 4000.0) AS INT64)`` in
    `GvsExtractAvroFilesForHail.wdl`, so superpartition numbering here lines up with the
    `vet_NNN` and `ref_ranges_NNN` table names.
    """
    if sample_id < 1:
        raise ValueError(f"sample_id must be positive, got {sample_id}")
    return (sample_id + superpartition_size - 1) // superpartition_size


def bin_start_for(position: int, bin_size: int = DEFAULT_BIN_SIZE) -> int:
    """Return the inclusive 1-based start of the bin containing ``position``."""
    return ((position - 1) // bin_size) * bin_size + 1


def bin_index_for(position: int, bin_size: int = DEFAULT_BIN_SIZE) -> int:
    """Return the 0-based index of the bin containing ``position``."""
    return (position - 1) // bin_size


def parse_sample_map(lines: Iterable[str]) -> dict[str, int]:
    """Parse a ``sample_name<TAB>sample_id`` table into a dict.

    Accepts comma-separated input too, since the documented ``bq query --format=csv``
    recipe produces commas if the caller forgets the ``tr``.
    """
    result: dict[str, int] = {}
    header_seen = False
    for lineno, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line:
            continue
        fields = line.split('\t') if '\t' in line else line.split(',')
        if len(fields) < 2:
            raise ValueError(f"sample map line {lineno}: expected 2 columns, got {len(fields)}")
        name, value = fields[0].strip(), fields[1].strip()
        if not header_seen:
            header_seen = True
            if value.lower() == 'sample_id':
                continue
            raise ValueError(
                "sample map must start with a 'sample_name<TAB>sample_id' header, found "
                f"{name!r}, {value!r}")
        try:
            result[name] = int(value)
        except ValueError as exc:
            raise ValueError(f"sample map line {lineno}: {exc}") from exc
    if not result:
        raise ValueError('sample map contained no rows')
    return result


def parse_interval_list(values: Sequence[str] | None) -> list[str]:
    """Flatten repeated and comma-separated interval arguments into one list.

    Repeatable ``--intervals`` is convenient by hand but cannot be expressed through the
    arguments JSON, which holds one value per key, so a single comma-separated string has
    to work too.
    """
    if not values:
        return []
    intervals: list[str] = []
    for value in values:
        intervals.extend(part.strip() for part in value.split(',') if part.strip())
    return intervals


def parse_contig_list(value: str | None) -> list[str]:
    """Parse a comma-separated contig list, defaulting to all GVS primary contigs."""
    if not value:
        return list(GVS_CONTIGS)
    contigs = [c.strip() for c in value.split(',') if c.strip()]
    unknown = [c for c in contigs if c not in GVS_CONTIGS]
    if unknown:
        raise ValueError(f"not GVS-encoded contig(s): {unknown}; expected from {GVS_CONTIGS}")
    return contigs


def format_summary_rows(
        totals: Mapping[tuple[str, int], Mapping[int, float]],
        bin_size: int,
        contig_order: Sequence[str] = GVS_CONTIGS,
) -> list[str]:
    """Render aggregated totals as summary TSV lines, skipping zero cells.

    Zero cells are omitted rather than written: they are the overwhelming majority in a
    region, and `vds_dropout_detect.py` restores them from the superpartition table.  Rows are ordered by contig then position so the output is
    stable and diffable across runs.
    """
    rank = {contig: i for i, contig in enumerate(contig_order)}
    rows: list[str] = []
    for contig, bin_start in sorted(totals, key=lambda k: (rank.get(k[0], len(rank)), k[1])):
        per_superpartition = totals[(contig, bin_start)]
        for superpartition in sorted(per_superpartition):
            observed = per_superpartition[superpartition]
            if observed:
                rows.append(
                    f'{contig}\t{bin_start}\t{bin_start + bin_size}\t{superpartition}\t'
                    f'{observed:.0f}'
                )
    return rows


def format_superpartition_rows(chosen: Mapping[int, Sequence[str]]) -> list[str]:
    """Render the superpartition -> sampled-count table as TSV lines."""
    return [f'{sp}\t{len(chosen[sp])}' for sp in sorted(chosen)]


def announce(message: str) -> None:
    """Print a progress line immediately.

    ``flush`` matters: stdout is block-buffered when it is not a terminal, so without it a
    stalled run can leave its most recent progress line sitting in the buffer, and the log
    then understates how far execution actually got.
    """
    print(f'[{datetime.datetime.now().isoformat(timespec="seconds")}] {message}', flush=True)


@contextlib.contextmanager
def step(description: str):
    """Announce a unit of work and report how long it took.

    Every Hail action is wrapped so a stall names the operation responsible. Without this
    the script printed its intent and then nothing until the whole aggregation finished,
    which meant a hang had to be localized by pulling a diagnostic bundle off the cluster
    and reading jstacks -- for information the job could simply have logged.
    """
    announce(f'{description}: starting')
    started = time.monotonic()
    try:
        yield
    except BaseException:
        announce(f'{description}: FAILED after {time.monotonic() - started:,.1f}s')
        raise
    else:
        announce(f'{description}: done in {time.monotonic() - started:,.1f}s')


def _require_hail() -> None:
    if hl is None:
        raise RuntimeError(
            'hail is not importable; the Hail-backed actions of this script must run on '
            'a Dataproc cluster. The pure helpers can be imported without it.'
        )


def _open_write(path: str) -> IO[str]:
    """Open a local or GCS path for writing text."""
    if hl is not None and path.startswith('gs://'):
        return hl.hadoop_open(path, 'w')
    return open(path, 'wt')


def _path_exists(path: str) -> bool:
    """Whether a local or GCS path exists."""
    if hl is not None and path.startswith('gs://'):
        return bool(hl.hadoop_exists(path))
    return os.path.exists(path)


def _open_read(path: str) -> IO[str]:
    """Open a local or GCS path for reading text.

    Do not pass ``buffer_size`` to `hl.hadoop_open` here, however tempting it looks for a
    bulk copy. Hail's reader asks the JVM for ``buffer_size`` bytes but assigns them into a
    destination sized by the layer above it, so anything much over the 8 KB default
    overruns on the first read: `hadoop_fs.py` raises ``ValueError: memoryview assignment:
    lvalue and rvalue have different structures``. That killed a completed 24-contig Delta
    scan at the merge, and no local test can catch it -- this branch is only taken for a
    `gs://` path.
    """
    if hl is not None and path.startswith('gs://'):
        return hl.hadoop_open(path, 'r')
    return open(path, 'rt')


WRITE_ATTEMPTS = 3
WRITE_RETRY_DELAY_SECONDS = 15


def write_lines(path: str, header: str, rows: Iterable[str]) -> int:
    """Write a header plus rows, returning the number of rows written.

    Retried, because this is the last step of a job that may have run for hours and a
    transient GCS error here would discard all of it. ``rows`` is materialized first so a
    retry cannot resume from a half-drained generator.
    """
    materialized = list(rows)
    last_error = None
    for attempt in range(WRITE_ATTEMPTS):
        try:
            with _open_write(path) as handle:
                handle.write(header + '\n')
                for row in materialized:
                    handle.write(row + '\n')
            return len(materialized)
        except Exception as exc:  # pragma: no cover - needs a real filesystem failure
            last_error = exc
            if attempt + 1 < WRITE_ATTEMPTS:
                announce(f'Writing {path} failed ({exc}); retrying in '
                         f'{WRITE_RETRY_DELAY_SECONDS}s')
                time.sleep(WRITE_RETRY_DELAY_SECONDS)
    raise RuntimeError(
        f'Could not write {path} after {WRITE_ATTEMPTS} attempts') from last_error


def read_sample_map(path: str) -> dict[str, int]:
    with _open_read(path) as handle:
        return parse_sample_map(handle)


def build_intervals(contigs: Sequence[str], interval_strings: Sequence[str]):
    """Build Hail locus intervals from explicit interval strings or whole contigs.

    Returns concrete ``hail.utils.Interval`` objects, not ``IntervalExpression``s.  That
    distinction is load-bearing: ``hl.parse_locus_interval`` yields an expression, but the
    ``intervals`` argument of ``read_vds`` / ``read_matrix_table`` requires evaluated
    intervals, and handing it expressions fails deep inside Hail with

        TypeError: __init__: parameter 'includes_start': expected bool,
                   found BooleanExpression

    which points at Hail internals rather than at the caller.

    The whole list is evaluated in a single ``hl.eval`` rather than one per interval,
    since the default is every primary contig and each ``hl.eval`` is its own job.
    """
    _require_hail()
    strings = list(interval_strings) if interval_strings else list(contigs)
    if not strings:
        return []
    return list(hl.eval(hl.array(
        [hl.parse_locus_interval(s, reference_genome='GRCh38') for s in strings])))


def total_partition_count(vds_path: str, mode: str) -> int:
    """Partitions in the whole VDS component, ignoring any interval subset.

    Metadata only, so it costs seconds. This is the denominator the cost model needs: a
    full run's wall clock is roughly per-partition time multiplied by this and divided by
    the tasks running concurrently. Without it, a timing from an interval subset cannot be
    turned into an estimate of anything.
    """
    _require_hail()
    vds = hl.vds.read_vds(vds_path)
    matrix = vds.variant_data if mode == 'variants' else vds.reference_data
    return int(matrix.n_partitions())


def read_and_subset_vds(vds_path: str, intervals):
    """Read a VDS and restrict it to ``intervals``, preserving its native partitioning.

    Deliberately **not** ``read_vds(vds_path, intervals=...)``. That argument is a
    *repartitioning* specification rather than a filter: N intervals produce exactly N
    partitions. Passing a single 10 Mb interval therefore collapsed roughly 380 of this
    VDS's native partitions into one, and since a partition is the unit of work, one task
    had to stream all of it end to end -- observed as a nine-hour RUNNABLE task with the
    rest of the cluster idle.

    It would have been worse for a real scan: ``--contigs`` supplies one interval per
    contig, so a genome-wide run would have used 24 partitions instead of the 119,189 the
    VDS actually has.

    Reading first and filtering afterwards keeps the native partitioning and prunes to the
    partitions overlapping the request, which is both correct and parallel.
    """
    _require_hail()
    vds = hl.vds.read_vds(vds_path)
    if intervals:
        # split_reference_blocks=False: reference blocks straddling an interval edge are
        # left whole. Splitting them costs work and would skew the covered-base metric.
        vds = hl.vds.filter_intervals(vds, intervals, split_reference_blocks=False)
    return vds


def sniff_delimiter(path: str) -> str:
    """Guess a table's delimiter from its header line.

    The documented recipe for the sample map allows either tab- or comma-separated input,
    so the Hail reader has to agree with the Python parser about which it got.
    """
    with _open_read(path) as handle:
        header = handle.readline()
    return '\t' if '\t' in header else ','


def superpartition_column_table(path: str,
                                superpartition_size: int = DEFAULT_SUPERPARTITION_SIZE):
    """Read the sample -> superpartition mapping as a distributed Hail Table.

    ``hl.import_table`` rather than ``hl.Table.parallelize`` or ``hl.literal``: both of the
    latter inline every row into the query IR, and for a full-width VDS that is ~535K rows
    embedded in the plan. The driver then stalls before any real work begins, which reads
    as a job hung on its first stage rather than as an obviously bad query.

    Works for both input files: the sample map and the recorded sample list both carry
    ``sample_name`` and ``sample_id``, and the superpartition is recomputed from the id in
    either case rather than trusting a column that only one of them has.
    """
    _require_hail()
    table = hl.import_table(
        path, delimiter=sniff_delimiter(path), types={'sample_id': hl.tint64})
    table = table.key_by(s=table.sample_name)
    # Integer ceiling division, matching superpartition_for and the Avro export WDL.
    return table.select(
        sp=hl.int32((table.sample_id + superpartition_size - 1) // superpartition_size))


def annotated_matrix(vds, mode: str, mapping, bin_size: int,
                     injections: Sequence[Injection] = ()):
    """Return the mode-appropriate matrix annotated with superpartition and bin.

    Columns absent from the mapping are dropped rather than carried with a missing
    superpartition, which would otherwise collect into a null group in the aggregation.
    That is also what restricts a scan to the recorded sample list when it is pointed at a
    full-width VDS.

    Returns the matrix and how many columns found no match, so the caller can decide
    whether that is expected. It is for a recorded sample list, which is deliberately a
    subset, and it is not for a sample map, which should cover the whole callset.
    """
    _require_hail()
    matrix = vds.variant_data if mode == 'variants' else vds.reference_data

    matrix = matrix.annotate_cols(_sp=mapping[matrix.s].sp)
    with step('Summarizing columns by superpartition'):
        n_unmatched, counts = superpartition_column_summary(matrix)
    matrix = matrix.filter_cols(hl.is_defined(matrix._sp))
    # After the column filter, so `_sp` is defined on every remaining column, and before
    # binning, so an injection is invisible to everything downstream -- which is the point.
    matrix = apply_injections(matrix, injections)
    matrix = matrix.annotate_rows(
        _bin=(matrix.locus.position - 1) // bin_size,
    )
    return matrix, n_unmatched, counts


class Injection(NamedTuple):
    """A synthetic dropout to remove before summarizing.  Diagnostic use only.

    This exists to close the one link the rest of the validation cannot reach.  The
    detector's response to a dropout is measured by ``perturb_r2_summary.py``, which
    perturbs a summary; the fidelity of the Hail pass that produces that summary is
    established for references mode by an exact covered-base tie-out against BigQuery.
    Neither watches a hole in a VDS become a depleted count, and r2 offers no reference
    dropout to watch, its two being variant-only.

    Be precise about what an injection does and does not demonstrate.  Entries are filtered
    *after* the VDS is read, so this exercises the aggregation and everything downstream of
    it on real production data -- but it cannot show that ``read_vds`` would not materialize
    blocks that a truncated Avro never wrote.  That residual is closed separately, and
    empirically: the tie-out window is only 98.41% covered, so it contains real holes, and
    the VDS reproduces them rather than filling them.  Closing it *here* would mean writing
    a holed VDS and reading it back, which costs far more than it settles.
    """

    contig: str
    start: int
    end: int
    superpartition: int

    def __str__(self) -> str:
        return f'{self.contig}:{self.start}-{self.end}:{self.superpartition}'


# contig:start-end:superpartition.  Separators are digits-only on either side of the dash
# so a contig name containing one cannot be misread, and underscores and commas are allowed
# in positions because a 10-digit coordinate typed by hand is otherwise easy to get wrong.
INJECTION_SPEC = re.compile(
    r'^(?P<contig>\w+):(?P<start>[\d_,]+)-(?P<end>[\d_,]+):(?P<superpartition>\d+)$')


def parse_injection(spec: str) -> Injection:
    """Parse one ``contig:start-end:superpartition`` injection spec."""
    match = INJECTION_SPEC.match(spec.strip())
    if not match:
        raise ValueError(
            f'--inject-dropout {spec!r} is not contig:start-end:superpartition, '
            'e.g. chr20:1000000-1100000:83')
    start = int(match['start'].replace(',', '').replace('_', ''))
    end = int(match['end'].replace(',', '').replace('_', ''))
    if start > end:
        raise ValueError(f'--inject-dropout {spec!r}: start {start:,} is past end {end:,}')
    return Injection(match['contig'], start, end, int(match['superpartition']))


def parse_injection_list(values: Iterable[str] | None) -> tuple[Injection, ...]:
    """Parse repeated ``--inject-dropout`` values.

    Not comma-separated, unlike ``--contigs`` and the other list-valued options: a position
    may be written with comma digit separators, so a comma does not reliably delimit specs.
    Repeat the flag instead. The WDL exposes a single string and so carries one window.
    """
    if not values:
        return ()
    return tuple(parse_injection(value) for value in values)


def apply_injections(matrix, injections: Sequence[Injection]):
    """Filter out entries for a superpartition over a window, simulating a lost shard.

    Filters on ``locus.position``, the block's *start*, for both modes.  That is the same
    rule the metric uses -- a reference block's whole length is attributed to the bin its
    start falls in, with no clipping -- and the same rule a real shard loss follows, since
    the Avro export orders by location and a truncation cuts between whole rows.  Filtering
    on overlap instead would remove a block starting before the window, which no truncation
    does.

    ``filter_entries`` rather than ``filter_rows``: a dropout takes one superpartition's
    samples and leaves the rest of the callset intact at those loci, so the rows survive.
    Hail excludes filtered entries from ``hl.agg``, which is what makes them stand in for
    entries that were never written.
    """
    if not injections:
        # The overwhelmingly common path, and kept ahead of the Hail check so it stays
        # exercisable off-cluster.
        return matrix
    _require_hail()
    for injection in injections:
        inside = (
            (matrix.locus.contig == injection.contig)
            & (matrix.locus.position >= injection.start)
            & (matrix.locus.position <= injection.end)
            & (matrix._sp == injection.superpartition)
        )
        matrix = matrix.filter_entries(~inside)
    return matrix


def reference_length_field(entry_fields) -> str:
    """Which spelling of reference-block length a VDS uses: ``'LEN'`` or ``'END'``.

    Forward compatibility rather than a live concern.  ``read_vds`` in 0.2.134 calls
    ``_add_len`` and then ``_add_end``, so both fields are present in memory whatever the VDS
    holds on disk, and either would work today.  ``LEN`` is preferred because that is the
    direction Hail is moving -- ``VariantDataset.write`` already stores ``LEN`` and drops
    ``END``, for VCF 4.5 alignment and better compression, and ``read_vds`` carries a
    ``_drop_end`` flag -- so a future version may stop synthesizing ``END``.

    Takes anything iterable over field names -- ``matrix.entry.dtype`` in production -- so
    the choice is testable without Hail, which the surrounding expression building is not.
    """
    fields = set(entry_fields)
    if 'LEN' in fields:
        return 'LEN'
    if 'END' in fields:
        return 'END'
    raise ValueError(
        'reference_data entries carry neither LEN nor END, so reference-block length '
        f'cannot be computed; found {sorted(fields)}'
    )


def _metric_expression(matrix, mode: str):
    """The per-entry quantity being summed for the given mode."""
    _require_hail()
    if mode == 'variants':
        # LA is defined for every present entry. LGT is deliberately not used: a variant
        # call with call_GQ = 0 is imported with a missing LGT but a present entry, and
        # it corresponds to a real vet row, so it must be counted.
        return hl.agg.count_where(hl.is_defined(matrix.LA))
    if reference_length_field(matrix.entry.dtype) == 'LEN':
        length = matrix.LEN
    else:
        length = matrix.END - matrix.locus.position + 1
    return hl.agg.sum(hl.if_else(hl.is_defined(length), length, 0))


def executor_summary() -> str:
    """Live executor count and total task slots, for the log.

    Worth recording because it is the term most likely to be wrong in a cost estimate and
    the hardest to reconstruct afterwards -- a cluster that starts with two primary workers
    and autoscales for the rest can spend a long time far below peak width, and an estimate
    assuming full width is then silently wrong by that ratio.

    Each figure is obtained independently, in its own `try`, so that diagnostics degrade one
    field at a time. Taking them together is easy to get wrong: py4j surfaces
    `getExecutorMemoryStatus().keySet()` as a Java object rather than a Python iterable, and
    a TypeError from iterating it would discard the task-slot count along with it.
    """
    _require_hail()
    try:
        context = hl.spark_context()
    except Exception as exc:  # pragma: no cover - needs a live Spark context
        return f'unavailable (no Spark context: {exc})'

    parts = []
    try:
        # A Scala Map, so ask it for its size rather than iterating it. The driver is
        # included in this map but runs no tasks, hence the -1.
        n_executors = int(context._jsc.sc().getExecutorMemoryStatus().size()) - 1
        parts.append(f'{max(n_executors, 0)} executor(s)')
    except Exception as exc:  # pragma: no cover - needs a live Spark context
        parts.append(f'executor count unavailable ({type(exc).__name__})')
    try:
        parts.append(f'~{context.defaultParallelism} task slot(s)')
    except Exception as exc:  # pragma: no cover - needs a live Spark context
        parts.append(f'task slots unavailable ({type(exc).__name__})')
    return ', '.join(parts)


@contextlib.contextmanager
def width_heartbeat(label: str, interval_seconds: int = WIDTH_HEARTBEAT_SECONDS):
    """Log cluster width periodically while a long aggregation runs.

    Endpoint readings describe an autoscaling cluster poorly. One taken before the
    aggregation starts is taken before autoscaling has seen the work, so it reflects
    whatever the previous unit of work left behind and is systematically low; one taken at
    the end misses a slow ramp entirely. Cost is per-partition time divided by concurrent
    tasks, so the term being divided by is a profile, not a number -- and reading a
    pre-ramp snapshot as though it were characteristic is an easy way to produce a cost
    estimate that is wrong by a large factor.

    A daemon thread, and every sample is wrapped, so nothing here can fail or delay the run
    it is describing.
    """
    stop = threading.Event()

    def sample() -> None:
        started = time.monotonic()
        while not stop.wait(interval_seconds):
            try:
                elapsed = (time.monotonic() - started) / 60
                announce(f'{label}: {elapsed:,.1f} min elapsed, {executor_summary()}')
            except Exception:  # pragma: no cover - never worth failing a run over
                pass

    thread = threading.Thread(target=sample, daemon=True, name='width-heartbeat')
    thread.start()
    try:
        yield
    finally:
        stop.set()
        thread.join(timeout=5)


def matrix_partition_count(matrix) -> int:
    """Partitions the aggregation will run over, i.e. its maximum parallelism.

    Worth reporting prominently. A narrow interval on an AoU-scale VDS can resolve to a
    single partition, and a single partition means a single task: the whole cluster idles
    while one thread streams tens of gigabytes of entry data. That is a property of the
    interval, not of the query, and it is invisible without asking.
    """
    _require_hail()
    return int(matrix.n_partitions())


def superpartition_column_summary(matrix) -> tuple[int, dict[int, int]]:
    """Unmatched column count and per-superpartition sample counts, in one pass.

    A single ``aggregate_cols`` rather than two. Reading the column table of a full-width
    AoU VDS is the one operation observed to be pathologically slow on these callsets, so
    it is done once and both numbers are taken from the same pass.

    Counted from the matrix rather than from the mapping file so the figures reflect the
    intersection actually screened -- a recorded sample list can name samples a later VDS
    no longer has.
    """
    _require_hail()
    result = matrix.aggregate_cols(hl.struct(
        unmatched=hl.agg.count_where(hl.is_missing(matrix._sp)),
        counts=hl.agg.filter(hl.is_defined(matrix._sp), hl.agg.counter(matrix._sp)),
    ))
    counts = {int(sp): int(n) for sp, n in result.counts.items() if sp is not None}
    return int(result.unmatched), counts


def aggregate_totals(matrix, mode: str,
                     bin_size: int) -> dict[tuple[str, int], dict[int, float]]:
    """Aggregate a VDS into per-(bin, superpartition) totals.

    Expressed as one nested ``hl.agg.group_by`` inside a single ``aggregate_entries`` so
    the whole thing is a streaming pass with no shuffle: the outer grouping folds bins and
    the inner folds superpartitions, and only the finished matrix -- on the order of
    25,000 x 134 numbers for the largest contig at the 10 kb default -- returns to the
    driver.  Each
    partition's accumulator stays small because a partition spans few bins.

    Size scales with bins, not with the genome: at the 10 kb default the largest contig
    yields about 25,000 x 134 numbers, and a single-shot genome-wide run about 310,000 x
    134.  Both are small against the master this runs on.
    """
    _require_hail()
    raw = matrix.aggregate_entries(
        hl.agg.group_by(
            hl.tuple([matrix.locus.contig, matrix._bin]),
            hl.agg.group_by(matrix._sp, _metric_expression(matrix, mode)),
        )
    )

    totals: dict[tuple[str, int], dict[int, float]] = {}
    for (contig, bin_index), per_superpartition in raw.items():
        key = (contig, bin_index * bin_size + 1)
        totals[key] = {int(sp): float(value) for sp, value in per_superpartition.items()}
    return totals


def per_superpartition_sample_counts(vds, mode: str, mapping,
                                     target: set[int] | None = None,
                                     ) -> dict[int, dict[str, float]]:
    """Per-sample totals within the VDS's current interval, grouped by superpartition.

    Column filtering happens in Hail, so only the requested superpartitions are ever
    brought back to the driver -- a few thousand samples rather than the whole callset.
    """
    _require_hail()
    matrix = vds.variant_data if mode == 'variants' else vds.reference_data
    matrix = matrix.annotate_cols(_sp=mapping[matrix.s].sp)
    matrix = matrix.filter_cols(hl.is_defined(matrix._sp))
    if target:
        matrix = matrix.filter_cols(hl.literal(sorted(target)).contains(matrix._sp))
    matrix = matrix.annotate_cols(_observed=_metric_expression(matrix, mode))

    grouped: dict[int, dict[str, float]] = defaultdict(dict)
    for row in matrix.cols().select('_sp', '_observed').collect():
        grouped[int(row._sp)][row.s] = float(row._observed)
    return dict(sorted(grouped.items()))


# ---------------------------------------------------------------------------
# Actions
# ---------------------------------------------------------------------------


def _screening_matrix(args, vds):
    """Annotated matrix and the per-superpartition sample counts it will be judged on."""
    source = args.sample_map_path
    with step(f'Reading sample mapping from {source}'):
        mapping = superpartition_column_table(source, args.superpartition_size)
    matrix, n_unmatched, counts = annotated_matrix(
        vds, args.mode, mapping, args.bin_size, args.injections)

    for injection in args.injections:
        announce(f'*** INJECTED DROPOUT {injection} -- diagnostic run, output is NOT a '
                 'valid scan of this VDS ***')

    if n_unmatched:
        raise ValueError(
            f'{n_unmatched:,} VDS sample(s) are absent from {source}. The map must cover '
            'every sample in the VDS: screening part of a superpartition would bias the '
            'peer comparison the detector relies on, and would do so silently. Regenerate '
            'the map from sample_info.')

    total = sum(counts.values())
    announce(f'Screening {total:,} samples across {len(counts):,} superpartitions')

    with step('Counting partitions'):
        n_partitions = matrix_partition_count(matrix)
    with step('Counting partitions in the whole VDS'):
        n_total = total_partition_count(args.vds_path, args.mode)
    announce(f'Aggregation will run over {n_partitions:,} of {n_total:,} total partition(s)')

    if n_partitions < MIN_HEALTHY_PARTITIONS:
        announce(
            f'WARNING: only {n_partitions} partition(s) selected, so at most '
            f'{n_partitions} task(s) will run and the rest of the cluster will idle. A '
            'partition is the unit of work -- one task streams one end to end -- and on an '
            'AoU-scale VDS a single partition can hold tens of GB of entry data.')
        announce(
            f'         Because there is no parallelism to divide by, this will take roughly '
            f'what a fully parallel scan of all {n_total:,} partitions would. It is not a '
            'cheap preview. Either widen the interval to span many partitions, or run the '
            'real scan and let the cluster work.')

    return matrix, counts


def shard_paths(summary_path: str, contig: str) -> tuple[str, str]:
    """Per-contig summary shard and its completion marker.

    A separate marker rather than the shard's own existence, because a write interrupted
    part-way leaves a shard that looks finished. The marker is written only after the shard
    closes cleanly, so a resume cannot mistake a truncated shard for a complete one -- which
    would silently drop a contig and report the remainder as clean, the exact class of
    omission this tool exists to find.
    """
    return f'{summary_path}.{contig}', f'{summary_path}.{contig}.done'


# The marker records what produced the shard, not just that something did. Without this, a
# re-run with the same output_prefix but a different --vds-path or --bin-size would silently
# reuse the previous run's shards and emit a summary describing the wrong VDS -- a silent
# wrong answer, which is the single outcome this tool cannot afford to produce.
MARKER_HEADER = 'contig\trows\tvds_path\tmode\tbin_size\tinjections'


def injection_provenance(injections: Sequence[Injection]) -> str:
    """How a run's injections are recorded in a shard marker.

    A shard produced with an injected dropout must never be picked up by a later clean
    resume: it holds a hole that no VDS has, and the summary built from it would look
    complete while describing data that never existed. Recording the spec makes that
    collision an abort rather than a silent wrong answer, which is the same reason
    `vds_path`, `mode` and `bin_size` are recorded.
    """
    return ','.join(str(injection) for injection in injections) if injections else '-'


def write_marker(marker_path: str, contig: str, n_rows: int, args) -> None:
    """Record the shard's provenance alongside its row count."""
    write_lines(marker_path, MARKER_HEADER,
                [f'{contig}\t{n_rows}\t{args.vds_path}\t{args.mode}\t{args.bin_size}\t'
                 f'{injection_provenance(args.injections)}'])


def verify_marker(marker_path: str, contig: str, args) -> None:
    """Check a checkpoint was produced by the run we are resuming.

    A mismatch aborts rather than warning: reusing another VDS's shard produces a summary
    that looks complete and describes the wrong data, and no downstream step could detect
    it. An unrecognized header is treated the same way, since a marker this cannot read is
    a marker whose provenance it cannot check.
    """
    with _open_read(marker_path) as handle:
        header = handle.readline().rstrip('\n')
        row = handle.readline().rstrip('\n')

    if header != MARKER_HEADER:
        raise RuntimeError(f'{marker_path}: unrecognized marker header {header!r}, so its '
                           f'provenance cannot be checked. Delete it and the matching '
                           f'shard, or choose a different --summary-path.')

    fields = row.split('\t')
    if len(fields) < 6:
        raise RuntimeError(f'{marker_path}: malformed marker row {row!r}; delete it and '
                           f're-run {contig}.')
    _, _, vds_path, mode, bin_size, injections = fields[:6]
    mismatches = []
    if injections != injection_provenance(args.injections):
        mismatches.append(
            f'injections {injections!r} != {injection_provenance(args.injections)!r}')
    if vds_path != args.vds_path:
        mismatches.append(f'vds_path {vds_path!r} != {args.vds_path!r}')
    if mode != args.mode:
        mismatches.append(f'mode {mode!r} != {args.mode!r}')
    if bin_size != str(args.bin_size):
        mismatches.append(f'bin_size {bin_size} != {args.bin_size}')
    if mismatches:
        raise RuntimeError(
            f'{marker_path} was produced by a different run ({"; ".join(mismatches)}). '
            'Reusing it would emit a summary describing the wrong data while appearing '
            f'complete. Delete {marker_path} and the matching shard, or choose a different '
            '--summary-path.')


# Bytes per `read` during the merge, which is a pure copy of a few hundred megabytes. Large
# because the merge moves undecoded bytes; see `_binary_openers` for why it does not go
# through Hail.
MERGE_CHUNK_BYTES = 32 * 1024 * 1024

# Only the fallback below still reads through Hail's text handles, and it is still capped by
# the 8 KB buffer `_open_read` describes, so a large value here buys it nothing.
MERGE_CHUNK_CHARACTERS = 8 * 1024 * 1024

# How often the copy reports progress. It is the longest non-Hail step in the job, and a step
# with no Hail UI behind it is exactly the one whose progress has to come from the script.
MERGE_PROGRESS_SECONDS = 60


def _gcs_bucket_and_name(path: str) -> tuple[str, str]:
    """Split a ``gs://bucket/object`` path into its two parts."""
    bucket, _, name = path[len('gs://'):].partition('/')
    if not bucket or not name:
        raise ValueError(f'{path!r} is not a gs://bucket/object path')
    return bucket, name


def _binary_openers(paths: Sequence[str]):
    """Binary open functions for these paths, or ``None`` if GCS is needed and unavailable.

    The merge is the only place in this script that moves bulk bytes, and it must not do so
    through `hl.hadoop_open`. That handle is capped at an 8 KB JVM read -- see `_open_read`
    for why the cap cannot be lifted -- which put the Foxtrot references merge at ~180 kB/s,
    just under two hours to copy 1.3 GB. That merge, not the ten-hour Hail scan ahead of it,
    is what ran the job past its cluster TTL and discarded the result.

    Local paths use the builtins, and are answered without touching GCS at all, so a run
    over local files neither imports `google.cloud.storage` nor constructs a client that
    would go looking for credentials. Unit tests therefore exercise this same merge code.

    `google.cloud.storage` is imported lazily rather than at module scope so that its
    absence degrades to the slow path instead of failing a job that has already done all
    the expensive work. `run_in_hail_cluster.py` asks `hailctl dataproc start` to install it.
    """
    if not any(path.startswith('gs://') for path in paths):
        return (lambda path: open(path, 'rb')), (lambda path: open(path, 'wb'))

    try:
        from google.cloud import storage
    except ImportError:
        return None

    client = storage.Client()

    def open_read(path: str) -> IO[bytes]:
        if not path.startswith('gs://'):
            return open(path, 'rb')
        bucket, name = _gcs_bucket_and_name(path)
        return client.bucket(bucket).blob(name).open('rb')

    def open_write(path: str) -> IO[bytes]:
        if not path.startswith('gs://'):
            return open(path, 'wb')
        bucket, name = _gcs_bucket_and_name(path)
        return client.bucket(bucket).blob(name).open('wb')

    return open_read, open_write


def _verify_shard_coverage(seen: set[str], expected_contigs: Sequence[str],
                           summary_path: str) -> None:
    """Raise unless every requested contig contributed rows to the merged summary."""
    missing = [c for c in expected_contigs if c not in seen]
    if missing:
        raise RuntimeError(
            f'The merged summary has no rows for {missing}. Every requested contig must be '
            'represented or the screen is incomplete while appearing clean. Delete the '
            f'stale .done markers matching {summary_path}.* and re-run those contigs.')


def concatenate_shards(shards: Sequence[str], summary_path: str,
                       expected_contigs: Sequence[str]) -> int:
    """Merge per-contig shards into the final summary, verifying coverage.

    The verification is the point: a resume that skipped a contig would otherwise produce a
    summary that looks complete. Introducing a silent omission by way of the resume logic
    would be a poor way to run a tool built to detect silent omissions.

    Copied as raw bytes in large chunks. The only per-row work is recording which contigs
    are present, and a shard holds exactly one contig, so its first data row answers that
    for the whole shard; everything after it is moved without being decoded or examined.

    That makes the coverage check cheap but not weaker in the direction that matters: a
    hypothetical shard holding rows for a contig other than its first would now be reported
    missing rather than passing. Erring toward raising is the right side to fail on here.

    The returned count is of newlines, which is the same as rows for any shard `write_lines`
    produced. A blank line in a hand-edited one would be counted, and `vds_dropout_detect.py`
    skips blank lines when it reads the summary back.
    """
    openers = _binary_openers([summary_path, *shards])
    if openers is None:
        announce(
            'WARNING: google-cloud-storage is not importable on this cluster, so the merge '
            "falls back to Hail's text handles. That path ran at ~180 kB/s on the Foxtrot "
            'references summary -- hours rather than minutes. Install google-cloud-storage '
            'to avoid it.')
        return _concatenate_shards_via_hail(shards, summary_path, expected_contigs)
    open_read, open_write = openers

    header = SUMMARY_HEADER.encode()
    seen: set[str] = set()
    total = 0
    started = time.monotonic()
    reported = started
    with open_write(summary_path) as out:
        out.write(header + b'\n')
        for index, shard in enumerate(shards, 1):
            with open_read(shard) as handle:
                first = handle.readline()
                if first.strip() != header:
                    found = first.strip().decode('utf-8', 'replace')
                    raise ValueError(
                        f'{shard}: expected header {SUMMARY_HEADER!r}, found {found!r}; '
                        f'the shard is corrupt or from another run')

                # One row read by hand, because its contig identifies the shard. The rest
                # is copied without being looked at.
                row = handle.readline()
                while row and not row.strip():
                    row = handle.readline()
                if not row:
                    # Header only. Its contig stays unseen, and is reported below if it was
                    # one of the expected ones.
                    announce(f'{shard} holds no data rows ({index} of {len(shards)})')
                    continue
                seen.add(row.split(b'\t', 1)[0].decode('utf-8', 'replace'))
                out.write(row if row.endswith(b'\n') else row + b'\n')
                total += 1

                last = b'\n'
                while True:
                    chunk = handle.read(MERGE_CHUNK_BYTES)
                    if not chunk:
                        break
                    out.write(chunk)
                    total += chunk.count(b'\n')
                    last = chunk[-1:]
                    now = time.monotonic()
                    if now - reported >= MERGE_PROGRESS_SECONDS:
                        rate = total / max(now - started, 1e-9)
                        announce(f'{shard}: {total:,} rows merged so far, {rate:,.0f} '
                                 f'rows/s ({index} of {len(shards)} shards)')
                        reported = now
                if last != b'\n':
                    # Without this a shard whose final row has no newline would glue that
                    # row onto the first row of the next shard, losing both.
                    out.write(b'\n')
                    total += 1
            announce(f'Merged {shard} ({index} of {len(shards)}): {total:,} rows so far')

    _verify_shard_coverage(seen, expected_contigs, summary_path)
    return total


def _concatenate_shards_via_hail(shards: Sequence[str], summary_path: str,
                                 expected_contigs: Sequence[str]) -> int:
    """The pre-`google-cloud-storage` merge, kept as a fallback.

    Identical in behavior to `concatenate_shards` and roughly two orders of magnitude
    slower, because every read and write crosses the Python/JVM boundary 8 KB at a time
    under `hl.hadoop_open`. Reached only when `google.cloud.storage` cannot be imported.
    """
    seen: set[str] = set()
    total = 0
    started = time.monotonic()
    reported = started
    with _open_write(summary_path) as out:
        out.write(SUMMARY_HEADER + '\n')
        for index, shard in enumerate(shards, 1):
            with _open_read(shard) as handle:
                first = handle.readline()
                if first.strip() != SUMMARY_HEADER:
                    raise ValueError(
                        f'{shard}: expected header {SUMMARY_HEADER!r}, found '
                        f'{first.strip()!r}; the shard is corrupt or from another run')

                row = handle.readline()
                while row and not row.strip():
                    row = handle.readline()
                if not row:
                    announce(f'{shard} holds no data rows ({index} of {len(shards)})')
                    continue
                seen.add(row.split('\t', 1)[0])
                out.write(row if row.endswith('\n') else row + '\n')
                total += 1

                last = '\n'
                while True:
                    chunk = handle.read(MERGE_CHUNK_CHARACTERS)
                    if not chunk:
                        break
                    out.write(chunk)
                    total += chunk.count('\n')
                    last = chunk[-1]
                    now = time.monotonic()
                    if now - reported >= MERGE_PROGRESS_SECONDS:
                        rate = total / max(now - started, 1e-9)
                        announce(f'{shard}: {total:,} rows merged so far, {rate:,.0f} '
                                 f'rows/s ({index} of {len(shards)} shards)')
                        reported = now
                if last != '\n':
                    out.write('\n')
                    total += 1
            announce(f'Merged {shard} ({index} of {len(shards)}): {total:,} rows so far')

    _verify_shard_coverage(seen, expected_contigs, summary_path)
    return total


def _scan_one_contig(args, contig: str):
    """Aggregate one contig into its shard. Returns (shard path, sample counts or None)."""
    shard, marker = shard_paths(args.summary_path, contig)
    if _path_exists(marker):
        verify_marker(marker, contig, args)
        announce(f'{contig}: already complete, skipping')
        return shard, None

    intervals = build_intervals([contig], [])
    with step(f'{contig}: opening VDS'):
        vds = read_and_subset_vds(args.vds_path, intervals)
    matrix, counts = _screening_matrix(args, vds)

    # Explicitly labeled: autoscaling has not seen this contig's work yet, so this is a
    # floor rather than the width the aggregation will actually run at.
    announce(f'{contig}: cluster width before autoscaling responds: {executor_summary()}')
    with step(f'{contig}: aggregating {args.mode} into {args.bin_size:,} bp bins'), \
            width_heartbeat(f'{contig}: aggregating'):
        totals = aggregate_totals(matrix, args.mode, args.bin_size)
    announce(f'{contig}: cluster width at completion: {executor_summary()}')

    n_rows = write_lines(shard, SUMMARY_HEADER,
                         format_summary_rows(totals, args.bin_size))
    write_marker(marker, contig, n_rows, args)
    announce(f'{contig}: wrote {n_rows:,} rows ({len(totals):,} bins)')
    return shard, counts


def _scan_intervals(args, requested_intervals: Sequence[str]) -> int:
    """Single-pass scan over an explicit interval list, without checkpointing.

    A targeted interval run is short enough not to need resuming, and has no natural
    per-contig sharding.
    """
    intervals = build_intervals(parse_contig_list(args.contigs), requested_intervals)
    announce(f'Scanning {len(intervals):,} interval(s) of {args.vds_path}')
    with step(f'Opening VDS {args.vds_path}'):
        vds = read_and_subset_vds(args.vds_path, intervals)
    matrix, counts = _screening_matrix(args, vds)

    announce(f'Cluster width before autoscaling responds: {executor_summary()}')
    with step(f'Aggregating {args.mode} into {args.bin_size:,} bp bins'), \
            width_heartbeat('Aggregating'):
        totals = aggregate_totals(matrix, args.mode, args.bin_size)
    announce(f'Cluster width at completion: {executor_summary()}')

    n_rows = write_lines(args.summary_path, SUMMARY_HEADER,
                         format_summary_rows(totals, args.bin_size))
    n_sp = write_lines(
        args.superpartitions_path, SUPERPARTITION_HEADER,
        [f'{sp}\t{counts[sp]}' for sp in sorted(counts)])
    announce(f'Wrote {n_rows:,} summary rows ({len(totals):,} bins) to {args.summary_path}')
    announce(f'Wrote {n_sp:,} superpartition rows to {args.superpartitions_path}')
    return 0


def action_scan(args) -> int:
    """Aggregate a VDS and write the summary and superpartition tables.

    Contigs are aggregated and written one at a time, each with a completion marker, so an
    interrupted run resumes instead of restarting. A genome-wide scan of an AoU-scale VDS
    runs for hours, and transient failures -- a capacity stockout, a metadata blip, a
    preempted driver -- are frequent enough that discarding all of it on the last one is not
    acceptable.
    """
    _require_hail()
    requested_intervals = parse_interval_list(args.intervals)
    if requested_intervals:
        return _scan_intervals(args, requested_intervals)

    contigs = parse_contig_list(args.contigs)
    announce(f'Scanning {len(contigs)} contig(s) of {args.vds_path}, checkpointing each')

    shards = []
    counts = None
    for contig in contigs:
        shard, contig_counts = _scan_one_contig(args, contig)
        shards.append(shard)
        if contig_counts is not None:
            counts = contig_counts

    if counts is None:
        # Every contig was already done, so nothing re-read the column table. The
        # superpartition table from the original run is still correct; leave it be.
        announce('All contigs were already complete; leaving the superpartition table as is')
    else:
        n_sp = write_lines(
            args.superpartitions_path, SUPERPARTITION_HEADER,
            [f'{sp}\t{counts[sp]}' for sp in sorted(counts)])
        announce(f'Wrote {n_sp:,} superpartition rows to {args.superpartitions_path}')

    with step('Merging per-contig shards'):
        n_rows = concatenate_shards(shards, args.summary_path, contigs)
    announce(f'Wrote {n_rows:,} summary rows to {args.summary_path}')
    return 0


def action_full_depth(args) -> int:
    """Count per-sample data in flagged windows, at full width."""
    _require_hail()
    requested_intervals = parse_interval_list(args.intervals)
    if not requested_intervals:
        raise ValueError('--intervals is required for --action full-depth')
    target = {int(s) for s in args.target_superpartitions.split(',')} \
        if args.target_superpartitions else None

    with step(f'Reading sample mapping from {args.sample_map_path}'):
        mapping = superpartition_column_table(args.sample_map_path, args.superpartition_size)

    rows: list[str] = []
    for interval_string in requested_intervals:
        intervals = build_intervals([], [interval_string])
        with step(f'Opening VDS for {interval_string}'):
            vds = read_and_subset_vds(args.vds_path, intervals)

        contig = intervals[0].start.contig
        start = intervals[0].start.position
        end = intervals[0].end.position

        # Columns are narrowed in Hail before anything is brought back, so the driver only
        # ever holds one superpartition's worth of samples rather than the whole callset.
        with step(f'Counting per-sample data in {interval_string}'):
            per_superpartition = per_superpartition_sample_counts(
                vds, args.mode, mapping, target)
        for superpartition, observed in per_superpartition.items():
            n_zero = sum(1 for value in observed.values() if value <= 0)
            total = sum(observed.values())
            names = observed
            fraction = n_zero / len(names) if names else 0.0
            rows.append(
                f'{contig}\t{start}\t{end}\t{superpartition}\t{len(names)}\t{n_zero}\t'
                f'{fraction:.6f}\t{total:.0f}'
            )
            print(
                f'{contig}:{start:,}-{end:,} sp {superpartition}: {n_zero:,} of '
                f'{len(names):,} samples carry no data ({fraction * 100:.2f}%), '
                f'{total:,.0f} total'
            )

    n_rows = write_lines(args.full_depth_path, FULL_DEPTH_HEADER, rows)
    print(f'Wrote {n_rows:,} rows to {args.full_depth_path}')
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser.

    Flags are kebab-case throughout because `run_in_hail_cluster.py` turns the keys of
    its arguments JSON straight into ``--key value`` pairs.
    """
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='Scan a GVS Hail VDS for rectangle data dropouts.',
    )
    parser.add_argument('--action', choices=ACTIONS, required=True,
                        help='Which step to run. See the module docstring.')
    parser.add_argument('--vds-path', required=True,
                        help='VDS to read.')
    parser.add_argument('--mode', choices=MODES, default='variants',
                        help='Metric to compute. Default: variants.')

    parser.add_argument('--sample-map-path', required=True,
                        help='TSV of sample_name and sample_id covering every sample in the '
                             'VDS. Superpartition membership is a function of sample_id, '
                             'which a VDS does not carry.')
    parser.add_argument('--superpartitions-path', default=None,
                        help='Where to write the superpartition/n_samples table.')
    parser.add_argument('--summary-path', default=None,
                        help='Where to write the summary table (scan).')
    parser.add_argument('--full-depth-path', default=None,
                        help='Where to write per-sample results (full-depth).')

    parser.add_argument('--superpartition-size', type=int, default=DEFAULT_SUPERPARTITION_SIZE,
                        help=f'Samples per GVS superpartition. Default: '
                             f'{DEFAULT_SUPERPARTITION_SIZE}.')
    parser.add_argument('--bin-size', type=int, default=DEFAULT_BIN_SIZE,
                        help=f'Genomic bin size in bp. Default: {DEFAULT_BIN_SIZE}.')

    parser.add_argument('--contigs', default=None,
                        help='Comma-separated contigs to scan. Default: all GVS primary contigs.')
    parser.add_argument('--intervals', action='append', default=None,
                        help='Hail locus interval. Repeatable, and each value may itself be a '
                             'comma-separated list. Overrides --contigs.')
    parser.add_argument('--target-superpartitions', default=None,
                        help='Comma-separated superpartitions to restrict full-depth to.')

    parser.add_argument('--inject-dropout', action='append', default=None,
                        metavar='CONTIG:START-END:SUPERPARTITION',
                        help='DIAGNOSTIC ONLY. Remove one superpartition\'s entries over a '
                             'window before summarizing, so that a known dropout can be '
                             'watched through the scan and the detector end to end. '
                             'Repeatable. A run using this does not describe the VDS it '
                             'reads, and its shards are marked so a later clean run '
                             'refuses them.')

    parser.add_argument('--temp-path', default=None, help='Hail temporary directory.')
    return parser


REQUIRED_BY_ACTION: dict[str, tuple[str, ...]] = {
    'scan': ('summary_path', 'superpartitions_path'),
    'full-depth': ('full_depth_path',),
}


def validate_args(args, parser: argparse.ArgumentParser) -> None:
    """Fail fast on missing per-action arguments, before any cluster time is spent."""
    for name in REQUIRED_BY_ACTION[args.action]:
        if getattr(args, name, None) is None:
            parser.error(f"--{name.replace('_', '-')} is required for --action {args.action}")
    if args.bin_size < 1:
        parser.error('--bin-size must be at least 1')
    try:
        args.injections = parse_injection_list(args.inject_dropout)
    except ValueError as exc:
        parser.error(str(exc))


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    _require_hail()
    hl.init(idempotent=True, tmp_dir=args.temp_path)
    hl.default_reference('GRCh38')

    handler = {
        'scan': action_scan,
        'full-depth': action_full_depth,
    }[args.action]
    return handler(args)


if __name__ == '__main__':
    sys.exit(main())
