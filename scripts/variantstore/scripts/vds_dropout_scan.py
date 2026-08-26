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

``--action probe``
    Run ``scan`` over a bounded interval and report measured throughput plus a genome-wide
    projection.  Useful for sizing a new callset.

    A partition is the unit of work -- one task streams one partition end to end -- so the
    interval must span enough partitions to exercise the cluster.  An interval narrow
    enough to resolve to a handful of partitions measures near-serial throughput and takes
    roughly as long as a fully parallel scan of the whole genome, which is the opposite of
    a cheap preview.  The report says so when it happens.

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
    Sums reference-block coverage, ``END - position + 1``.  Each block is attributed
    entirely to the bin holding its start locus; with 50 kb bins and a 1000 bp maximum
    block length the edge error is under 2%, and it is a uniform bias that cancels in
    the superpartition-versus-peers comparison the detector performs.

    Reference dropouts leave blocks *absent* rather than GQ 0.  AoU ingests with
    ``drop_state = "ZERO"`` so GQ0 blocks never reach BigQuery, and nothing in VDS
    construction back-fills gaps -- `import_gvs.py` builds `reference_data` purely from
    the Avro rows, and `update_reference_data_ploidy` only annotates `GT` on entries that
    already exist.

Every sample is screened
------------------------
There is no sampling, of samples or of loci.  An earlier design screened a stratified
subset of ~100 samples per superpartition, on the premise that a full-width pass was
expensive enough to need amortizing over a downsampled copy.  Measurement retired that
premise: with the reader pruning to native partitions, a genome-wide variant scan of a
535K-sample VDS runs in about an hour and a reference scan in a couple of hours.  Screening
everything is therefore simpler, strictly more sensitive -- a sampled screen cannot see a
handful of individually lost samples -- and needs none of the machinery that made sampling
safe to compare across VDSes.

Sample map
----------
Superpartition membership comes from a TSV of sample name to sample ID, since a VDS knows
only names.  ``GvsValidateVdsCompleteness.wdl`` generates that map itself when given a
BigQuery project and dataset; ``scripts/variantstore/bq/vds_dropout_sample_map.sql`` is
the hand-run copy, for building one outside the workflow.  Either tab- or comma-separated
input is accepted, so ``bq query --format=csv`` output works as-is.

The map must cover every sample in the VDS.  A VDS sample missing from it is a hard error
rather than a skip, because screening part of a superpartition biases the peer comparison
the detector relies on, and would do so silently.  Samples in the map but absent from the
VDS are ignored, so a map generated after a withdrawal is still usable.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import IO, Iterable, Mapping, Sequence

try:
    import hail as hl
except ModuleNotFoundError:  # pragma: no cover - exercised only off-cluster
    hl = None

DEFAULT_SUPERPARTITION_SIZE = 4000
DEFAULT_BIN_SIZE = 50_000
# Below this many partitions the aggregation cannot use the cluster, which on an AoU-scale
# VDS turns a nominally small interval into hours of single-threaded streaming.
MIN_HEALTHY_PARTITIONS = 8
# A 10 Mb window on chr20 that is 100% callable: comfortably clear of the centromere gap
# (chr20:26,386,232-30,088,349) and of the telomere, and free of the smaller assembly gaps
# near 31 Mb. Verified against
# wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list. An interval overlapping
# a gap would make the probe's timing optimistic and emit confusing empty bins.
DEFAULT_PROBE_INTERVAL = 'chr20:37000000-47000000'
# End of the chr20 centromere gap, from the same interval list. Used to keep the default honest.
CHR20_CENTROMERE_END = 30_088_349

MODES = ('variants', 'references')
ACTIONS = ('scan', 'probe', 'full-depth')

# GRCh38 primary contigs GVS encodes. Alt and decoy contigs carry no GVS location
# encoding and are not screened.
GVS_CONTIGS = tuple([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])

# Approximate GRCh38 primary assembly length, used only to extrapolate probe timings.
GENOME_LENGTH = 3_088_269_832

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


def extrapolate_runtime(
        elapsed_seconds: float,
        scanned_bases: int,
        genome_length: int = GENOME_LENGTH,
) -> float:
    """Scale a probe's wall time up by genomic span. Kept for reference only.

    Misleading on its own, because it assumes the probe used the cluster as well as a
    genome-wide run would. A narrow interval resolves to few partitions and therefore few
    tasks, so it measures near-serial throughput and overstates the real cost, often by
    the width of the cluster. Prefer ``extrapolate_by_partitions``.
    """
    if scanned_bases <= 0:
        raise ValueError('scanned_bases must be positive')
    return elapsed_seconds * (genome_length / scanned_bases)


def extrapolate_by_partitions(elapsed_seconds: float, n_partitions_probed: int,
                              n_partitions_total: int, concurrent_tasks: int) -> float:
    """Estimate a full run from per-partition cost, in seconds.

    The unit of work is a partition, not a base pair: one task streams one partition end
    to end. So the honest model is time-per-partition, multiplied by the partitions a full
    run covers, divided by how many run at once -- which is what makes a genome-wide scan
    cheaper per partition than a narrow probe rather than more expensive.
    """
    for name, value in (('n_partitions_probed', n_partitions_probed),
                        ('n_partitions_total', n_partitions_total),
                        ('concurrent_tasks', concurrent_tasks)):
        if value <= 0:
            raise ValueError(f'{name} must be positive, got {value}')
    per_partition = elapsed_seconds / n_partitions_probed
    return per_partition * n_partitions_total / concurrent_tasks


def format_probe_report(
        elapsed_seconds: float,
        scanned_bases: int,
        n_cells: int,
        n_bins: int,
        n_samples: int | None = None,
        n_partitions: int | None = None,
        n_total_partitions: int | None = None,
        parallelism: int | None = None,
        full_parallelism: int | None = None,
) -> str:
    """Human-readable probe result.

    Reports per-partition cost rather than a bp-scaled figure, because the partition is
    the unit of work and a narrow probe is close to serial.
    """
    lines = [
        'Probe result',
        f'  interval span:        {scanned_bases:,} bp',
        f'  bins produced:        {n_bins:,}',
        f'  non-zero cells:       {n_cells:,}',
    ]
    if n_samples is not None:
        lines.append(f'  samples screened:     {n_samples:,}')
    if n_partitions:
        lines.append(f'  partitions read:      {n_partitions:,}')
    lines.append(f'  elapsed:              {elapsed_seconds:,.1f} s')

    if n_partitions:
        per_partition = elapsed_seconds / n_partitions
        lines.append(f'  per partition:        {per_partition:,.1f} s')
        if parallelism:
            lines.append(f'  concurrent tasks:     ~{parallelism:,}')

        if n_total_partitions:
            same_width = elapsed_seconds * (n_total_partitions / n_partitions)
            lines += [
                '',
                f'Genome-wide estimate ({n_total_partitions:,} partitions, '
                f'{100 * n_partitions / n_total_partitions:.2f}% covered here):',
                f'  at this cluster width: {same_width / 3600:,.1f} h',
            ]
            if parallelism and full_parallelism and full_parallelism > parallelism:
                scaled = same_width * parallelism / full_parallelism
                lines += [
                    f'  at ~{full_parallelism:,} concurrent tasks: {scaled / 3600:,.1f} h',
                    '',
                    'The wider figure is the one to plan against: a probe covering a small',
                    'share of partitions often completes before autoscaling ramps, so its',
                    'observed concurrency understates a full run and its per-partition cost',
                    'is conservative.',
                ]
        lines += [
            '',
            'Do not scale by base pairs. A partition is the unit of work, one task streams',
            'one partition end to end, and a narrow probe touching few partitions is close',
            'to serial, so span-scaling overstates a full run by the width of the cluster.',
        ]
        if n_partitions < MIN_HEALTHY_PARTITIONS:
            lines += [
                '',
                f'CAUTION: only {n_partitions} partition(s) were read, so this timing is',
                'essentially serial and says little about a run that can use the whole',
                'cluster. Treat the per-partition number as the useful output here.',
            ]
    lines += [
        '',
        'What this covers: reading and aggregating. A scan writes only a small summary, so',
        'output cost is negligible and this is the whole picture.',
    ]
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Hail plumbing
# ---------------------------------------------------------------------------


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


def _open_read(path: str) -> IO[str]:
    """Open a local or GCS path for reading text."""
    if hl is not None and path.startswith('gs://'):
        return hl.hadoop_open(path, 'r')
    return open(path, 'rt')


def write_lines(path: str, header: str, rows: Iterable[str]) -> int:
    """Write a header plus rows, returning the number of rows written."""
    count = 0
    with _open_write(path) as handle:
        handle.write(header + '\n')
        for row in rows:
            handle.write(row + '\n')
            count += 1
    return count


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


def interval_span(intervals) -> int:
    """Total base span of a list of Hail locus intervals.

    ``build_intervals`` returns evaluated intervals, so the endpoints are plain ``Locus``
    objects with integer positions and need no further ``hl.eval``.
    """
    total = 0
    for interval in intervals:
        total += max(0, interval.end.position - interval.start.position)
    return total


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


def annotated_matrix(vds, mode: str, mapping, bin_size: int):
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
    matrix = matrix.annotate_rows(
        _bin=(matrix.locus.position - 1) // bin_size,
    )
    return matrix, n_unmatched, counts


def _metric_expression(matrix, mode: str):
    """The per-entry quantity being summed for the given mode."""
    _require_hail()
    if mode == 'variants':
        # LA is defined for every present entry. LGT is deliberately not used: a variant
        # call with call_GQ = 0 is imported with a missing LGT but a present entry, and
        # it corresponds to a real vet row, so it must be counted.
        return hl.agg.count_where(hl.is_defined(matrix.LA))
    return hl.agg.sum(
        hl.if_else(hl.is_defined(matrix.END), matrix.END - matrix.locus.position + 1, 0)
    )


def observed_parallelism() -> int | None:
    """Roughly how many tasks can run at once, or None if it cannot be determined.

    ``defaultParallelism`` tracks total executor cores, which is the right order of
    magnitude for turning a probe's wall clock into an estimate. Best-effort: it is only
    used to annotate a report, so any failure degrades the message rather than the run.
    """
    _require_hail()
    try:
        return int(hl.spark_context().defaultParallelism)
    except Exception:  # pragma: no cover - depends on a live Spark context
        return None


def executor_summary() -> str:
    """Live executor count and total task slots, for the log.

    Worth recording because it is the term most likely to be wrong in a cost estimate,
    and the hardest to reconstruct afterwards. A cluster that starts with two primary
    workers and depends on autoscaling for the rest can spend a long time far below its
    peak width, and nothing in the job output would otherwise say so -- an estimate that
    assumes full width is then wrong by the ratio, silently.
    """
    _require_hail()
    try:
        context = hl.spark_context()
        # Excludes the driver, which appears in this map but runs no tasks.
        executors = [
            executor for executor in
            context._jsc.sc().getExecutorMemoryStatus().keySet()
            if 'driver' not in str(executor).lower()
        ]
        return (f'{len(executors)} executor(s), '
                f'~{context.defaultParallelism} task slot(s)')
    except Exception as exc:  # pragma: no cover - needs a live Spark context
        return f'could not determine executor count ({exc})'


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
    62,000 x 134 numbers for a genome-wide 50 kb scan -- returns to the driver.  Each
    partition's accumulator stays small because a partition spans few bins.
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


@dataclass(frozen=True)
class ScanGeometry:
    """Partition counts and cluster width, the inputs a cost estimate needs."""

    n_partitions: int
    n_total_partitions: int
    parallelism: int | None
    full_parallelism: int | None


def _screening_matrix(args, vds):
    """Annotated matrix, per-superpartition sample counts, and the scan geometry."""
    source = args.sample_map_path
    with step(f'Reading sample mapping from {source}'):
        mapping = superpartition_column_table(source, args.superpartition_size)
    matrix, n_unmatched, counts = annotated_matrix(vds, args.mode, mapping, args.bin_size)

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

    geometry = ScanGeometry(
        n_partitions=n_partitions,
        n_total_partitions=n_total,
        parallelism=observed_parallelism(),
        # Peak width once the autoscaling policy has added its secondary workers.
        full_parallelism=(args.max_secondary_workers * args.worker_cores
                          if args.max_secondary_workers else None),
    )
    return matrix, counts, geometry


def action_scan(args) -> int:
    """Aggregate a VDS and write the summary and superpartition tables."""
    _require_hail()
    intervals = build_intervals(parse_contig_list(args.contigs),
                                parse_interval_list(args.intervals))
    announce(f'Scanning {len(intervals):,} interval(s) of {args.vds_path}')
    with step(f'Opening VDS {args.vds_path}'):
        vds = read_and_subset_vds(args.vds_path, intervals)
    matrix, counts, _ = _screening_matrix(args, vds)

    announce(f'Cluster width at start of aggregation: {executor_summary()}')
    with step(f'Aggregating {args.mode} into {args.bin_size:,} bp bins'):
        totals = aggregate_totals(matrix, args.mode, args.bin_size)
    # Logged again because autoscaling ramps: the two readings bracket the width the
    # run actually had, which is what a cost estimate needs and what one reading misses.
    announce(f'Cluster width at end of aggregation:   {executor_summary()}')
    rows = format_summary_rows(totals, args.bin_size)
    n_rows = write_lines(args.summary_path, SUMMARY_HEADER, rows)
    n_sp = write_lines(
        args.superpartitions_path, SUPERPARTITION_HEADER,
        [f'{sp}\t{counts[sp]}' for sp in sorted(counts)])

    print(f'Wrote {n_rows:,} summary rows ({len(totals):,} bins) to {args.summary_path}')
    print(f'Wrote {n_sp:,} superpartition rows to {args.superpartitions_path}')
    return 0


def action_probe(args) -> int:
    """Time a small scan and extrapolate to the whole genome."""
    _require_hail()
    interval_string = args.probe_interval
    intervals = build_intervals([], [interval_string])
    span = interval_span(intervals)
    announce(f'Probing {interval_string} ({span:,} bp) in {args.mode} mode')

    with step(f'Opening VDS {args.vds_path}'):
        vds = read_and_subset_vds(args.vds_path, intervals)
    matrix, counts, geometry = _screening_matrix(args, vds)
    n_samples = sum(counts.values())

    started = time.monotonic()
    with step(f'Aggregating {args.mode} over {interval_string}'):
        totals = aggregate_totals(matrix, args.mode, args.bin_size)
    elapsed = time.monotonic() - started

    n_cells = sum(len(v) for v in totals.values())
    print()
    print(format_probe_report(
        elapsed, span, n_cells, len(totals), n_samples,
        geometry.n_partitions, geometry.n_total_partitions,
        geometry.parallelism, geometry.full_parallelism))
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
    parser.add_argument('--max-secondary-workers', type=int, default=None,
                        help='Peak secondary workers the autoscaling policy may add. Used '
                             'only to annotate the probe report with the concurrency a full '
                             'run could reach.')
    parser.add_argument('--worker-cores', type=int, default=4,
                        help='Cores per worker, for the same estimate. Default: 4.')
    parser.add_argument('--probe-interval', default=DEFAULT_PROBE_INTERVAL,
                        help='Interval used by --action probe. Default: '
                             f'{DEFAULT_PROBE_INTERVAL}, a fully callable 10 Mb window on '
                             'chr20. Pick a gap-free interval; one overlapping a centromere '
                             'reads faster than average and skews the extrapolation.')

    parser.add_argument('--temp-path', default=None, help='Hail temporary directory.')
    return parser


REQUIRED_BY_ACTION: dict[str, tuple[str, ...]] = {
    'scan': ('summary_path', 'superpartitions_path'),
    'probe': (),
    'full-depth': ('full_depth_path',),
}


def validate_args(args, parser: argparse.ArgumentParser) -> None:
    """Fail fast on missing per-action arguments, before any cluster time is spent."""
    for name in REQUIRED_BY_ACTION[args.action]:
        if getattr(args, name, None) is None:
            parser.error(f"--{name.replace('_', '-')} is required for --action {args.action}")
    if args.bin_size < 1:
        parser.error('--bin-size must be at least 1')


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    _require_hail()
    hl.init(idempotent=True, tmp_dir=args.temp_path)
    hl.default_reference('GRCh38')

    handler = {
        'scan': action_scan,
        'probe': action_probe,
        'full-depth': action_full_depth,
    }[args.action]
    return handler(args)


if __name__ == '__main__':
    sys.exit(main())
