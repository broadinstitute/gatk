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
``--action materialize``
    Read the source VDS, keep a stratified sample of columns, and write a narrow copy.
    Costs one full-width read of the callset, which is unavoidable -- picking 100 of
    535K samples still requires looking at all of them -- and is the single large bill.
    Everything afterwards runs against the narrow copy at a small fraction of the cost.
    Also writes the chosen sample list and the superpartition table, so later runs and
    other VDSes can be screened against exactly the same samples.

``--action scan``
    Aggregate a VDS into per-(bin, superpartition) totals and write the summary and
    superpartition TSVs that `vds_dropout_detect.py` consumes.  Normally pointed at a
    narrow copy, but works on any VDS.

``--action probe``
    Run ``scan`` over a single small interval and report measured throughput plus a
    genome-wide extrapolation.  Cheap way to size the real runs before committing to
    them, and to check the Hail log for shuffle stages before paying for a big one.

    Which VDS this is pointed at decides which question it answers, and the more valuable
    use is the less obvious one:

    * Against the **source, full-width VDS** it measures full-width read throughput, which
      is what ``materialize`` costs.  That is the project's one large bill, so this is the
      run worth doing *first* -- before materializing, not after.
    * Against the **narrow copy** it estimates ongoing screening cost, which is small.

    Either way the estimate covers reading and aggregating, not writing, so it is a floor
    for ``materialize`` rather than a full prediction of it.

``--action full-depth``
    For named intervals and superpartitions, count per-sample data at full width and
    report how many samples carry nothing.  Meant for the handful of windows the screen
    flags, not for the genome.  Point this at the *source* VDS, not the narrow copy.

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

Sample map and sample list
--------------------------
Superpartition membership comes from a TSV of sample name to sample ID, since a VDS knows
only names.  ``GvsValidateVdsCompleteness.wdl`` generates that map itself when given a
BigQuery project and dataset; ``scripts/variantstore/bq/vds_dropout_sample_map.sql`` is
the hand-run copy, for building one outside the workflow.  Either tab- or comma-separated
input is accepted, so ``bq query --format=csv`` output works as-is.

``--action materialize`` writes the ``--sample-list-path`` file as a side effect, recording
which samples the stratified draw chose.  Later actions should be pointed at that file
rather than at the map, so every VDS is screened on identical samples -- which is what
makes figures comparable between r2 and r3.

Samples in the map but absent from the VDS are ignored, so a map generated after a
withdrawal is still usable.  A VDS sample *missing* from the map is a hard error instead,
because silently excluding it would bias the superpartition totals.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime
import hashlib
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
DEFAULT_SAMPLES_PER_SUPERPARTITION = 100
DEFAULT_STRIDE = 1
# Below this many partitions the aggregation cannot use the cluster, which on an AoU-scale
# VDS turns a nominally small interval into hours of single-threaded streaming.
MIN_HEALTHY_PARTITIONS = 8
DEFAULT_SEED = 'vs-1998'
# A 10 Mb window on chr20 that is 100% callable: comfortably clear of the centromere gap
# (chr20:26,386,232-30,088,349) and of the telomere, and free of the smaller assembly gaps
# near 31 Mb. Verified against
# wgs_calling_regions.hg38.noCentromeres.noTelomeres.interval_list. An interval overlapping
# a gap would make the probe's timing optimistic and emit confusing empty bins.
DEFAULT_PROBE_INTERVAL = 'chr20:37000000-47000000'
# End of the chr20 centromere gap, from the same interval list. Used to keep the default honest.
CHR20_CENTROMERE_END = 30_088_349

MODES = ('variants', 'references')
ACTIONS = ('materialize', 'scan', 'probe', 'full-depth')

# GRCh38 primary contigs GVS encodes. Alt and decoy contigs carry no GVS location
# encoding and are not screened.
GVS_CONTIGS = tuple([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])

# Approximate GRCh38 primary assembly length, used only to extrapolate probe timings.
GENOME_LENGTH = 3_088_269_832

SUMMARY_HEADER = 'contig\tbin_start\tbin_end\tsuperpartition\tobserved'
SUPERPARTITION_HEADER = 'superpartition\tn_samples'
SAMPLE_LIST_HEADER = 'sample_name\tsample_id\tsuperpartition'
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


def sample_sort_key(sample_name: str, seed: str = DEFAULT_SEED) -> str:
    """Stable pseudo-random ordering key for a sample.

    SHA-256 of ``seed:sample_name`` rather than ``hash()``, which is salted per process,
    or a Hail RNG, which would depend on partitioning.  The point is that re-running the
    selection -- next week, on another VDS, on another cluster -- yields the same samples,
    because the r3 comparison is only meaningful against the same sampled set.
    """
    return hashlib.sha256(f'{seed}:{sample_name}'.encode('utf-8')).hexdigest()


def stratified_sample(
        sample_to_superpartition: Mapping[str, int],
        samples_per_superpartition: int = DEFAULT_SAMPLES_PER_SUPERPARTITION,
        seed: str = DEFAULT_SEED,
) -> dict[int, list[str]]:
    """Choose up to ``samples_per_superpartition`` samples from each superpartition.

    Stratified rather than a global random draw: a global 2.5% draw averages the right
    count per superpartition but with binomial scatter, giving unequal detection power
    across superpartitions.  Taking a fixed count from each gives every superpartition
    the same power, and superpartitions smaller than the depth -- the last one in a
    callset -- simply contribute everything they have.

    Returns superpartition -> sorted list of chosen sample names.
    """
    if samples_per_superpartition < 1:
        raise ValueError(
            f"samples_per_superpartition must be positive, got {samples_per_superpartition}")

    grouped: dict[int, list[str]] = defaultdict(list)
    for sample_name, superpartition in sample_to_superpartition.items():
        grouped[superpartition].append(sample_name)

    chosen: dict[int, list[str]] = {}
    for superpartition, names in grouped.items():
        # Tie-break on the name itself rather than leaning on sort stability, which would
        # make the result depend on the order samples happen to arrive in. The selection has
        # to be reproducible across runs and across VDSes for the r2/r3 comparison to mean
        # anything, so it must not depend on dict iteration order.
        names.sort(key=lambda name: (sample_sort_key(name, seed), name))
        chosen[superpartition] = sorted(names[:samples_per_superpartition])
    return chosen


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


def assign_superpartitions(
        sample_map: Mapping[str, int],
        present_samples: Iterable[str],
        superpartition_size: int = DEFAULT_SUPERPARTITION_SIZE,
) -> dict[str, int]:
    """Restrict the sample map to samples present in the VDS and assign superpartitions.

    Samples in the VDS but missing from the map are a hard error: they would be silently
    excluded from every superpartition total, quietly biasing the peer comparison the
    detector relies on.
    """
    present = list(present_samples)
    missing = [s for s in present if s not in sample_map]
    if missing:
        raise ValueError(
            f"{len(missing)} sample(s) in the VDS are absent from the sample map, "
            f"starting with {missing[:5]}; regenerate the map from sample_info"
        )
    return {s: superpartition_for(sample_map[s], superpartition_size) for s in present}


def stride_keeps_bin(bin_index: int, stride: int) -> bool:
    """Whether a strided scan reads the bin at ``bin_index``.

    A stride of 1 reads everything.  Larger strides read one bin in every ``stride``,
    which guarantees that any dropout wider than ``stride * bin_size`` fully contains at
    least one scanned bin, while saying nothing about the unscanned remainder.
    """
    if stride < 1:
        raise ValueError(f"stride must be at least 1, got {stride}")
    return bin_index % stride == 0


def strided_intervals(contig: str, start: int, end: int, bin_size: int,
                      stride: int) -> list[str]:
    """Express a strided scan as intervals, one per retained bin.

    Striding has to reach ``read_vds`` as intervals to save anything. Expressed instead as
    a row predicate on the bin index, Hail cannot consult the row index and streams every
    byte of every partition it opens, discarding most of it after the fact -- so the scan
    costs the same as reading everything and merely reports less. Intervals let the reader
    seek.
    """
    if stride < 1:
        raise ValueError(f"stride must be at least 1, got {stride}")
    intervals = []
    first = bin_index_for(start, bin_size)
    last = bin_index_for(end - 1, bin_size)
    for index in range(first, last + 1):
        if not stride_keeps_bin(index, stride):
            continue
        bin_start = max(start, index * bin_size + 1)
        bin_end = min(end, (index + 1) * bin_size + 1)
        if bin_end > bin_start:
            intervals.append(f'{contig}:{bin_start}-{bin_end}')
    return intervals


def parse_bool(value: str | bool) -> bool:
    """Parse an explicit boolean.

    Boolean flags are spelled with a value rather than as ``store_true`` because
    `run_in_hail_cluster.py` renders every key of its arguments JSON as ``--key value``;
    a bare ``store_true`` flag would receive an argument and fail to parse.
    """
    if isinstance(value, bool):
        return value
    lowered = value.strip().lower()
    if lowered in ('true', 't', 'yes', '1'):
        return True
    if lowered in ('false', 'f', 'no', '0'):
        return False
    raise ValueError(f"expected a boolean, got {value!r}")


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
    strided or sparse region, and `vds_dropout_detect.py` restores them from the
    superpartition table.  Rows are ordered by contig then position so the output is
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


def format_sample_list_rows(
        chosen: Mapping[int, Sequence[str]],
        sample_map: Mapping[str, int],
) -> list[str]:
    """Render the chosen sample list as TSV lines, for reuse across VDSes."""
    rows: list[str] = []
    for superpartition in sorted(chosen):
        for name in chosen[superpartition]:
            rows.append(f'{name}\t{sample_map[name]}\t{superpartition}')
    return rows


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
        'What this covers: reading and aggregating, not writing. Against a full-width VDS',
        'it is therefore a floor for the materialize pass, not a full prediction of it.',
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


def apply_stride(interval_strings: Sequence[str], bin_size: int, stride: int) -> list[str]:
    """Expand intervals into the strided subset actually to be read.

    Returns them unchanged when stride is 1. See ``strided_intervals`` for why this has to
    happen here, at read time, rather than as a row filter later.
    """
    if stride <= 1:
        return list(interval_strings)
    expanded: list[str] = []
    for text in interval_strings:
        contig, _, span = text.partition(':')
        if not span:
            raise ValueError(
                f'striding needs bounded intervals; {text!r} covers a whole contig. '
                'Pass explicit intervals, or use stride 1.')
        start_text, _, end_text = span.partition('-')
        expanded.extend(strided_intervals(
            contig, int(start_text.replace(',', '')), int(end_text.replace(',', '')),
            bin_size, stride))
    return expanded


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


def vds_sample_names(vds) -> list[str]:
    """Collect the VDS column keys.

    Reads only the column metadata, not entries, so this is cheap even on a full-width
    AoU VDS.
    """
    _require_hail()
    return vds.variant_data.cols().s.collect()


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


def action_materialize(args) -> int:
    """Write a stratified narrow copy of the VDS, plus the sample and superpartition tables."""
    _require_hail()
    with step(f'Opening VDS {args.vds_path}'):
        vds = hl.vds.read_vds(args.vds_path)
    with step(f'Reading sample map {args.sample_map_path}'):
        sample_map = read_sample_map(args.sample_map_path)
    # The one remaining driver-side collect: stratified selection needs the full name list.
    with step('Collecting VDS sample names'):
        present = vds_sample_names(vds)
    print(f'VDS holds {len(present):,} samples; map holds {len(sample_map):,}')

    assigned = assign_superpartitions(sample_map, present, args.superpartition_size)
    chosen = stratified_sample(assigned, args.samples_per_superpartition, args.random_seed)
    selected = sorted(name for names in chosen.values() for name in names)
    print(
        f'Selected {len(selected):,} samples across {len(chosen):,} superpartitions '
        f'(depth {args.samples_per_superpartition}, seed {args.random_seed!r})'
    )

    keep = hl.Table.parallelize(
        [{'s': name} for name in selected], schema=hl.tstruct(s=hl.tstr),
    ).key_by('s')
    narrow = hl.vds.filter_samples(vds, keep, keep=True)

    # Rows carrying no data for any retained sample are dead weight. Rare variants
    # dominate the callset, so most loci have no sampled carrier and drop out here,
    # shrinking the copy well beyond the column ratio alone. Safe for this purpose:
    # the detector counts entries per bin, and an all-missing row contributes zero
    # either way.
    variant_data = narrow.variant_data.filter_rows(
        hl.agg.any(hl.is_defined(narrow.variant_data.LA)))
    reference_data = narrow.reference_data.filter_rows(
        hl.agg.any(hl.is_defined(narrow.reference_data.END)))
    with step(f'Writing narrow VDS to {args.output_path}'):
        hl.vds.VariantDataset(reference_data, variant_data).write(
            args.output_path, overwrite=args.overwrite)

    n_sp = write_lines(
        args.superpartitions_path, SUPERPARTITION_HEADER, format_superpartition_rows(chosen))
    n_samples = write_lines(
        args.sample_list_path, SAMPLE_LIST_HEADER, format_sample_list_rows(chosen, sample_map))
    print(f'Wrote {n_sp:,} superpartition rows to {args.superpartitions_path}')
    print(f'Wrote {n_samples:,} sample rows to {args.sample_list_path}')
    return 0


def sample_source_path(args) -> str:
    """File the superpartition mapping is read from, preferring the recorded sample list.

    The recorded list is preferred so that every VDS is screened on identical samples,
    which is what makes figures comparable between r2 and r3.
    """
    return args.sample_list_path or args.sample_map_path


@dataclass(frozen=True)
class ScanGeometry:
    """Partition counts and cluster width, the inputs a cost estimate needs."""

    n_partitions: int
    n_total_partitions: int
    parallelism: int | None
    full_parallelism: int | None


def _screening_matrix(args, vds):
    """Annotated matrix, per-superpartition sample counts, and the scan geometry."""
    source = sample_source_path(args)
    with step(f'Reading sample mapping from {source}'):
        mapping = superpartition_column_table(source, args.superpartition_size)
    matrix, n_unmatched, counts = annotated_matrix(vds, args.mode, mapping, args.bin_size)

    if n_unmatched:
        detail = f'{n_unmatched:,} VDS sample(s) are absent from {source}'
        if args.sample_list_path:
            # A recorded list is deliberately a subset of the callset.
            print(f'Note: {detail}; they will not be screened.')
        else:
            raise ValueError(
                f'{detail}. A sample map is expected to cover every sample in the VDS, and '
                'screening only part of a superpartition would bias the peer comparison the '
                'detector relies on. Regenerate the map from sample_info, or pass a recorded '
                'sample list if screening a subset is intended.')

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
    requested = parse_interval_list(args.intervals)
    strided = apply_stride(requested, args.bin_size, args.stride) if requested else requested
    if args.stride > 1 and not requested:
        raise ValueError('--stride requires --intervals; a whole contig cannot be strided '
                         'without bounds. Pass explicit intervals or use --stride 1.')
    intervals = build_intervals(parse_contig_list(args.contigs), strided)
    announce(f'Scanning {len(intervals):,} interval(s) of {args.vds_path}')
    with step(f'Opening VDS {args.vds_path}'):
        vds = read_and_subset_vds(args.vds_path, intervals)
    matrix, counts, _ = _screening_matrix(args, vds)

    with step(f'Aggregating {args.mode} into {args.bin_size:,} bp bins'):
        totals = aggregate_totals(matrix, args.mode, args.bin_size)
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

    with step(f'Reading sample mapping from {sample_source_path(args)}'):
        mapping = superpartition_column_table(sample_source_path(args), args.superpartition_size)

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
                        help='VDS to read. The source VDS for materialize and full-depth; '
                             'normally the narrow copy for scan and probe.')
    parser.add_argument('--mode', choices=MODES, default='variants',
                        help='Metric to compute. Default: variants.')

    parser.add_argument('--sample-map-path', default=None,
                        help='TSV of sample_name and sample_id from sample_info.')
    parser.add_argument('--sample-list-path', default=None,
                        help='Chosen-sample TSV. Written by materialize, and preferred by '
                             'later actions so every VDS is screened on the same samples.')
    parser.add_argument('--superpartitions-path', default=None,
                        help='Where to write the superpartition/n_samples table.')
    parser.add_argument('--summary-path', default=None,
                        help='Where to write the summary table (scan).')
    parser.add_argument('--full-depth-path', default=None,
                        help='Where to write per-sample results (full-depth).')
    parser.add_argument('--output-path', default=None,
                        help='Destination VDS path (materialize).')

    parser.add_argument('--samples-per-superpartition', type=int,
                        default=DEFAULT_SAMPLES_PER_SUPERPARTITION,
                        help=f'Stratified sampling depth. Default: '
                             f'{DEFAULT_SAMPLES_PER_SUPERPARTITION}.')
    parser.add_argument('--random-seed', default=DEFAULT_SEED,
                        help=f'Seed for sample selection. Default: {DEFAULT_SEED!r}. Keep it '
                             'fixed across VDSes so results stay comparable.')
    parser.add_argument('--superpartition-size', type=int, default=DEFAULT_SUPERPARTITION_SIZE,
                        help=f'Samples per GVS superpartition. Default: '
                             f'{DEFAULT_SUPERPARTITION_SIZE}.')
    parser.add_argument('--bin-size', type=int, default=DEFAULT_BIN_SIZE,
                        help=f'Genomic bin size in bp. Default: {DEFAULT_BIN_SIZE}.')
    parser.add_argument('--stride', type=int, default=DEFAULT_STRIDE,
                        help='Read one bin in every N. Default 1 (read everything). Larger '
                             'values cost about 1/N but say nothing about unscanned bins.')

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
    parser.add_argument('--overwrite', type=parse_bool, default=False,
                        help='Overwrite the destination VDS (materialize). Takes an explicit '
                             'true/false because the cluster runner renders every argument '
                             'as --key value.')
    return parser


REQUIRED_BY_ACTION: dict[str, tuple[str, ...]] = {
    'materialize': ('output_path', 'superpartitions_path', 'sample_list_path'),
    'scan': ('summary_path', 'superpartitions_path'),
    'probe': (),
    'full-depth': ('full_depth_path',),
}


def validate_args(args, parser: argparse.ArgumentParser) -> None:
    """Fail fast on missing per-action arguments, before any cluster time is spent."""
    for name in REQUIRED_BY_ACTION[args.action]:
        if getattr(args, name, None) is None:
            parser.error(f"--{name.replace('_', '-')} is required for --action {args.action}")
    if args.action == 'materialize' and not args.sample_map_path:
        parser.error('--sample-map-path is required for --action materialize')
    if args.action != 'materialize' and not (args.sample_map_path or args.sample_list_path):
        parser.error(
            f'--action {args.action} needs either --sample-list-path (preferred) or '
            '--sample-map-path'
        )
    if args.stride < 1:
        parser.error('--stride must be at least 1')
    if args.bin_size < 1:
        parser.error('--bin-size must be at least 1')
    if args.samples_per_superpartition < 1:
        parser.error('--samples-per-superpartition must be at least 1')


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    _require_hail()
    hl.init(idempotent=True, tmp_dir=args.temp_path)
    hl.default_reference('GRCh38')

    handler = {
        'materialize': action_materialize,
        'scan': action_scan,
        'probe': action_probe,
        'full-depth': action_full_depth,
    }[args.action]
    return handler(args)


if __name__ == '__main__':
    sys.exit(main())
