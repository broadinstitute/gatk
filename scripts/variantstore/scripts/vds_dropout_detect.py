#!/usr/bin/env python3
"""
Detect "rectangle" data dropouts in a GVS Hail VDS from a pre-computed summary table.

A rectangle dropout is a contiguous genomic window in which one GVS superpartition has
little or no data while every other superpartition has the usual amount.  This is the
shape produced when a single Avro export shard is lost, truncated, or never read: the
`EXPORT DATA ... ORDER BY location` in `GvsExtractAvroFilesForHail.wdl` means each
numbered output file holds a contiguous location range for exactly one superpartition.

This module does no Hail work and reads no VDS.  It consumes the small summary table
emitted by `vds_dropout_scan.py` so that thresholds can be re-tuned offline, for free,
without re-reading a multi-terabyte VDS -- and so that the judging logic can be unit
tested without Hail or a real VDS.

Input files
-----------
Summary TSV (``--summary``), one row per non-empty (bin, superpartition) cell::

    contig  bin_start  bin_end   superpartition  observed
    chr4    56550001   56600001  83              138

``bin_start`` is inclusive and ``bin_end`` exclusive, both 1-based to match Hail locus
positions.  ``observed`` is the metric being screened: defined variant-data entries in
``variants`` mode, covered reference bases in ``references`` mode.  Cells absent from the
file are taken to be zero, which is why the superpartition universe has to be declared
separately.

Superpartition TSV (``--superpartitions``)::

    superpartition  n_samples
    83              100

``n_samples`` is how many samples that superpartition contributes. Usually 4000, but the
final superpartition of a callset holds fewer, and withdrawals thin the others -- so cells
are converted to per-sample rates before anything is compared.

Method
------
1. Convert each cell to a per-sample rate, ``observed / n_samples``, so superpartitions
   contributing unequal sample counts are comparable.
2. Divide by a high-quantile rate across superpartitions within the bin.  This removes
   regional variation in variant density, and neutralizes regions that are hard for
   every sample -- a centromere depresses every superpartition equally, so the ratio
   stays near 1.

   The baseline is the 75th percentile rather than the median deliberately.  A median
   baseline collapses to zero once more than half the superpartitions in a bin are
   empty, and such a bin is then discarded as uninformative -- which would silently
   swallow a dropout that happened to hit many superpartitions at once.  The 75th
   percentile stays usable until three quarters of them are empty, and costs only about
   3% inflation of the expected value in a clean region, applied uniformly.
3. Divide by that superpartition's median ratio across all bins.  This removes
   superpartition-level scale differences.  Superpartitions are assigned by ingest order
   and are not ancestry-balanced, so a superpartition can legitimately carry 30% more
   variants than its peers everywhere; that is a global offset, not a regional dropout.
4. Flag cells whose residual is low, subject to an evidence floor (see below).
5. Merge flagged cells that are adjacent in the same superpartition into rectangles.

Step 3 has a blind spot worth being explicit about: a superpartition depleted across
*most* of the genome would have its own median dragged down and its residuals would come
back near 1.  That case is caught separately by comparing each superpartition's global
scale against its peers, reported as a distinct finding type.

Severity is graded rather than boolean.  A threshold tuned to "residual is essentially
zero" would miss the Foxtrot r2 state, where the affected windows were thinned rather
than emptied, and would equally miss an incomplete repair.  Every finding therefore
carries a ratio and a depletion score, and the report is ranked by severity.

The depletion score is ``(expected - observed) / sqrt(expected)``, the Poisson-flavoured
size of the shortfall.  It is a heuristic ranking and evidence floor, not a calibrated
p-value -- variant counts are overdispersed across samples, so the nominal
interpretation would be optimistic.  It exists to stop thin, low-count bins from firing
on noise; correctness of an individual finding is settled by BigQuery adjudication, not
by this score.

Usage
-----
::

    python3 vds_dropout_detect.py \\
        --summary         gs://.../summary.tsv \\
        --superpartitions gs://.../superpartitions.tsv \\
        --mode            variants \\
        --project-id      aou-genomics-curation-prod \\
        --dataset-name    foxtrot \\
        --report-path     findings.tsv \\
        --sql-path        adjudicate.sql
"""

from __future__ import annotations

import argparse
import gzip
import math
import statistics
import sys
from array import array
from dataclasses import dataclass, field
from typing import IO, Iterable, Sequence

# GVS encodes a genomic location as chromosome_index * 1e12 + position.  See
# SchemaUtils.java; chrX is 23 and chrY is 24.
LOCATION_MULTIPLIER = 1_000_000_000_000

_CONTIG_INDEX: dict[str, int] = {f'chr{i}': i for i in range(1, 23)}
_CONTIG_INDEX['chrX'] = 23
_CONTIG_INDEX['chrY'] = 24

# Defaults.  All are cheap to revisit because this stage runs offline against a small
# summary table rather than against the VDS.
DEFAULT_BASELINE_QUANTILE = 0.75
DEFAULT_RATIO_THRESHOLD = 0.5
DEFAULT_SCORE_THRESHOLD = 8.0
DEFAULT_MIN_EXPECTED = 30.0
# Superpartitions are assigned by ingest order and are not ancestry- or batch-balanced,
# so real genome-wide variation between them can reach 20-30%. This is therefore the
# threshold most likely to need tuning against real data; it is deliberately loose.
DEFAULT_SUPERPARTITION_SCALE_THRESHOLD = 0.7

# Adjudication queries are generated for the worst candidates only. A genome-wide scan
# examines millions of cells, and if thresholds turn out loose the report could name
# thousands of rectangles -- more queries than anyone will run, burying the real findings.
# The cap is on generated SQL only; every rectangle still appears in the report, and the
# truncation is stated explicitly rather than left to be noticed.
DEFAULT_MAX_SQL_QUERIES = 50

SUMMARY_COLUMNS = ('contig', 'bin_start', 'bin_end', 'superpartition', 'observed')
SUPERPARTITION_COLUMNS = ('superpartition', 'n_samples')

MODES = ('variants', 'references')

# AoU callsets use the compressed reference schema, so that is the default. The compressed
# ref_ranges tables hold only (packed_ref_data, sample_id) -- there is no `location` column --
# and packed_ref_data is the clustering field (GvsCreateTables.wdl:32).
REFERENCE_SCHEMAS = ('compressed', 'uncompressed')
DEFAULT_REFERENCE_SCHEMA = 'compressed'

# Bit layout of packed_ref_data. The canonical producer is
# SchemaUtils#encodeCompressedRefBlock (SchemaUtils.java:114); UnpackRefRangeInfo in
# GvsExtractAvroFilesForHail.wdl is the SQL-side consumer of the same layout.
#   bits 48-63  chromosome index  (16 bits)
#   bits 16-47  position          (32 bits)
#   bits  4-15  block length      (12 bits)
#   bits  0-3   state             (4 bits)
#
# Because chromosome and position occupy the top 48 bits, numeric ordering of the packed
# value is identical to ordering by (chromosome, position) -- length and state only break
# ties within a single position. That is why clustering on packed_ref_data prunes a
# location range just as well as a clustered `location` column would, and why a range
# predicate has to be written against the packed value rather than a decoded expression.
PACKED_CHROMOSOME_SHIFT = 48
PACKED_POSITION_SHIFT = 16
PACKED_CHROMOSOME_MASK = 0xFFFF
PACKED_POSITION_MASK = 0xFFFFFFFF
# Length and state together occupy the low 16 bits.
PACKED_LOW_BITS_MASK = 0xFFFF


class SummaryFormatError(ValueError):
    """Raised when an input table is malformed, with the offending line identified."""


def contig_index(contig: str) -> int:
    """Return the GVS chromosome index for a contig name.

    Raises ``KeyError`` for contigs GVS does not encode, which includes the alt and
    decoy contigs; callers screening a whole-genome VDS should expect to skip those.
    """
    return _CONTIG_INDEX[contig]


def encode_location(contig: str, position: int) -> int:
    """Encode a contig and 1-based position into a GVS location."""
    return contig_index(contig) * LOCATION_MULTIPLIER + position


def decode_packed_location(packed: int) -> int:
    """Decode the GVS location out of a ``packed_ref_data`` value.

    Deliberately a transcription of ``UnpackRefRangeInfo`` in
    `GvsExtractAvroFilesForHail.wdl`, so the two can be compared by eye.
    """
    chromosome = (packed >> PACKED_CHROMOSOME_SHIFT) & PACKED_CHROMOSOME_MASK
    position = (packed >> PACKED_POSITION_SHIFT) & PACKED_POSITION_MASK
    return LOCATION_MULTIPLIER * chromosome + position


def packed_ref_data_bounds(contig: str, start_position: int, end_position: int) -> tuple[int, int]:
    """Inclusive ``packed_ref_data`` bounds for blocks *starting* in a position range.

    Because the chromosome and position occupy the high bits, a contiguous position range
    maps to a contiguous packed range: the low 16 bits hold only length and state, so the
    bounds are the range endpoints with those bits cleared and set respectively.

    Filtering on the packed value rather than on a decoded location expression is what
    keeps the query cheap -- ``packed_ref_data`` is the clustering field, and a predicate
    on a computed expression would prune nothing and scan the whole table.

    These are *start* positions, matching how the scan attributes each reference block
    entirely to the bin holding its start locus. A block beginning just before the window
    and extending into it is therefore excluded on both sides, so the two figures stay
    comparable; with blocks capped at 1000 bp and windows around 500 kb the effect is at
    most one block per edge.
    """
    if end_position < start_position:
        raise ValueError(
            f"end_position {end_position} precedes start_position {start_position}")
    chromosome = contig_index(contig) << PACKED_CHROMOSOME_SHIFT
    low = chromosome | (start_position << PACKED_POSITION_SHIFT)
    high = chromosome | (end_position << PACKED_POSITION_SHIFT) | PACKED_LOW_BITS_MASK
    return low, high


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------


@dataclass
class Summary:
    """A parsed summary table, held as a dense matrix of per-bin rates.

    ``observed[i][j]`` is the metric for bin ``bins[i]`` and superpartition
    ``superpartitions[j]``.  Stored as ``array('d')`` rows rather than nested dicts
    because a genome-wide 50 kb scan is on the order of 8 million cells, which is
    tolerable as packed doubles and wasteful as Python objects.
    """

    superpartitions: list[int]
    n_samples: dict[int, int]
    bins: list[tuple[str, int, int]] = field(default_factory=list)
    observed: list[array] = field(default_factory=list)

    def superpartition_position(self, superpartition: int) -> int:
        return self.superpartitions.index(superpartition)


def _open_text(path: str) -> IO[str]:
    """Open a possibly-gzipped local path for reading text."""
    if path.endswith('.gz') or path.endswith('.bgz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


def _split_header(line: str, expected: Sequence[str], path: str) -> list[str]:
    header = line.rstrip('\n').split('\t')
    if tuple(header[:len(expected)]) != tuple(expected):
        raise SummaryFormatError(
            f"{path}: expected header {list(expected)}, found {header[:len(expected)]}"
        )
    return header


def load_superpartitions(path: str) -> dict[int, int]:
    """Load the superpartition -> sampled-sample-count map."""
    result: dict[int, int] = {}
    with _open_text(path) as handle:
        first = handle.readline()
        if not first:
            raise SummaryFormatError(f"{path}: file is empty")
        _split_header(first, SUPERPARTITION_COLUMNS, path)
        for lineno, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 2:
                raise SummaryFormatError(f"{path}:{lineno}: expected 2 columns, found {len(fields)}")
            try:
                superpartition = int(fields[0])
                n_samples = int(fields[1])
            except ValueError as exc:
                raise SummaryFormatError(f"{path}:{lineno}: {exc}") from exc
            if superpartition in result:
                raise SummaryFormatError(f"{path}:{lineno}: duplicate superpartition {superpartition}")
            result[superpartition] = n_samples
    if not result:
        raise SummaryFormatError(f"{path}: no superpartitions declared")
    return result


def load_summary(summary_path: str, n_samples: dict[int, int]) -> Summary:
    """Stream the summary TSV into a dense matrix, filling absent cells with zero.

    Rows for superpartitions not present in ``n_samples`` are a hard error rather than a
    silent skip: it means the summary and the superpartition declaration came from
    different runs, and every downstream median would be computed over the wrong
    universe.
    """
    superpartitions = sorted(n_samples)
    position = {sp: i for i, sp in enumerate(superpartitions)}
    width = len(superpartitions)

    summary = Summary(superpartitions=superpartitions, n_samples=dict(n_samples))
    bin_position: dict[tuple[str, int, int], int] = {}

    with _open_text(summary_path) as handle:
        first = handle.readline()
        if not first:
            raise SummaryFormatError(f"{summary_path}: file is empty")
        _split_header(first, SUMMARY_COLUMNS, summary_path)
        for lineno, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 5:
                raise SummaryFormatError(
                    f"{summary_path}:{lineno}: expected 5 columns, found {len(fields)}"
                )
            contig = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
                superpartition = int(fields[3])
                observed = float(fields[4])
            except ValueError as exc:
                raise SummaryFormatError(f"{summary_path}:{lineno}: {exc}") from exc

            if superpartition not in position:
                raise SummaryFormatError(
                    f"{summary_path}:{lineno}: superpartition {superpartition} is not declared in "
                    "the superpartition table; the two inputs are from different runs"
                )
            if end <= start:
                raise SummaryFormatError(
                    f"{summary_path}:{lineno}: bin_end {end} must exceed bin_start {start}"
                )

            key = (contig, start, end)
            index = bin_position.get(key)
            if index is None:
                index = len(summary.bins)
                bin_position[key] = index
                summary.bins.append(key)
                summary.observed.append(array('d', [0.0]) * width if width else array('d'))
            summary.observed[index][position[superpartition]] = observed

    return summary


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Cell:
    """One flagged (bin, superpartition) cell."""

    contig: str
    start: int
    end: int
    superpartition: int
    n_samples: int
    observed: float
    expected: float

    @property
    def ratio(self) -> float:
        return self.observed / self.expected if self.expected > 0 else 1.0

    @property
    def score(self) -> float:
        return depletion_score(self.observed, self.expected)


@dataclass(frozen=True)
class Rectangle:
    """A run of adjacent flagged cells in one superpartition on one contig."""

    contig: str
    start: int
    end: int
    superpartition: int
    n_samples: int
    n_bins: int
    observed: float
    expected: float

    @property
    def ratio(self) -> float:
        return self.observed / self.expected if self.expected > 0 else 1.0

    @property
    def score(self) -> float:
        return depletion_score(self.observed, self.expected)

    @property
    def span(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class SuperpartitionScale:
    """A superpartition whose global data volume is out of line with its peers.

    This is the finding type that covers the blind spot in per-bin normalization: a
    superpartition depleted across most of the genome would look unremarkable bin by bin.
    """

    superpartition: int
    n_samples: int
    # Data volume relative to a typical superpartition, where 1.0 is typical. Normalized
    # against the median superpartition rather than the raw per-bin baseline, so the
    # value means the same thing regardless of --baseline-quantile.
    relative_scale: float


@dataclass
class Report:
    rectangles: list[Rectangle] = field(default_factory=list)
    superpartition_scales: list[SuperpartitionScale] = field(default_factory=list)
    n_bins_considered: int = 0
    n_bins_skipped_empty: int = 0
    n_cells_flagged: int = 0

    @property
    def clean(self) -> bool:
        return not self.rectangles and not self.superpartition_scales


def depletion_score(observed: float, expected: float) -> float:
    """Poisson-flavoured size of a shortfall: ``(expected - observed) / sqrt(expected)``.

    Zero when there is no shortfall.  A heuristic ranking and evidence floor rather than
    a calibrated statistic; see the module docstring.
    """
    if expected <= 0:
        return 0.0
    shortfall = expected - observed
    if shortfall <= 0:
        return 0.0
    return shortfall / math.sqrt(expected)


def _median(values: Iterable[float]) -> float:
    collected = list(values)
    if not collected:
        return 0.0
    return statistics.median(collected)


def _quantile(values: Iterable[float], q: float) -> float:
    """Linearly interpolated quantile over ``values``, or 0.0 when empty.

    Written out rather than delegated to ``statistics.quantiles`` so the behavior at
    small sample sizes and at the extremes is explicit and does not vary by Python
    version.
    """
    ordered = sorted(values)
    if not ordered:
        return 0.0
    if q <= 0.0:
        return ordered[0]
    if q >= 1.0:
        return ordered[-1]
    position = q * (len(ordered) - 1)
    lower = math.floor(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def bin_baseline_rates(
        summary: Summary,
        quantile: float = DEFAULT_BASELINE_QUANTILE,
) -> list[float]:
    """Baseline per-sample rate across superpartitions, for each bin.

    Superpartitions contributing no samples are excluded.  Zero-valued cells are
    *included*: excluding them would hide exactly the signal being screened for.

    See the module docstring for why this is a high quantile rather than the median --
    briefly, a median baseline goes to zero as soon as half the superpartitions in a bin
    are empty, and the bin is then dropped as uninformative.
    """
    counts = [summary.n_samples[sp] for sp in summary.superpartitions]
    baselines: list[float] = []
    for row in summary.observed:
        rates = [row[j] / counts[j] for j in range(len(counts)) if counts[j] > 0]
        baselines.append(_quantile(rates, quantile))
    return baselines


def superpartition_scales(summary: Summary, bin_baselines: Sequence[float]) -> list[float]:
    """Median normalized ratio per superpartition, across all informative bins.

    A value near 1 means the superpartition carries a typical amount of data overall.
    Bins with no usable baseline carry no information and are excluded.

    This stays a median across bins: a superpartition dropped in a handful of windows is
    unaffected, while one dropped across most of the genome has its median pulled down,
    which is what makes the scale check catch the case per-bin residuals cannot.
    """
    counts = [summary.n_samples[sp] for sp in summary.superpartitions]
    per_sp: list[list[float]] = [[] for _ in summary.superpartitions]
    for i, row in enumerate(summary.observed):
        baseline_rate = bin_baselines[i]
        if baseline_rate <= 0:
            continue
        for j, count in enumerate(counts):
            if count <= 0:
                continue
            per_sp[j].append((row[j] / count) / baseline_rate)
    return [_median(values) if values else 0.0 for values in per_sp]


def flag_cells(
        summary: Summary,
        bin_baselines: Sequence[float],
        scales: Sequence[float],
        ratio_threshold: float = DEFAULT_RATIO_THRESHOLD,
        score_threshold: float = DEFAULT_SCORE_THRESHOLD,
        min_expected: float = DEFAULT_MIN_EXPECTED,
) -> tuple[list[Cell], int, int]:
    """Flag depleted cells.

    Returns the flagged cells plus counts of bins considered and bins skipped for
    carrying no data at all.

    The ``min_expected`` floor is what stops thin bins from firing on noise: a cell
    expecting five entries and observing none is unremarkable, while one expecting ten
    thousand and observing five hundred is not.
    """
    counts = [summary.n_samples[sp] for sp in summary.superpartitions]
    flagged: list[Cell] = []
    considered = 0
    skipped = 0

    for i, row in enumerate(summary.observed):
        baseline_rate = bin_baselines[i]
        if baseline_rate <= 0:
            skipped += 1
            continue
        considered += 1
        contig, start, end = summary.bins[i]
        for j, superpartition in enumerate(summary.superpartitions):
            count = counts[j]
            if count <= 0:
                continue
            scale = scales[j]
            if scale <= 0:
                # Handled as a superpartition-level finding; per-bin residuals are
                # meaningless when the superpartition has no usable global scale.
                continue
            expected = baseline_rate * scale * count
            if expected < min_expected:
                continue
            observed = row[j]
            if expected > 0 and observed / expected >= ratio_threshold:
                continue
            if depletion_score(observed, expected) < score_threshold:
                continue
            flagged.append(Cell(
                contig=contig,
                start=start,
                end=end,
                superpartition=superpartition,
                n_samples=count,
                observed=observed,
                expected=expected,
            ))

    return flagged, considered, skipped


def merge_cells(cells: Sequence[Cell]) -> list[Rectangle]:
    """Merge cells adjacent in the same superpartition and contig into rectangles.

    Two cells merge when they share a contig and superpartition and the later one starts
    no further along than the earlier one ends, so touching bins join and a gap breaks
    the run.
    """
    ordered = sorted(cells, key=lambda c: (c.superpartition, c.contig, c.start, c.end))
    rectangles: list[Rectangle] = []

    run: list[Cell] = []

    def flush() -> None:
        if not run:
            return
        rectangles.append(Rectangle(
            contig=run[0].contig,
            start=run[0].start,
            end=run[-1].end,
            superpartition=run[0].superpartition,
            n_samples=run[0].n_samples,
            n_bins=len(run),
            observed=sum(c.observed for c in run),
            expected=sum(c.expected for c in run),
        ))
        run.clear()

    for cell in ordered:
        if run and (
                cell.superpartition == run[-1].superpartition
                and cell.contig == run[-1].contig
                and cell.start <= run[-1].end
        ):
            run.append(cell)
        else:
            flush()
            run.append(cell)
    flush()

    rectangles.sort(key=lambda r: r.score, reverse=True)
    return rectangles


def flag_superpartition_scales(
        summary: Summary,
        scales: Sequence[float],
        scale_threshold: float = DEFAULT_SUPERPARTITION_SCALE_THRESHOLD,
) -> list[SuperpartitionScale]:
    """Flag superpartitions whose overall data volume is far below their peers.

    Scales are re-normalized against the median superpartition before comparison, so the
    threshold reads as "carries less than this fraction of what a typical superpartition
    carries" and does not shift when --baseline-quantile changes.
    """
    contributing = [
        j for j, sp in enumerate(summary.superpartitions) if summary.n_samples[sp] > 0
    ]
    typical = _median([scales[j] for j in contributing if scales[j] > 0])
    if typical <= 0:
        return []

    findings: list[SuperpartitionScale] = []
    for j in contributing:
        superpartition = summary.superpartitions[j]
        relative = scales[j] / typical
        if relative < scale_threshold:
            findings.append(SuperpartitionScale(
                superpartition=superpartition,
                n_samples=summary.n_samples[superpartition],
                relative_scale=relative,
            ))
    findings.sort(key=lambda f: f.relative_scale)
    return findings


def analyze(
        summary: Summary,
        ratio_threshold: float = DEFAULT_RATIO_THRESHOLD,
        score_threshold: float = DEFAULT_SCORE_THRESHOLD,
        min_expected: float = DEFAULT_MIN_EXPECTED,
        scale_threshold: float = DEFAULT_SUPERPARTITION_SCALE_THRESHOLD,
        baseline_quantile: float = DEFAULT_BASELINE_QUANTILE,
) -> Report:
    """Run the full detection pipeline over a parsed summary."""
    bin_baselines = bin_baseline_rates(summary, baseline_quantile)
    scales = superpartition_scales(summary, bin_baselines)
    cells, considered, skipped = flag_cells(
        summary,
        bin_baselines,
        scales,
        ratio_threshold=ratio_threshold,
        score_threshold=score_threshold,
        min_expected=min_expected,
    )
    return Report(
        rectangles=merge_cells(cells),
        superpartition_scales=flag_superpartition_scales(summary, scales, scale_threshold),
        n_bins_considered=considered,
        n_bins_skipped_empty=skipped,
        n_cells_flagged=len(cells),
    )


# ---------------------------------------------------------------------------
# Adjudication SQL
# ---------------------------------------------------------------------------


def adjudication_sql(
        rectangle: Rectangle,
        project_id: str,
        dataset_name: str,
        mode: str = 'variants',
        sample_table: str = 'sample_info',
        reference_schema: str = DEFAULT_REFERENCE_SCHEMA,
) -> str:
    """Generate the BigQuery query that settles whether a rectangle is a real dropout.

    The screen produces candidates; this query produces the answer.  If BigQuery holds
    rows in the window and the VDS does not, the data was lost between the two, which is
    proof rather than inference.  The tables are clustered by location, so restricting to
    the window keeps each check to a trivial scan.

    Both sides now cover the whole superpartition, so the counts are directly comparable
    rather than merely indicative.
    """
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    if reference_schema not in REFERENCE_SCHEMAS:
        raise ValueError(
            f"reference_schema must be one of {REFERENCE_SCHEMAS}, got {reference_schema!r}")

    table_index = f'{rectangle.superpartition:03d}'
    start_location = encode_location(rectangle.contig, rectangle.start)
    # bin_end is exclusive; the last covered position is one before it.
    end_location = encode_location(rectangle.contig, rectangle.end - 1)

    header = (
        f"-- Candidate dropout: {rectangle.contig}:{rectangle.start:,}-{rectangle.end:,} "
        f"superpartition {rectangle.superpartition}\n"
        f"-- VDS ({rectangle.n_samples:,} samples): observed "
        f"{rectangle.observed:,.0f} vs expected {rectangle.expected:,.0f} "
        f"(ratio {rectangle.ratio:.4f}, score {rectangle.score:,.1f})\n"
        f"-- Expect a non-trivial row count below if the data exists in BigQuery.\n"
        f"-- Add --dry_run to confirm the clustering prunes as expected before paying for it;\n"
        f"-- a scan far larger than the window means the predicate is not pruning.\n"
    )

    if mode == 'variants' or reference_schema == 'uncompressed':
        table = f'vet_{table_index}' if mode == 'variants' else f'ref_ranges_{table_index}'
        location_expr = 'v.location'
        # `location` is the clustering field for these tables, so this prunes.
        filter_clause = f'v.location BETWEEN {start_location} AND {end_location}'
        note = ''
    else:
        table = f'ref_ranges_{table_index}'
        location_expr = (
            f'{LOCATION_MULTIPLIER} * ((v.packed_ref_data >> {PACKED_CHROMOSOME_SHIFT}) & '
            f'0x{PACKED_CHROMOSOME_MASK:X}) + ((v.packed_ref_data >> {PACKED_POSITION_SHIFT}) & '
            f'0x{PACKED_POSITION_MASK:X})'
        )
        low, high = packed_ref_data_bounds(rectangle.contig, rectangle.start, rectangle.end - 1)
        # Filter on the packed value, not on the decoded location: packed_ref_data is the
        # clustering field, and a predicate on the decoded expression would prune nothing.
        filter_clause = f'v.packed_ref_data BETWEEN {low} AND {high}'
        note = (
            f"-- Compressed reference schema: filtering on packed_ref_data, the clustering\n"
            f"-- field, over the range encoding {rectangle.contig}:{rectangle.start:,}-"
            f"{rectangle.end - 1:,}.\n"
        )

    return (
        f"{header}{note}"
        f"SELECT\n"
        f"  COUNT(*) AS bq_rows,\n"
        f"  COUNT(DISTINCT v.sample_id) AS bq_samples,\n"
        f"  MIN({location_expr}) AS min_location,\n"
        f"  MAX({location_expr}) AS max_location\n"
        f"FROM `{project_id}.{dataset_name}.{table}` v\n"
        f"INNER JOIN `{project_id}.{dataset_name}.{sample_table}` s ON s.sample_id = v.sample_id\n"
        f"WHERE {filter_clause}\n"
        f"  AND s.withdrawn IS NULL\n"
        f"  AND s.is_control = false;\n"
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

REPORT_COLUMNS = (
    'contig', 'start', 'end', 'span', 'superpartition', 'n_bins',
    'n_samples', 'observed', 'expected', 'ratio', 'score',
)


def write_report(report: Report, handle: IO[str]) -> None:
    """Write flagged rectangles as a TSV."""
    handle.write('\t'.join(REPORT_COLUMNS) + '\n')
    for r in report.rectangles:
        handle.write('\t'.join([
            r.contig, str(r.start), str(r.end), str(r.span), str(r.superpartition),
            str(r.n_bins), str(r.n_samples), f'{r.observed:.0f}', f'{r.expected:.0f}',
            f'{r.ratio:.6f}', f'{r.score:.2f}',
        ]) + '\n')


def format_summary(report: Report, mode: str) -> str:
    """Human-readable digest for stdout."""
    lines = [
        f"Mode: {mode}",
        f"Bins considered: {report.n_bins_considered:,}",
        f"Bins skipped (no usable baseline): {report.n_bins_skipped_empty:,}",
        f"Cells flagged: {report.n_cells_flagged:,}",
        f"Rectangles after merging: {len(report.rectangles):,}",
    ]

    if report.superpartition_scales:
        lines.append('')
        lines.append('Superpartitions globally depleted relative to peers:')
        for finding in report.superpartition_scales:
            lines.append(
                f"  superpartition {finding.superpartition}: relative scale "
                f"{finding.relative_scale:.4f} over {finding.n_samples} sampled samples"
            )

    if report.rectangles:
        lines.append('')
        lines.append('Candidate dropouts, most severe first:')
        for r in report.rectangles:
            lines.append(
                f"  {r.contig}:{r.start:,}-{r.end:,} sp {r.superpartition} "
                f"({r.span:,} bp, {r.n_bins} bin(s)) observed {r.observed:,.0f} of "
                f"expected {r.expected:,.0f} -> {r.ratio * 100:.2f}% present, score {r.score:,.1f}"
            )
    elif not report.superpartition_scales:
        lines.append('')
        lines.append('No candidate dropouts found.')

    return '\n'.join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            'Detect rectangle data dropouts in a GVS VDS from the summary table emitted '
            'by vds_dropout_scan.py. Runs offline with no Hail dependency.'
        ),
    )
    parser.add_argument('--summary', required=True,
                        help='Path to the summary TSV (optionally gzipped).')
    parser.add_argument('--superpartitions', required=True,
                        help='Path to the superpartition/n_samples TSV (optionally gzipped).')
    parser.add_argument('--mode', choices=MODES, default='variants',
                        help='Which metric the summary holds. Default: variants.')
    parser.add_argument('--report-path', default=None,
                        help='Write flagged rectangles here as TSV. Omit for stdout digest only.')
    parser.add_argument('--sql-path', default=None,
                        help='Write BigQuery adjudication queries here. Requires --project-id '
                             'and --dataset-name.')
    parser.add_argument('--project-id', default=None,
                        help='BigQuery project for generated adjudication SQL.')
    parser.add_argument('--dataset-name', default=None,
                        help='BigQuery dataset for generated adjudication SQL.')
    parser.add_argument('--sample-table', default='sample_info',
                        help='Sample table or view used by adjudication SQL. Default: sample_info.')
    parser.add_argument('--reference-schema', choices=REFERENCE_SCHEMAS,
                        default=DEFAULT_REFERENCE_SCHEMA,
                        help='Which ref_ranges schema the dataset uses. AoU callsets use '
                             f'compressed. Default: {DEFAULT_REFERENCE_SCHEMA}.')
    parser.add_argument('--baseline-quantile', type=float, default=DEFAULT_BASELINE_QUANTILE,
                        help='Quantile across superpartitions used as the per-bin baseline. '
                             f'Default: {DEFAULT_BASELINE_QUANTILE}. Higher values stay robust '
                             'when a dropout hits many superpartitions at once.')
    parser.add_argument('--ratio-threshold', type=float, default=DEFAULT_RATIO_THRESHOLD,
                        help=f'Flag cells below this fraction of expected. '
                             f'Default: {DEFAULT_RATIO_THRESHOLD}.')
    parser.add_argument('--score-threshold', type=float, default=DEFAULT_SCORE_THRESHOLD,
                        help=f'Minimum depletion score to flag. Default: {DEFAULT_SCORE_THRESHOLD}.')
    parser.add_argument('--min-expected', type=float, default=DEFAULT_MIN_EXPECTED,
                        help=f'Skip cells expecting less than this much data. '
                             f'Default: {DEFAULT_MIN_EXPECTED}.')
    parser.add_argument('--scale-threshold', type=float,
                        default=DEFAULT_SUPERPARTITION_SCALE_THRESHOLD,
                        help='Flag superpartitions whose global relative scale falls below '
                             f'this. Default: {DEFAULT_SUPERPARTITION_SCALE_THRESHOLD}.')
    parser.add_argument('--max-sql-queries', type=int, default=DEFAULT_MAX_SQL_QUERIES,
                        help='Generate adjudication SQL for at most this many rectangles, '
                             'worst first. Every rectangle still appears in the report. '
                             f'Default: {DEFAULT_MAX_SQL_QUERIES}. Use 0 for no limit.')
    parser.add_argument('--fail-on-findings', action='store_true', default=False,
                        help='Exit non-zero when any candidate is found, for use as a gate.')

    args = parser.parse_args(argv)

    if args.sql_path and not (args.project_id and args.dataset_name):
        parser.error('--sql-path requires both --project-id and --dataset-name')

    n_samples = load_superpartitions(args.superpartitions)
    summary = load_summary(args.summary, n_samples)
    report = analyze(
        summary,
        ratio_threshold=args.ratio_threshold,
        score_threshold=args.score_threshold,
        min_expected=args.min_expected,
        scale_threshold=args.scale_threshold,
        baseline_quantile=args.baseline_quantile,
    )

    print(format_summary(report, args.mode))

    if args.report_path:
        with open(args.report_path, 'wt') as handle:
            write_report(report, handle)
        print(f"\nReport written to: {args.report_path}")

    if args.sql_path:
        limit = args.max_sql_queries or len(report.rectangles)
        selected = report.rectangles[:limit]
        with open(args.sql_path, 'wt') as handle:
            if not report.rectangles:
                handle.write('-- No candidate dropouts to adjudicate.\n')
            elif len(selected) < len(report.rectangles):
                omitted = len(report.rectangles) - len(selected)
                notice = (f'-- Showing the {len(selected)} most severe of '
                          f'{len(report.rectangles)} rectangles; {omitted} omitted.\n'
                          f'-- Raise --max-sql-queries to generate the rest. All '
                          f'{len(report.rectangles)} appear in the report.\n\n')
                handle.write(notice)
                print(f'\nNote: generated SQL for the {len(selected)} most severe of '
                      f'{len(report.rectangles)} rectangles; {omitted} omitted. '
                      'Raise --max-sql-queries to generate the rest.')
            for rectangle in selected:
                handle.write(adjudication_sql(
                    rectangle,
                    project_id=args.project_id,
                    dataset_name=args.dataset_name,
                    mode=args.mode,
                    sample_table=args.sample_table,
                    reference_schema=args.reference_schema,
                ))
                handle.write('\n')
        print(f"Adjudication SQL written to: {args.sql_path}")

    if args.fail_on_findings and not report.clean:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
