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
only names.  Two queries in ``scripts/variantstore/bq/`` produce the inputs:

``vds_dropout_sample_map.sql``
    The ``--sample-map-path`` input: every non-withdrawn, non-control sample, using the
    same filter GvsExtractAvroFilesForHail.wdl applies when exporting Avro.  Either tab-
    or comma-separated is accepted, so ``bq query --format=csv`` output works as-is.

``vds_dropout_sample_list.sql``
    The ``--sample-list-path`` input, reproducing ``stratified_sample`` below in SQL.
    ``materialize`` already writes this file, so the query is for inspecting the selection
    before spending cluster time, fixing the sample set before an r3 VDS exists, and
    cross-checking the Python.  That last one matters: the recorded list is what makes r2
    and r3 figures comparable, so two independent implementations agreeing is the evidence
    that reproducibility actually holds.  Output format and row order match what
    ``materialize`` writes, so the two can be diffed directly.

Samples in the map but absent from the VDS are ignored, so a map generated after a
withdrawal is still usable.  A VDS sample *missing* from the map is a hard error instead,
because silently excluding it would bias the superpartition totals.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import sys
import time
from collections import defaultdict
from typing import IO, Iterable, Mapping, Sequence

try:
    import hail as hl
except ModuleNotFoundError:  # pragma: no cover - exercised only off-cluster
    hl = None

DEFAULT_SUPERPARTITION_SIZE = 4000
DEFAULT_BIN_SIZE = 50_000
DEFAULT_SAMPLES_PER_SUPERPARTITION = 100
DEFAULT_STRIDE = 1
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
        # make the result depend on input order. This also makes the selection provably
        # identical to the SQL in bq/vds_dropout_sample_list.sql, whose ORDER BY needs an
        # explicit second key; see that file's header for the equivalence argument.
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
    least one scanned bin, at roughly ``1 / stride`` the cost -- while saying nothing
    about the unscanned remainder.
    """
    if stride < 1:
        raise ValueError(f"stride must be at least 1, got {stride}")
    return bin_index % stride == 0


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
    """Scale a probe's wall time up to a whole-genome estimate, in seconds."""
    if scanned_bases <= 0:
        raise ValueError('scanned_bases must be positive')
    return elapsed_seconds * (genome_length / scanned_bases)


def format_probe_report(
        elapsed_seconds: float,
        scanned_bases: int,
        n_cells: int,
        n_bins: int,
        n_samples: int | None = None,
) -> str:
    """Human-readable probe result, including the genome-wide extrapolation."""
    projected = extrapolate_runtime(elapsed_seconds, scanned_bases)
    lines = [
        'Probe result',
        f'  interval span:        {scanned_bases:,} bp',
        f'  bins produced:        {n_bins:,}',
        f'  non-zero cells:       {n_cells:,}',
    ]
    if n_samples is not None:
        lines.append(f'  samples screened:     {n_samples:,}')
    lines += [
        f'  elapsed:              {elapsed_seconds:,.1f} s',
        f'  genome-wide estimate: {projected / 3600:,.2f} h '
        f'({projected / 60:,.1f} min) at this throughput',
        '',
        'What this estimate covers: reading and aggregating, not writing. Against a',
        'full-width VDS it is therefore a floor for the materialize pass, not a full',
        'prediction of it. Against a narrow copy it estimates ongoing screening cost.',
        '',
        'Extrapolation assumes uniform density genome-wide, which is approximate: this',
        'interval is fully callable, while the real genome includes centromeres and other',
        'gaps that read faster than average.',
        '',
        'Check the Hail log for shuffle stages before committing to a genome-wide run.',
        'A single streaming pass is expected; a shuffle means the aggregation did not',
        'fuse, and the two-stage fallback should be used instead.',
    ]
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Hail plumbing
# ---------------------------------------------------------------------------


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


def read_sample_list(path: str) -> dict[str, int]:
    """Read a previously written sample list back as a name -> sample_id map.

    Reusing the recorded list is what keeps the r3 comparison honest: the headline
    figures are only comparable across VDSes if they cover the same samples.
    """
    result: dict[str, int] = {}
    with _open_read(path) as handle:
        for lineno, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith('sample_name\t'):
                continue
            fields = line.split('\t')
            if len(fields) < 2:
                raise ValueError(f"sample list line {lineno}: expected at least 2 columns")
            result[fields[0]] = int(fields[1])
    if not result:
        raise ValueError(f'{path} contained no samples')
    return result


def build_intervals(contigs: Sequence[str], interval_strings: Sequence[str]):
    """Build Hail locus intervals from explicit interval strings or whole contigs."""
    _require_hail()
    if interval_strings:
        return [hl.parse_locus_interval(s, reference_genome='GRCh38') for s in interval_strings]
    return [hl.parse_locus_interval(c, reference_genome='GRCh38') for c in contigs]


def interval_span(intervals) -> int:
    """Total base span of a list of Hail locus intervals."""
    _require_hail()
    total = 0
    for interval in intervals:
        start = hl.eval(interval.start.position)
        end = hl.eval(interval.end.position)
        total += max(0, end - start)
    return total


def read_and_subset_vds(vds_path: str, intervals):
    """Read a VDS, restricted to ``intervals`` when any are given."""
    _require_hail()
    if intervals:
        return hl.vds.read_vds(vds_path, intervals=intervals)
    return hl.vds.read_vds(vds_path)


def vds_sample_names(vds) -> list[str]:
    """Collect the VDS column keys.

    Reads only the column metadata, not entries, so this is cheap even on a full-width
    AoU VDS.
    """
    _require_hail()
    return vds.variant_data.cols().s.collect()


def _annotated_matrix(vds, mode: str, sample_to_superpartition: Mapping[str, int],
                      bin_size: int, stride: int):
    """Return the mode-appropriate matrix annotated with superpartition and bin."""
    _require_hail()
    matrix = vds.variant_data if mode == 'variants' else vds.reference_data

    lookup = hl.literal(dict(sample_to_superpartition), dtype=hl.tdict(hl.tstr, hl.tint32))
    matrix = matrix.annotate_cols(_sp=lookup[matrix.s])
    matrix = matrix.annotate_rows(
        _bin=(matrix.locus.position - 1) // bin_size,
    )
    if stride > 1:
        matrix = matrix.filter_rows(matrix._bin % stride == 0)
    return matrix


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


def aggregate_totals(vds, mode: str, sample_to_superpartition: Mapping[str, int],
                     bin_size: int, stride: int) -> dict[tuple[str, int], dict[int, float]]:
    """Aggregate a VDS into per-(bin, superpartition) totals.

    Expressed as one nested ``hl.agg.group_by`` inside a single ``aggregate_entries`` so
    the whole thing is a streaming pass with no shuffle: the outer grouping folds bins and
    the inner folds superpartitions, and only the finished matrix -- on the order of
    62,000 x 134 numbers for a genome-wide 50 kb scan -- returns to the driver.  Each
    partition's accumulator stays small because a partition spans few bins.
    """
    _require_hail()
    matrix = _annotated_matrix(vds, mode, sample_to_superpartition, bin_size, stride)

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


def per_sample_counts(vds, mode: str, samples: Sequence[str]) -> dict[str, float]:
    """Full-width per-sample totals over whatever the VDS is already restricted to."""
    _require_hail()
    matrix = vds.variant_data if mode == 'variants' else vds.reference_data
    matrix = matrix.filter_cols(hl.literal(set(samples)).contains(matrix.s))
    matrix = matrix.annotate_cols(_observed=_metric_expression(matrix, mode))
    collected = matrix.cols().select('_observed').collect()
    return {row.s: float(row._observed) for row in collected}


# ---------------------------------------------------------------------------
# Actions
# ---------------------------------------------------------------------------


def action_materialize(args) -> int:
    """Write a stratified narrow copy of the VDS, plus the sample and superpartition tables."""
    _require_hail()
    vds = hl.vds.read_vds(args.vds_path)
    sample_map = read_sample_map(args.sample_map_path)
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
    hl.vds.VariantDataset(reference_data, variant_data).write(
        args.output_path, overwrite=args.overwrite)
    print(f'Narrow VDS written to {args.output_path}')

    n_sp = write_lines(
        args.superpartitions_path, SUPERPARTITION_HEADER, format_superpartition_rows(chosen))
    n_samples = write_lines(
        args.sample_list_path, SAMPLE_LIST_HEADER, format_sample_list_rows(chosen, sample_map))
    print(f'Wrote {n_sp:,} superpartition rows to {args.superpartitions_path}')
    print(f'Wrote {n_samples:,} sample rows to {args.sample_list_path}')
    return 0


def _resolve_sample_assignment(args, vds) -> tuple[dict[str, int], dict[int, list[str]]]:
    """Work out each retained sample's superpartition, preferring a recorded sample list."""
    present = vds_sample_names(vds)
    if args.sample_list_path:
        sample_map = read_sample_list(args.sample_list_path)
        # The recorded list may name samples a later VDS no longer has; screen against
        # the intersection rather than failing, and report the shortfall.
        usable = [s for s in present if s in sample_map]
        if len(usable) < len(present):
            print(
                f'Note: {len(present) - len(usable):,} VDS sample(s) are absent from the '
                'recorded sample list and will not be screened'
            )
        assigned = {
            s: superpartition_for(sample_map[s], args.superpartition_size) for s in usable
        }
    else:
        sample_map = read_sample_map(args.sample_map_path)
        assigned = assign_superpartitions(sample_map, present, args.superpartition_size)

    grouped: dict[int, list[str]] = defaultdict(list)
    for name, superpartition in assigned.items():
        grouped[superpartition].append(name)
    return assigned, {sp: sorted(names) for sp, names in grouped.items()}


def action_scan(args) -> int:
    """Aggregate a VDS and write the summary and superpartition tables."""
    _require_hail()
    intervals = build_intervals(parse_contig_list(args.contigs),
                                parse_interval_list(args.intervals))
    vds = read_and_subset_vds(args.vds_path, intervals)
    assigned, grouped = _resolve_sample_assignment(args, vds)
    print(f'Screening {len(assigned):,} samples across {len(grouped):,} superpartitions')

    totals = aggregate_totals(vds, args.mode, assigned, args.bin_size, args.stride)
    rows = format_summary_rows(totals, args.bin_size)
    n_rows = write_lines(args.summary_path, SUMMARY_HEADER, rows)
    n_sp = write_lines(
        args.superpartitions_path, SUPERPARTITION_HEADER, format_superpartition_rows(grouped))

    print(f'Wrote {n_rows:,} summary rows ({len(totals):,} bins) to {args.summary_path}')
    print(f'Wrote {n_sp:,} superpartition rows to {args.superpartitions_path}')
    return 0


def action_probe(args) -> int:
    """Time a small scan and extrapolate to the whole genome."""
    _require_hail()
    interval_string = args.probe_interval
    intervals = build_intervals([], [interval_string])
    span = interval_span(intervals)
    print(f'Probing {interval_string} ({span:,} bp) in {args.mode} mode')

    vds = read_and_subset_vds(args.vds_path, intervals)
    assigned, _ = _resolve_sample_assignment(args, vds)

    started = time.monotonic()
    totals = aggregate_totals(vds, args.mode, assigned, args.bin_size, args.stride)
    elapsed = time.monotonic() - started

    n_cells = sum(len(v) for v in totals.values())
    print()
    print(format_probe_report(elapsed, span, n_cells, len(totals), len(assigned)))
    return 0


def action_full_depth(args) -> int:
    """Count per-sample data in flagged windows, at full width."""
    _require_hail()
    requested_intervals = parse_interval_list(args.intervals)
    if not requested_intervals:
        raise ValueError('--intervals is required for --action full-depth')
    target = {int(s) for s in args.target_superpartitions.split(',')} \
        if args.target_superpartitions else None

    sample_map = (read_sample_list(args.sample_list_path) if args.sample_list_path
                  else read_sample_map(args.sample_map_path))

    rows: list[str] = []
    for interval_string in requested_intervals:
        intervals = build_intervals([], [interval_string])
        vds = read_and_subset_vds(args.vds_path, intervals)
        present = [s for s in vds_sample_names(vds) if s in sample_map]
        assigned = {
            s: superpartition_for(sample_map[s], args.superpartition_size) for s in present
        }

        grouped: dict[int, list[str]] = defaultdict(list)
        for name, superpartition in assigned.items():
            if target is None or superpartition in target:
                grouped[superpartition].append(name)

        contig = hl.eval(intervals[0].start.contig)
        start = hl.eval(intervals[0].start.position)
        end = hl.eval(intervals[0].end.position)

        for superpartition in sorted(grouped):
            names = grouped[superpartition]
            counts = per_sample_counts(vds, args.mode, names)
            n_zero = sum(1 for name in names if counts.get(name, 0.0) <= 0)
            total = sum(counts.get(name, 0.0) for name in names)
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
