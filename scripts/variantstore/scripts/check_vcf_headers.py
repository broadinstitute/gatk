# -*- coding: utf-8 -*-
"""VS-1966: fail-fast validation of ingested VCF headers.

Runs after the *headers-only* GVS ingest phase (``load_vcf_headers = true``) and before the
expensive vet/ref data ingest, as an early sanity check of the input gVCFs. It reads the header
tables populated by ingest -- ``vcf_header_lines`` and ``sample_vcf_header`` -- joined against
``sample_info``, and validates cohort-wide consistency so a bad cohort is caught before any
vet/ref compute is spent. This automates the manual DRAGEN-version query in ``AOU_DELIVERABLES.md``
and implements the section 8 checks of the header-loading design doc
(``scripts/variantstore/docs/parquet/header_loading_design.md``).

Header model recap (``CreateVariantIngestFiles.buildAllVcfLineHeaders``): each header line whose
key contains "CommandLine" is stored individually with ``is_expected_unique = TRUE`` (varies per
sample); all remaining lines are joined into one ``is_expected_unique = FALSE`` blob (shared across
samples). So both signals we validate live in ``is_expected_unique = TRUE`` chunks:

  * reblocking  -> ``##GATKCommandLine=<ID=ReblockGVCF, ...`` (the ``##source=ReblockGVCF`` line, by
    contrast, is a non-CommandLine line and falls into the shared FALSE blob, so it is not a
    per-sample signal);
  * DRAGEN version -> ``##DRAGENCommandLine=<ID=dragen, ...Version="SW: ...3.7.8">`` (filtered on
    ``ID=dragen,`` to skip the ``ID=HashTableBuild`` DRAGENCommandLine, which carries a different
    SW version).

Checks (all fatal unless noted):
  * Integrity: every non-control, non-withdrawn sample in ``sample_info`` has header data
    (>= 1 chunk) and >= 1 ``is_expected_unique = TRUE`` chunk; no orphan hashes.
  * Reblocking: every such sample has an ``is_expected_unique = TRUE`` chunk mentioning
    ``ReblockGVCF``.
  * DRAGEN version: all samples share a single DRAGEN version triplet (e.g. ``3.7.8``); if
    ``--expected_dragen_version`` is given as a single value every sample must match it exactly
    (AoU), or as a range every sample's triplet must fall within it. Ranges accept interval
    notation -- ``3.4.12-3.7.8`` (both inclusive), ``[3.7.8-3.8)`` (inclusive-exclusive),
    ``(3.7-3.8)`` (both exclusive). A DRAGEN command line whose SW version cannot be reduced to a
    numeric triplet fails the check in every mode -- such rows are never silently dropped. A *mixed*
    cohort, where some non-control samples carry a DRAGEN command line and some do not, also fails --
    even with no ``--expected_dragen_version`` -- because mixing DRAGEN and non-DRAGEN provenance
    risks batch effects. A non-DRAGEN sample yields no breakdown row, so the shortfall is caught
    by cross-referencing the expected cohort count. A cohort with *no* DRAGEN command lines at all
    is not mixed, so it is fine and reported informationally.
  * Shared-blob distribution (informational): how many distinct ``is_expected_unique = FALSE``
    blobs there are and how many samples carry each (> 1 indicates distinct delivery batches).

By default every check covers the whole dataset: all non-control, non-withdrawn samples in
``sample_info``. ``--sample_names_file`` restricts them instead to the named samples -- the current
ingest batch, as ``GvsBulkIngestGenomes`` passes it -- so an incremental ingest is judged on its own
samples rather than being failed by a pre-existing problem in an earlier batch. The names are staged
in a short-lived BigQuery table and semi-joined into each query, which keeps the restriction working
for AoU-sized batches; an inline ``IN`` list of a few hundred thousand names would blow BigQuery's
1 MB query-text limit. Restricting also adds a fatal check of its own: every requested name must
actually be an eligible sample in ``sample_info``, so a name that ingest never assigned an id fails
rather than silently shrinking the cohort under test.

Whether a failing result aborts the pipeline is decided by the caller
(``GvsValidateVcfHeaders.wdl`` ``fail_on_validation_errors``, default true); this script always
computes an overall pass/fail and writes a human-readable report.
"""
import argparse
import contextlib
import datetime
import os
import re
import tempfile
import uuid
from collections import namedtuple

# add labels for DSP Cloud Cost Control Labeling and Reporting
query_labels_map = {'service': 'gvs', 'team': 'variants', 'managedby': 'gvs_validate_vcf_headers'}

# Number of offending sample names to include in the report per check (the full lists can be huge).
EXAMPLE_LIMIT = 20

# The SW version string in a DRAGEN command line, e.g. "SW: 05.021.604.3.7.8". A '.' inside a
# character class is a literal, so no escaping is needed.
SW_VERSION_REGEX = r'SW: [0-9.]+'

# name:  short identifier for the check
# passed: bool -- did the check pass
# fatal: bool -- does a failure of this check make the overall validation fail
# lines: list[str] -- human-readable detail for the report
CheckResult = namedtuple('CheckResult', ['name', 'passed', 'fatal', 'lines'])

# A parsed DRAGEN version range. low/high are int-tuple version keys (e.g. (3, 7, 8)); the
# *_inclusive flags come from interval notation -- '[' / ']' are inclusive, '(' / ')' exclusive.
VersionRange = namedtuple('VersionRange', ['low', 'low_inclusive', 'high', 'high_inclusive'])


# --- SQL builders (pure functions, no BigQuery client needed -- unit-testable) ---------------------

def _sample_name_restriction(sample_names_table, column, indent=14):
    """SQL fragment restricting ``column`` (a sample_name expression) to the staged batch.

    Returns '' when ``sample_names_table`` is None, so the unrestricted whole-dataset query is
    byte-for-byte what it was before this option existed. A semi-join is used rather than an inline
    ``IN`` list because the batch can hold hundreds of thousands of names; see
    ``stage_sample_names_table``. ``indent`` only lines the fragment up with the WHERE clause it is
    appended to.
    """
    if sample_names_table is None:
        return ""
    _validate_bq_table_path(sample_names_table)
    return f"\n{' ' * indent}AND {column} IN (SELECT sample_name FROM `{sample_names_table}`)"


def _sample_id_restriction(project_id, dataset_name, sample_names_table, column):
    """As ``_sample_name_restriction``, but for queries that have only a sample_id to filter on.

    ``sample_vcf_header`` carries no sample_name, so the batch's names are resolved through
    ``sample_info``. The is_control / withdrawn filter is applied here too, matching the cohort the
    name-keyed queries use.
    """
    if sample_names_table is None:
        return ""
    _validate_bq_table_path(sample_names_table)
    return f"""
          AND {column} IN (
              SELECT si.sample_id
              FROM `{project_id}.{dataset_name}.sample_info` si
              WHERE si.is_control = FALSE AND si.withdrawn IS NULL
                AND si.sample_name IN (SELECT sample_name FROM `{sample_names_table}`)
          )"""


def per_sample_summary_sql(project_id, dataset_name, sample_names_table=None):
    """One-pass per-sample rollup over the expected (non-control, non-withdrawn) cohort.

    Returns a single summary row: the expected sample count, how many have header data, and the
    counts + example sample names for each fatal integrity / reblocking condition. ``COUNTIF`` over
    a LEFT JOIN treats a sample with no header rows as all-zero counts, so missing-header samples
    surface as ``chunk_count = 0``. Each example array is capped at ``EXAMPLE_LIMIT + 1`` -- one more
    than we display -- so the report can tell "exactly N offenders" from "more than N" (see
    ``_examples``).

    When ``sample_names_table`` is given the cohort is narrowed to the samples named in it (the
    current ingest batch); otherwise it is the whole dataset.
    """
    return f"""
        WITH expected AS (
            SELECT sample_id, sample_name
            FROM `{project_id}.{dataset_name}.sample_info`
            WHERE is_control = FALSE AND withdrawn IS NULL{_sample_name_restriction(sample_names_table, 'sample_name')}
        ),
        per_sample AS (
            SELECT
                e.sample_id,
                e.sample_name,
                COUNTIF(vhl.vcf_header_lines_hash IS NOT NULL) AS chunk_count,
                COUNTIF(vhl.is_expected_unique) AS unique_chunk_count,
                COUNTIF(vhl.is_expected_unique
                        AND CONTAINS_SUBSTR(vhl.vcf_header_lines, 'ReblockGVCF')) AS reblock_chunk_count
            FROM expected e
            LEFT JOIN `{project_id}.{dataset_name}.sample_vcf_header` svh ON svh.sample_id = e.sample_id
            LEFT JOIN `{project_id}.{dataset_name}.vcf_header_lines` vhl
                ON vhl.vcf_header_lines_hash = svh.vcf_header_lines_hash
            GROUP BY e.sample_id, e.sample_name
        )
        SELECT
            COUNT(*) AS expected_samples,
            COUNTIF(chunk_count > 0) AS samples_with_headers,
            COUNTIF(chunk_count = 0) AS samples_missing_headers,
            COUNTIF(unique_chunk_count = 0) AS samples_missing_unique_chunk,
            COUNTIF(reblock_chunk_count = 0) AS samples_not_reblocked,
            ARRAY_AGG(IF(chunk_count = 0, sample_name, NULL) IGNORE NULLS LIMIT {EXAMPLE_LIMIT + 1}) AS example_missing_headers,
            ARRAY_AGG(IF(unique_chunk_count = 0, sample_name, NULL) IGNORE NULLS LIMIT {EXAMPLE_LIMIT + 1}) AS example_missing_unique_chunk,
            ARRAY_AGG(IF(reblock_chunk_count = 0, sample_name, NULL) IGNORE NULLS LIMIT {EXAMPLE_LIMIT + 1}) AS example_not_reblocked
        FROM per_sample
    """


def dragen_version_breakdown_sql(project_id, dataset_name, sample_names_table=None):
    """Per full-SW-version sample counts, restricted to ``ID=dragen`` command lines.

    Mirrors the manual AoU query (``AOU_DELIVERABLES.md``) but joins ``sample_info`` to restrict to
    non-control, non-withdrawn samples and counts distinct samples per version. The Python side
    reduces each full version to its final triplet for the consistency check.
    """
    return f"""
        SELECT
            REGEXP_EXTRACT(vhl.vcf_header_lines, r'{SW_VERSION_REGEX}') AS sw_version,
            COUNT(DISTINCT svh.sample_id) AS n_samples
        FROM `{project_id}.{dataset_name}.vcf_header_lines` vhl
        JOIN `{project_id}.{dataset_name}.sample_vcf_header` svh USING (vcf_header_lines_hash)
        JOIN `{project_id}.{dataset_name}.sample_info` si ON si.sample_id = svh.sample_id
        WHERE vhl.is_expected_unique = TRUE
          AND CONTAINS_SUBSTR(vhl.vcf_header_lines, 'DRAGENCommandLine=<ID=dragen,')
          AND si.is_control = FALSE AND si.withdrawn IS NULL{_sample_name_restriction(sample_names_table, 'si.sample_name', indent=10)}
        GROUP BY sw_version
        ORDER BY n_samples DESC
    """


def orphan_hash_sql(project_id, dataset_name, sample_names_table=None):
    """Referential integrity: associations in ``sample_vcf_header`` whose hash has no
    ``vcf_header_lines`` row. Should always be zero for a correct load.

    Restricted to the batch when ``sample_names_table`` is given: an orphan left behind by an
    earlier batch is not this ingest's to fail on."""
    return f"""
        SELECT
            COUNT(*) AS orphan_associations,
            COUNT(DISTINCT svh.sample_id) AS affected_samples
        FROM `{project_id}.{dataset_name}.sample_vcf_header` svh
        LEFT JOIN `{project_id}.{dataset_name}.vcf_header_lines` vhl USING (vcf_header_lines_hash)
        WHERE vhl.vcf_header_lines_hash IS NULL{_sample_id_restriction(project_id, dataset_name, sample_names_table, 'svh.sample_id')}
    """


def shared_blob_distribution_sql(project_id, dataset_name, sample_names_table=None):
    """Distinct shared (``is_expected_unique = FALSE``) blobs and how many samples carry each.

    Restricted to the batch when ``sample_names_table`` is given, so the blob count reports the
    delivery batches within *this* ingest rather than across the dataset's whole history."""
    return f"""
        SELECT
            svh.vcf_header_lines_hash AS blob_hash,
            COUNT(DISTINCT svh.sample_id) AS n_samples
        FROM `{project_id}.{dataset_name}.sample_vcf_header` svh
        JOIN `{project_id}.{dataset_name}.vcf_header_lines` vhl USING (vcf_header_lines_hash)
        WHERE vhl.is_expected_unique = FALSE{_sample_id_restriction(project_id, dataset_name, sample_names_table, 'svh.sample_id')}
        GROUP BY blob_hash
        ORDER BY n_samples DESC
    """


# --- Pure evaluation logic (no BigQuery client needed -- unit-testable) ----------------------------

def dragen_version_triplet(sw_version):
    """Reduce a DRAGEN SW version string to its final numeric triplet.

    AoU cares that the *last triplet* of the version is the expected one (e.g. ``3.7.8``): the
    leading components mix hardware / build identifiers that legitimately differ between delivery
    batches. E.g. both ``SW: 05.021.604.3.7.8`` and ``SW: 07.021.604.3.7.8`` -> ``3.7.8``.
    Returns ``None`` if fewer than three numeric components are present.
    """
    if sw_version is None:
        return None
    nums = re.findall(r'[0-9]+', sw_version)
    if len(nums) < 3:
        return None
    return '.'.join(nums[-3:])


def _version_key(version):
    """Turn a dotted numeric version like ``3.7.8`` into a tuple of ints for ordering.

    Returns ``None`` if there are no numeric components. Comparisons use tuple ordering over exactly
    the components present, so bounds are matched *positionally*: a short bound is NOT shorthand for
    "all of that minor". An inclusive upper bound of ``3.8`` -> ``(3, 8)`` excludes ``3.8.1`` ->
    ``(3, 8, 1)`` because ``(3, 8, 1) > (3, 8)``. For "everything below 3.8" use the exclusive form
    ``...-3.8)``; for "through all of 3.8.x" use an exclusive next-minor bound, ``...-3.9)``.
    """
    if version is None:
        return None
    nums = re.findall(r'[0-9]+', version)
    if not nums:
        return None
    return tuple(int(n) for n in nums)


# A well-formed, user-supplied version token: plain dotted-numeric, e.g. '3', '3.7', '3.7.8'.
_STRICT_VERSION_RE = re.compile(r'^[0-9]+(?:\.[0-9]+)*$')


def _require_version(token, original):
    """Strictly validate a *user-supplied* version token and return its int-tuple key.

    Unlike ``_version_key`` (lenient ``findall``, used on trusted sample data), this rejects
    anything that is not a plain dotted-numeric version -- e.g. '3.7,8', '3.7.8x', embedded
    spaces -- so a typo in ``--expected_dragen_version`` fails loudly here rather than silently
    never matching any sample. ``original`` is the full spec, quoted in the error for context.
    """
    t = token.strip()
    if not _STRICT_VERSION_RE.match(t):
        raise ValueError(
            f"invalid DRAGEN version '{original}': '{t}' is not a dotted-numeric version (e.g. '3.7.8')")
    return tuple(int(n) for n in t.split('.'))


def parse_expected_dragen_spec(expected_dragen_version):
    """Interpret the ``--expected_dragen_version`` value as either an exact version or a range.

    A single value (no ``-``) means an exact triplet match, as AoU requires. A value of the form
    ``LOW-HIGH`` means a range; interval-notation brackets set each end's inclusivity (``[``/``]``
    inclusive, ``(``/``)`` exclusive), and a bare bound defaults to inclusive. Examples:
      * ``3.4.12-3.7.8``  -> low and high both inclusive;
      * ``[3.7.8-3.8)``   -> low inclusive, high exclusive (3.7.8 <= t < 3.8);
      * ``(3.7-3.8)``     -> low and high both exclusive.
    In range mode samples may legitimately span multiple triplets inside the range.

    Returns one of:
      * ``(None, None)``            -- no expected version supplied;
      * ``('exact', 'X.Y.Z')``      -- exact-match mode;
      * ``('range', VersionRange)`` -- range mode.

    Raises ``ValueError`` if the value is malformed (bracket without a range, wrong shape,
    a bound that is not a plain dotted-numeric version -- e.g. a typo like '3.7,8' or stray
    characters/whitespace -- or an empty/inverted interval). Versions are dot-separated digits, so
    the ``-`` separator and the bracket characters are unambiguous.
    """
    if not expected_dragen_version:
        return None, None
    spec = expected_dragen_version.strip()
    if not spec:  # whitespace-only -- treat as unspecified rather than a never-matching exact value
        return None, None

    has_bracket = any(c in '[]()' for c in spec)
    if '-' not in spec:
        if has_bracket:
            raise ValueError(
                f"invalid DRAGEN version '{expected_dragen_version}': looks like an interval but has "
                f"no 'LOW-HIGH' separator, e.g. '[3.7.8-3.8)'")
        _require_version(spec, expected_dragen_version)  # reject typos (e.g. '3.7,8') loudly
        return 'exact', spec

    # Peel off optional interval-notation brackets to determine each end's inclusivity.
    low_inclusive = True
    high_inclusive = True
    body = spec
    if body[0] in '[(':
        low_inclusive = body[0] == '['
        body = body[1:]
    if body and body[-1] in '])':
        high_inclusive = body[-1] == ']'
        body = body[:-1]

    parts = body.split('-')
    if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
        raise ValueError(
            f"invalid DRAGEN version range '{expected_dragen_version}' "
            f"(expected 'LOW-HIGH', e.g. '3.4.12-3.7.8' or '[3.7.8-3.8)')")
    low = _require_version(parts[0], expected_dragen_version)
    high = _require_version(parts[1], expected_dragen_version)
    if low > high or (low == high and not (low_inclusive and high_inclusive)):
        raise ValueError(
            f"invalid DRAGEN version range '{expected_dragen_version}': empty or inverted interval")
    return 'range', VersionRange(low, low_inclusive, high, high_inclusive)


def _examples(names):
    """Render a (possibly truncated) list of example sample names for the report.

    The SQL fetches one more name than ``EXAMPLE_LIMIT`` (see ``per_sample_summary_sql``), so
    ``len(names) > EXAMPLE_LIMIT`` distinguishes "more than EXAMPLE_LIMIT offenders" from "exactly
    EXAMPLE_LIMIT"; either way only the first ``EXAMPLE_LIMIT`` are shown.
    """
    names = list(names or [])
    if not names:
        return ""
    shown = ", ".join(names[:EXAMPLE_LIMIT])
    suffix = f" (showing first {EXAMPLE_LIMIT})" if len(names) > EXAMPLE_LIMIT else ""
    return f"    e.g. {shown}{suffix}"


def evaluate_integrity(summary, requested_samples=None):
    """Integrity checks derived from the per-sample summary row. Returns a list of CheckResult.

    ``requested_samples``, when the run was restricted to a sample names file, is how many distinct
    names that file held. The cohort the queries found is compared against it so a name that never
    became an eligible ``sample_info`` row fails the run instead of quietly shrinking the cohort
    every other check is measured over.
    """
    expected = summary.get('expected_samples') or 0
    with_headers = summary.get('samples_with_headers') or 0
    missing_headers = summary.get('samples_missing_headers') or 0
    missing_unique = summary.get('samples_missing_unique_chunk') or 0

    results = []

    # There must be samples to validate at all.
    scope = "in sample_info" if requested_samples is None else "matching the requested sample names"
    results.append(CheckResult(
        name='samples_present',
        passed=expected > 0,
        fatal=True,
        lines=[f"Non-control, non-withdrawn samples {scope}: {expected}"] +
              ([] if expected > 0 else ["    No samples to validate -- did the header ingest run against the right dataset?"]),
    ))

    # When restricted to a sample names file, every requested name must be an eligible sample.
    if requested_samples is not None:
        unmatched = requested_samples - expected
        lines = [f"Requested sample names: {requested_samples}",
                 f"Matched as non-control, non-withdrawn samples in sample_info: {expected}"]
        if unmatched > 0:
            lines.append(f"    {unmatched} requested name(s) did not match an eligible sample -- "
                         "not assigned a sample id, or marked control / withdrawn.")
        results.append(CheckResult(
            name='requested_samples_present',
            passed=unmatched <= 0,
            fatal=True,
            lines=lines,
        ))

    # Every expected sample must have header data.
    lines = [f"Samples with header data: {with_headers} / {expected}",
             f"Samples missing header data: {missing_headers}"]
    ex = _examples(summary.get('example_missing_headers'))
    if ex:
        lines.append(ex)
    results.append(CheckResult(
        name='all_samples_have_headers',
        passed=missing_headers == 0,
        fatal=True,
        lines=lines,
    ))

    # Every expected sample must have at least one per-sample (is_expected_unique) chunk.
    lines = [f"Samples missing an is_expected_unique chunk: {missing_unique}"]
    ex = _examples(summary.get('example_missing_unique_chunk'))
    if ex:
        lines.append(ex)
    results.append(CheckResult(
        name='all_samples_have_unique_chunk',
        passed=missing_unique == 0,
        fatal=True,
        lines=lines,
    ))

    return results


def evaluate_referential_integrity(orphan):
    orphan_associations = orphan.get('orphan_associations') or 0
    affected = orphan.get('affected_samples') or 0
    lines = [f"Orphan sample->hash associations (no matching vcf_header_lines row): {orphan_associations}"]
    if orphan_associations:
        lines.append(f"    affected samples: {affected}")
    return CheckResult(
        name='referential_integrity',
        passed=orphan_associations == 0,
        fatal=True,
        lines=lines,
    )


def evaluate_reblocking(summary, require_reblocking=True):
    """Reblocking check. Fatal by default (AoU wants to fail fast on non-reblocked input); pass
    ``require_reblocking=False`` to downgrade it to informational for cohorts where reblocking is
    not required."""
    expected = summary.get('expected_samples') or 0
    not_reblocked = summary.get('samples_not_reblocked') or 0
    lines = [f"Samples with a reblocking (ReblockGVCF) command line: {expected - not_reblocked} / {expected}",
             f"Samples NOT reblocked: {not_reblocked}"]
    ex = _examples(summary.get('example_not_reblocked'))
    if ex:
        lines.append(ex)
    if not require_reblocking:
        lines.append("    (reblocking not required for this run; reported for information only)")
    return CheckResult(
        name='reblocking',
        passed=not_reblocked == 0,
        fatal=require_reblocking,
        lines=lines,
    )


def evaluate_dragen_version(dragen_rows, expected_dragen_version=None, expected_samples=None):
    """Consistency (and optionally an exact match or inclusive range) of the DRAGEN version.

    ``dragen_rows``: iterable of rows with ``sw_version`` and ``n_samples`` -- one row per distinct
    SW version, covering only samples that carry a ``DRAGENCommandLine=<ID=dragen, ...>`` line.
    ``expected_dragen_version`` accepts either a single triplet (``'3.7.8'`` -- exact match, as AoU
    requires) or a range with optional interval-notation brackets (``'3.4.12-3.7.8'``,
    ``'[3.7.8-3.8)'``, ``'(3.7-3.8)'``); see ``parse_expected_dragen_spec``. ``expected_samples`` is
    the size of the expected (non-control, non-withdrawn) cohort; it is compared against the number
    of samples that actually carry a DRAGEN line so a partial cohort cannot slip through.
    Failures:
      * expected version given but no DRAGEN command lines found;
      * a mixed cohort -- some samples carry a DRAGEN command line and some do not -- which is fatal
        even with no expected version, because mixing DRAGEN and non-DRAGEN provenance risks batch
        effects. A non-DRAGEN sample yields no breakdown row, so this is caught by comparing the
        DRAGEN sample count against ``expected_samples``;
      * any DRAGEN command line whose SW version cannot be reduced to a numeric triplet (fatal in
        every mode -- these rows are never silently dropped from the comparison);
      * exact / consistency-only mode: more than one distinct version triplet across the cohort;
      * exact mode: any triplet differs from the expected version;
      * range mode: any triplet falls outside the interval (multiple triplets inside the range are
        allowed -- relaxing a range but still requiring a single triplet would be pointless);
      * ``expected_dragen_version`` is a malformed range.
    A cohort with NO DRAGEN command lines and no expected version is fine (a legitimately non-DRAGEN
    cohort) and reported informationally; a cohort where only SOME samples are DRAGEN is a mixed
    cohort and always fails.
    """
    # Aggregate sample counts by full SW version and by triplet.
    version_counts = {}
    triplet_counts = {}
    for row in dragen_rows:
        sw = row.get('sw_version')
        n = row.get('n_samples') or 0
        version_counts[sw] = version_counts.get(sw, 0) + n
        triplet = dragen_version_triplet(sw)
        triplet_counts[triplet] = triplet_counts.get(triplet, 0) + n

    lines = ["DRAGEN version breakdown (SW string -> sample count):"]
    if version_counts:
        for sw, n in sorted(version_counts.items(), key=lambda kv: (-kv[1], str(kv[0]))):
            lines.append(f"    {sw} (triplet {dragen_version_triplet(sw)}): {n}")
    else:
        lines.append("    (no DRAGENCommandLine=<ID=dragen, ...> lines found)")

    # Interpret the expected-version spec up front; a malformed range is itself a fatal failure.
    try:
        kind, payload = parse_expected_dragen_spec(expected_dragen_version)
    except ValueError as e:
        # Name the offending input so the WDL log points the operator straight at what to fix.
        return CheckResult(name='dragen_version', passed=False, fatal=True,
                           lines=lines + [f"    FAIL: malformed --expected_dragen_version: {e}"])

    distinct_triplets = sorted(t for t in triplet_counts if t is not None)
    # DRAGEN rows whose SW string we could not reduce to a numeric triplet: fewer than three numeric
    # components, or a NULL sw_version -- dragen_version_breakdown_sql matches ID=dragen command lines
    # but extracts the 'SW: X.Y.Z' field, which a differently-shaped DRAGEN line can leave NULL.
    unparseable = sum(n for t, n in triplet_counts.items() if t is None)

    if not version_counts:
        if kind is not None:
            return CheckResult(
                name='dragen_version',
                passed=False,
                fatal=True,
                lines=lines + [f"    FAIL: expected DRAGEN version {expected_dragen_version} but no DRAGEN command lines were found"],
            )
        return CheckResult(
            name='dragen_version',
            passed=True,
            fatal=False,  # informational: nothing to check for a non-DRAGEN cohort
            lines=lines + ["    No DRAGEN command lines found; skipping DRAGEN version check"],
        )

    # Invariant past the guard above: version_counts is non-empty, and every DRAGEN row (a GROUP BY
    # row, so n_samples >= 1) lands in exactly one of distinct_triplets (parseable) or `unparseable`
    # (not). The two therefore cannot both be empty, so the ``> 1`` and ``elif distinct_triplets``
    # guards below never leave the *returned rows* silently unchecked. Samples with no DRAGEN row at
    # all are invisible to those rows; they are handled by the cohort-coverage check just below.
    passed = True

    # Cohort coverage: dragen_version_breakdown_sql returns a row only for samples that
    # carry a DRAGENCommandLine=<ID=dragen, ...> line, so a sample with no DRAGEN line contributes no
    # row and is invisible to the triplet checks. Cross-reference the number of DRAGEN samples against
    # the full expected cohort. We are past the `not version_counts` guard, so at least one sample IS
    # DRAGEN; any shortfall here therefore means a MIXED cohort (some DRAGEN, some not). A mixed
    # cohort is a hard error even with no expected version, because silently mixing DRAGEN and
    # non-DRAGEN provenance risks batch effects -- the operator must split the cohort by provenance
    # or otherwise confirm intent. (A cohort with no DRAGEN lines at all is not mixed; it was
    # handled by the informational path above.) This may be revisited later (e.g. an explicit
    # opt-in flag).
    dragen_sample_count = sum(version_counts.values())
    if expected_samples is not None and dragen_sample_count < expected_samples:
        missing = expected_samples - dragen_sample_count
        passed = False
        lines.append(f"Samples with a DRAGEN command line: {dragen_sample_count} / {expected_samples}")
        if kind is not None:
            lines.append(f"    FAIL: {missing} non-control sample(s) have no DRAGEN command line, but "
                         f"--expected_dragen_version {expected_dragen_version} requires every sample to be DRAGEN")
        else:
            lines.append(f"    FAIL: mixed cohort -- {missing} of {expected_samples} non-control "
                         f"sample(s) have no DRAGEN command line while others do; mixing DRAGEN and "
                         f"non-DRAGEN samples is not allowed (batch effects). Split the cohort by "
                         f"provenance.")

    # An unparseable DRAGEN version is a hole in the check, not a pass. Filtering these rows out of the
    # comparison would let real version drift (or a NULL SW field) slip through the very check meant to
    # catch it -- in exact mode a handful of matching samples could otherwise mask thousands of
    # unreadable ones. Fatal in every mode.
    if unparseable:
        passed = False
        lines.append(f"    FAIL: {unparseable} sample(s) have a DRAGEN command line whose SW version "
                     f"could not be parsed to a numeric triplet")

    if kind == 'range':
        # Range mode: every triplet must lie within the interval (inclusivity per the brackets);
        # multiple triplets are fine.
        r = payload

        def _in_range(key):
            low_ok = key >= r.low if r.low_inclusive else key > r.low
            high_ok = key <= r.high if r.high_inclusive else key < r.high
            return low_ok and high_ok

        out_of_range = [t for t in distinct_triplets if not _in_range(_version_key(t))]
        if out_of_range:
            passed = False
            lines.append(f"    FAIL: expected DRAGEN version within range {expected_dragen_version}, "
                         f"out-of-range triplet(s): {', '.join(out_of_range)}")
        elif distinct_triplets:  # empty means every row was unparseable (already failed above)
            lines.append(f"    OK: all triplet(s) within range {expected_dragen_version}: "
                         f"{', '.join(distinct_triplets)}")
    else:
        # Exact or consistency-only mode: all samples must share a single triplet. An empty
        # distinct_triplets means every row was unparseable (already failed above), so guard against
        # emitting the 'span multiple triplets:' line with nothing after the colon.
        if len(distinct_triplets) > 1:
            passed = False
            lines.append(f"    FAIL: samples span multiple DRAGEN version triplets: {', '.join(distinct_triplets)}")
        if kind == 'exact':
            mismatched = [t for t in distinct_triplets if t != payload]
            if mismatched:
                passed = False
                lines.append(f"    FAIL: expected DRAGEN version triplet {payload}, also saw: {', '.join(mismatched)}")

    return CheckResult(name='dragen_version', passed=passed, fatal=True, lines=lines)


def evaluate_shared_blob_distribution(blob_rows):
    """Informational: distinct shared blobs and their sample counts."""
    rows = list(blob_rows)
    lines = [f"Distinct shared (is_expected_unique = FALSE) blobs: {len(rows)}"]
    if len(rows) > 1:
        lines.append("    (more than one shared blob usually means samples came from distinct delivery batches)")
    for row in rows[:EXAMPLE_LIMIT]:
        lines.append(f"    {row.get('blob_hash')}: {row.get('n_samples')} samples")
    return CheckResult(name='shared_blob_distribution', passed=True, fatal=False, lines=lines)


def evaluate(summary, dragen_rows, orphan, blob_rows, expected_dragen_version=None,
             require_reblocking=True, requested_samples=None):
    """Run every check and return ``(overall_passed, [CheckResult, ...])``.

    ``overall_passed`` is true iff every *fatal* check passed; informational checks never affect it.
    ``requested_samples`` is the size of the requested batch when the run was restricted to a sample
    names file, and None otherwise; see ``evaluate_integrity``.
    """
    checks = []
    checks.extend(evaluate_integrity(summary, requested_samples))
    checks.append(evaluate_reblocking(summary, require_reblocking))
    checks.append(evaluate_dragen_version(dragen_rows, expected_dragen_version,
                                          summary.get('expected_samples')))
    checks.append(evaluate_referential_integrity(orphan))
    checks.append(evaluate_shared_blob_distribution(blob_rows))

    overall_passed = all(c.passed for c in checks if c.fatal)
    return overall_passed, checks


def compose_report(project_id, dataset_name, expected_dragen_version, overall_passed, checks,
                   requested_samples=None):
    """Render the human-readable report."""
    scope = ("whole dataset (all non-control, non-withdrawn samples)" if requested_samples is None
             else f"{requested_samples} requested sample name(s)")
    header = [
        "GVS VCF Header Validation Report",
        "==========================================",
        f"Dataset: {project_id}.{dataset_name}",
        f"Scope: {scope}",
        f"Expected DRAGEN version: {expected_dragen_version or '(not specified -- consistency only)'}",
        f"OVERALL: {'PASS' if overall_passed else 'FAIL'}",
        "",
    ]
    body = []
    for c in checks:
        if not c.fatal:
            status = "INFO"
        else:
            status = "PASS" if c.passed else "FAIL"
        body.append(f"[{status}] {c.name}")
        body.extend(f"  {line}" for line in c.lines)
        body.append("")
    return "\n".join(header + body)


# --- Execution (imports the BigQuery client lazily so the SQL builders / evaluators stay
#     dependency-free for unit tests) -----------------------------------------------------------

# project_id and dataset_name are interpolated directly into the SQL builders, so reject anything
# outside a conservative allowlist before they reach a query (guards against malformed identifiers
# and SQL injection, mirroring _validate_table_prefixes in parse_and_group_files.py). GCP project
# ids may contain letters, digits, hyphens; BigQuery dataset names allow only letters, digits, and
# underscores.
_PROJECT_ID_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*$')
_DATASET_NAME_RE = re.compile(r'^[A-Za-z0-9_]+$')
# The staged sample-names table is interpolated into the SQL builders as a whole `project.dataset.table`
# path, so it gets its own allowlist rather than being decomposed. We generate this name ourselves, so
# this is a backstop against a future caller passing one in.
_TABLE_PATH_RE = re.compile(r'^[A-Za-z0-9][A-Za-z0-9._-]*\.[A-Za-z0-9_]+\.[A-Za-z0-9_]+$')


def _validate_bq_identifiers(project_id, dataset_name):
    """Raise ``ValueError`` if project_id / dataset_name are unsafe to embed in SQL."""
    if not _PROJECT_ID_RE.match(project_id or ''):
        raise ValueError(
            f"invalid project_id '{project_id}': must match {_PROJECT_ID_RE.pattern}")
    if not _DATASET_NAME_RE.match(dataset_name or ''):
        raise ValueError(
            f"invalid dataset_name '{dataset_name}': must match {_DATASET_NAME_RE.pattern}")


def _validate_bq_table_path(table_path):
    """Raise ``ValueError`` if a fully-qualified ``project.dataset.table`` path is unsafe to embed."""
    if not _TABLE_PATH_RE.match(table_path or ''):
        raise ValueError(
            f"invalid table path '{table_path}': must match {_TABLE_PATH_RE.pattern}")


def read_sample_names(path):
    """Read a one-name-per-line sample names file into a de-duplicated list, order preserved.

    Blank lines and surrounding whitespace are ignored so a trailing newline -- which every FOFN
    written by ``GvsBulkIngestGenomes`` has -- does not become an empty sample name. A file with no
    names at all raises ``ValueError``: silently validating nothing would let a broken FOFN pass as
    a clean cohort, which is the opposite of what this fail-fast check is for.
    """
    with open(path) as names_file:
        names = [line.strip() for line in names_file]
    seen = set()
    deduped = []
    for name in names:
        if name and name not in seen:
            seen.add(name)
            deduped.append(name)
    if not deduped:
        raise ValueError(f"sample names file '{path}' contains no sample names")
    return deduped


@contextlib.contextmanager
def stage_sample_names_table(client, project_id, dataset_name, sample_names):
    """Stage ``sample_names`` in a short-lived BigQuery table and yield its full path.

    Yields ``None`` when ``sample_names`` is None, so callers can wrap the unrestricted case in the
    same ``with`` block. The table is created in the dataset under validation (which the ingest that
    just ran already has write access to), loaded from a temporary local CSV, and dropped on the way
    out. A 12-hour expiration is also set at creation time so a process killed mid-run -- a Cromwell
    preemption, say -- cannot leave the table behind permanently.
    """
    if sample_names is None:
        yield None
        return

    from google.cloud import bigquery
    table_id = f"{project_id}.{dataset_name}.header_validation_sample_names_{uuid.uuid4().hex[:12]}"
    _validate_bq_table_path(table_id)
    table = bigquery.Table(table_id, schema=[
        bigquery.SchemaField('sample_name', 'STRING', mode='REQUIRED'),
    ])
    table.expires = datetime.datetime.now(datetime.timezone.utc) + datetime.timedelta(hours=12)
    client.create_table(table)
    try:
        # A load job rather than INSERT DML: the batch can be hundreds of thousands of names, which
        # is far past what is sane to send as statement text.
        with tempfile.TemporaryDirectory() as staging_dir:
            staging_path = os.path.join(staging_dir, 'sample_names.csv')
            with open(staging_path, 'w') as staging_file:
                for name in sample_names:
                    staging_file.write(f"{name}\n")
            job_config = bigquery.LoadJobConfig(
                source_format=bigquery.SourceFormat.CSV,
                schema=table.schema,
                write_disposition=bigquery.WriteDisposition.WRITE_TRUNCATE,
            )
            # Passing the size lets the client use a single multipart upload for anything under its
            # 5 MB threshold instead of a resumable handshake -- one fewer round trip for a typical
            # batch, and the resumable path is what a BigQuery emulator tends not to implement.
            staging_size = os.path.getsize(staging_path)
            with open(staging_path, 'rb') as staging_file:
                client.load_table_from_file(staging_file, table_id, size=staging_size,
                                            job_config=job_config).result()
        print(f"Staged {len(sample_names)} sample names in {table_id}")
        yield table_id
    finally:
        client.delete_table(table_id, not_found_ok=True)


def _make_client(project_id):
    from google.cloud import bigquery
    from google.cloud.bigquery.job import QueryJobConfig
    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE",
                                    use_query_cache=True, use_legacy_sql=False)
    return bigquery.Client(project=project_id, default_query_job_config=default_config)


def run_queries(client, project_id, dataset_name, sample_names_table=None):
    """Run the four check queries, returning ``(summary, dragen_rows, orphan, blob_rows)``.

    Split out from ``run_checks`` so the queries can be exercised against an already-staged names
    table. ``stage_sample_names_table`` needs a BigQuery load job, which the emulator the unit tests
    run against does not implement; the integration test builds the table with plain SQL and calls
    this directly, so the restricted SQL itself is still covered.
    """
    import utils
    summary = list(utils.execute_with_retry(
        client, "per sample header summary",
        per_sample_summary_sql(project_id, dataset_name, sample_names_table))['results'])[0]
    dragen_rows = list(utils.execute_with_retry(
        client, "dragen version breakdown",
        dragen_version_breakdown_sql(project_id, dataset_name, sample_names_table))['results'])
    orphan = list(utils.execute_with_retry(
        client, "orphan header hashes",
        orphan_hash_sql(project_id, dataset_name, sample_names_table))['results'])[0]
    blob_rows = list(utils.execute_with_retry(
        client, "shared blob distribution",
        shared_blob_distribution_sql(project_id, dataset_name, sample_names_table))['results'])
    return summary, dragen_rows, orphan, blob_rows


def run_checks(project_id, dataset_name, expected_dragen_version=None, require_reblocking=True,
               client=None, sample_names=None):
    """Execute all queries and evaluate them. Returns ``(overall_passed, [CheckResult, ...])``.

    ``project_id`` and ``dataset_name`` are interpolated directly into the SQL, so they are first
    validated with ``_validate_bq_identifiers`` (project ids: letters/digits and ``. _ -``; dataset
    names: letters/digits/``_`` only). A value outside that allowlist raises ``ValueError`` before
    any query runs. Example: ``run_checks('broad-dsde-methods', 'aou_wgs_ingest', '3.7.8')``.

    ``sample_names``, when given, is a list of sample names (see ``read_sample_names``) restricting
    every check to those samples; it is staged in a temporary table for the duration of the run.
    """
    _validate_bq_identifiers(project_id, dataset_name)
    if client is None:
        client = _make_client(project_id)

    with stage_sample_names_table(client, project_id, dataset_name, sample_names) as names_table:
        summary, dragen_rows, orphan, blob_rows = run_queries(
            client, project_id, dataset_name, names_table)

    return evaluate(summary, dragen_rows, orphan, blob_rows, expected_dragen_version,
                    require_reblocking,
                    requested_samples=None if sample_names is None else len(sample_names))


def write_output_files(overall_passed, report_text, pass_file_output, report_file_output):
    with open(pass_file_output, 'w') as pass_file:
        pass_file.write('true' if overall_passed else 'false')
    with open(report_file_output, 'w') as report_file:
        report_file.write(report_text)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Validate ingested GVS VCF headers.')
    parser.add_argument('--project_id', type=str, required=True,
                        help='Google project for the GVS dataset')
    parser.add_argument('--dataset_name', type=str, required=True,
                        help='BigQuery dataset name holding the header tables')
    parser.add_argument('--sample_names_file', type=str, required=False, default=None,
                        help='Optional file of sample names, one per line, restricting validation to '
                             'those samples (e.g. the current ingest batch). Omit to validate every '
                             'non-control, non-withdrawn sample in the dataset.')
    parser.add_argument('--expected_dragen_version', type=str, required=False, default=None,
                        help="Expected DRAGEN version. A single triplet (e.g. '3.7.8') requires every "
                             "sample to match it exactly; a range requires every sample's triplet to "
                             "fall within it. Ranges accept interval notation: '3.4.12-3.7.8' (both "
                             "inclusive), '[3.7.8-3.8)' (inclusive-exclusive), '(3.7-3.8)' (both "
                             "exclusive). Bounds compare positionally, so an inclusive two-component "
                             "upper bound is not 'all of that minor': '3.7.8-3.8' excludes 3.8.1. For "
                             "'everything below 3.8' use the exclusive form '3.7.8-3.8)'. If unset, "
                             "only cross-sample consistency is checked.")
    parser.add_argument('--allow_non_reblocked', action='store_true',
                        help='Downgrade the reblocking check from fatal to informational. By default '
                             '(flag absent) a sample without a ReblockGVCF command line fails validation.')
    parser.add_argument('--pass_file_output', type=str, required=False, default='pass.txt',
                        help="Location to write 'true'/'false' pass indicator")
    parser.add_argument('--report_file_output', type=str, required=False, default='report.txt',
                        help='Location to write the human-readable validation report')

    args = parser.parse_args()

    sample_names = read_sample_names(args.sample_names_file) if args.sample_names_file else None

    overall_passed, checks = run_checks(args.project_id, args.dataset_name, args.expected_dragen_version,
                                        require_reblocking=not args.allow_non_reblocked,
                                        sample_names=sample_names)
    report = compose_report(args.project_id, args.dataset_name, args.expected_dragen_version,
                            overall_passed, checks,
                            requested_samples=None if sample_names is None else len(sample_names))
    write_output_files(overall_passed, report, args.pass_file_output, args.report_file_output)
    print(report)
