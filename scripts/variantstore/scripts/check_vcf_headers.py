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
    numeric triplet fails the check in every mode -- such rows are never silently dropped.
  * Shared-blob distribution (informational): how many distinct ``is_expected_unique = FALSE``
    blobs there are and how many samples carry each (> 1 indicates distinct delivery batches).

Whether a failing result aborts the pipeline is decided by the caller
(``GvsValidateVcfHeaders.wdl`` ``fail_on_validation_errors``, default true); this script always
computes an overall pass/fail and writes a human-readable report.
"""
import argparse
import re
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

def per_sample_summary_sql(project_id, dataset_name):
    """One-pass per-sample rollup over the expected (non-control, non-withdrawn) cohort.

    Returns a single summary row: the expected sample count, how many have header data, and the
    counts + example sample names for each fatal integrity / reblocking condition. ``COUNTIF`` over
    a LEFT JOIN treats a sample with no header rows as all-zero counts, so missing-header samples
    surface as ``chunk_count = 0``. Each example array is capped at ``EXAMPLE_LIMIT + 1`` -- one more
    than we display -- so the report can tell "exactly N offenders" from "more than N" (see
    ``_examples``).
    """
    return f"""
        WITH expected AS (
            SELECT sample_id, sample_name
            FROM `{project_id}.{dataset_name}.sample_info`
            WHERE is_control = FALSE AND withdrawn IS NULL
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


def dragen_version_breakdown_sql(project_id, dataset_name):
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
          AND si.is_control = FALSE AND si.withdrawn IS NULL
        GROUP BY sw_version
        ORDER BY n_samples DESC
    """


def orphan_hash_sql(project_id, dataset_name):
    """Referential integrity: associations in ``sample_vcf_header`` whose hash has no
    ``vcf_header_lines`` row. Should always be zero for a correct load."""
    return f"""
        SELECT
            COUNT(*) AS orphan_associations,
            COUNT(DISTINCT svh.sample_id) AS affected_samples
        FROM `{project_id}.{dataset_name}.sample_vcf_header` svh
        LEFT JOIN `{project_id}.{dataset_name}.vcf_header_lines` vhl USING (vcf_header_lines_hash)
        WHERE vhl.vcf_header_lines_hash IS NULL
    """


def shared_blob_distribution_sql(project_id, dataset_name):
    """Distinct shared (``is_expected_unique = FALSE``) blobs and how many samples carry each."""
    return f"""
        SELECT
            svh.vcf_header_lines_hash AS blob_hash,
            COUNT(DISTINCT svh.sample_id) AS n_samples
        FROM `{project_id}.{dataset_name}.sample_vcf_header` svh
        JOIN `{project_id}.{dataset_name}.vcf_header_lines` vhl USING (vcf_header_lines_hash)
        WHERE vhl.is_expected_unique = FALSE
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


def evaluate_integrity(summary):
    """Integrity checks derived from the per-sample summary row. Returns a list of CheckResult."""
    expected = summary.get('expected_samples') or 0
    with_headers = summary.get('samples_with_headers') or 0
    missing_headers = summary.get('samples_missing_headers') or 0
    missing_unique = summary.get('samples_missing_unique_chunk') or 0

    results = []

    # There must be samples to validate at all.
    results.append(CheckResult(
        name='samples_present',
        passed=expected > 0,
        fatal=True,
        lines=[f"Non-control, non-withdrawn samples in sample_info: {expected}"] +
              ([] if expected > 0 else ["    No samples to validate -- did the header ingest run against the right dataset?"]),
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


def evaluate_dragen_version(dragen_rows, expected_dragen_version=None):
    """Consistency (and optionally an exact match or inclusive range) of the DRAGEN version.

    ``dragen_rows``: iterable of rows with ``sw_version`` and ``n_samples``.
    ``expected_dragen_version`` accepts either a single triplet (``'3.7.8'`` -- exact match, as AoU
    requires) or a range with optional interval-notation brackets (``'3.4.12-3.7.8'``,
    ``'[3.7.8-3.8)'``, ``'(3.7-3.8)'``); see ``parse_expected_dragen_spec``. Failures:
      * expected version given but no DRAGEN command lines found;
      * any DRAGEN command line whose SW version cannot be reduced to a numeric triplet (fatal in
        every mode -- these rows are never silently dropped from the comparison);
      * exact / consistency-only mode: more than one distinct version triplet across the cohort;
      * exact mode: any triplet differs from the expected version;
      * range mode: any triplet falls outside the interval (multiple triplets inside the range are
        allowed -- relaxing a range but still requiring a single triplet would be pointless);
      * ``expected_dragen_version`` is a malformed range.
    With no DRAGEN lines and no expected version, this is informational (a non-DRAGEN cohort).
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
    # guards below never leave the cohort silently unchecked. That edge only reopens if the shape of
    # dragen_version_breakdown_sql changes to admit an all-zero/NULL-count row.
    passed = True

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
             require_reblocking=True):
    """Run every check and return ``(overall_passed, [CheckResult, ...])``.

    ``overall_passed`` is true iff every *fatal* check passed; informational checks never affect it.
    """
    checks = []
    checks.extend(evaluate_integrity(summary))
    checks.append(evaluate_reblocking(summary, require_reblocking))
    checks.append(evaluate_dragen_version(dragen_rows, expected_dragen_version))
    checks.append(evaluate_referential_integrity(orphan))
    checks.append(evaluate_shared_blob_distribution(blob_rows))

    overall_passed = all(c.passed for c in checks if c.fatal)
    return overall_passed, checks


def compose_report(project_id, dataset_name, expected_dragen_version, overall_passed, checks):
    """Render the human-readable report (the VS-1215 "small report")."""
    header = [
        "GVS VCF Header Validation Report (VS-1966)",
        "==========================================",
        f"Dataset: {project_id}.{dataset_name}",
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


def _validate_bq_identifiers(project_id, dataset_name):
    """Raise ``ValueError`` if project_id / dataset_name are unsafe to embed in SQL."""
    if not _PROJECT_ID_RE.match(project_id or ''):
        raise ValueError(
            f"invalid project_id '{project_id}': must match {_PROJECT_ID_RE.pattern}")
    if not _DATASET_NAME_RE.match(dataset_name or ''):
        raise ValueError(
            f"invalid dataset_name '{dataset_name}': must match {_DATASET_NAME_RE.pattern}")


def _make_client(project_id):
    from google.cloud import bigquery
    from google.cloud.bigquery.job import QueryJobConfig
    default_config = QueryJobConfig(labels=query_labels_map, priority="INTERACTIVE",
                                    use_query_cache=True, use_legacy_sql=False)
    return bigquery.Client(project=project_id, default_query_job_config=default_config)


def run_checks(project_id, dataset_name, expected_dragen_version=None, require_reblocking=True,
               client=None):
    """Execute all queries and evaluate them. Returns ``(overall_passed, [CheckResult, ...])``.

    ``project_id`` and ``dataset_name`` are interpolated directly into the SQL, so they are first
    validated with ``_validate_bq_identifiers`` (project ids: letters/digits and ``. _ -``; dataset
    names: letters/digits/``_`` only). A value outside that allowlist raises ``ValueError`` before
    any query runs. Example: ``run_checks('broad-dsde-methods', 'aou_wgs_ingest', '3.7.8')``.
    """
    import utils
    _validate_bq_identifiers(project_id, dataset_name)
    if client is None:
        client = _make_client(project_id)

    summary = list(utils.execute_with_retry(
        client, "per sample header summary", per_sample_summary_sql(project_id, dataset_name))['results'])[0]
    dragen_rows = list(utils.execute_with_retry(
        client, "dragen version breakdown", dragen_version_breakdown_sql(project_id, dataset_name))['results'])
    orphan = list(utils.execute_with_retry(
        client, "orphan header hashes", orphan_hash_sql(project_id, dataset_name))['results'])[0]
    blob_rows = list(utils.execute_with_retry(
        client, "shared blob distribution", shared_blob_distribution_sql(project_id, dataset_name))['results'])

    return evaluate(summary, dragen_rows, orphan, blob_rows, expected_dragen_version, require_reblocking)


def write_output_files(overall_passed, report_text, pass_file_output, report_file_output):
    with open(pass_file_output, 'w') as pass_file:
        pass_file.write('true' if overall_passed else 'false')
    with open(report_file_output, 'w') as report_file:
        report_file.write(report_text)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Validate ingested GVS VCF headers (VS-1966).')
    parser.add_argument('--project_id', type=str, required=True,
                        help='Google project for the GVS dataset')
    parser.add_argument('--dataset_name', type=str, required=True,
                        help='BigQuery dataset name holding the header tables')
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

    overall_passed, checks = run_checks(args.project_id, args.dataset_name, args.expected_dragen_version,
                                        require_reblocking=not args.allow_non_reblocked)
    report = compose_report(args.project_id, args.dataset_name, args.expected_dragen_version,
                            overall_passed, checks)
    write_output_files(overall_passed, report, args.pass_file_output, args.report_file_output)
    print(report)
