#!/usr/bin/env python3
"""
Independent, row-count-based structural checks for Parquet ingest verification (VS-1989).

These checks deliberately do NOT reuse ``get_already_loaded_tables_and_sample_ids`` -- the predicate
``parse_and_group_files.py`` uses to decide which files to skip. That predicate answers only one
question, "does this partition have ``total_logical_bytes > 0``?", so a verifier built on the same
call can never contradict the loader: by construction it agrees. It is therefore blind to the two
failure classes that matter most before the source Parquet is deleted:

  * a partial or lost load -- a ``(table, sample_id)`` partition exists (bytes > 0) but holds fewer
    rows than the source, and
  * a duplicated load -- a partition holds *more* rows than the source, which presence-testing can
    never see.

The functions here ask a different question -- "how many rows?" -- through two cheap BigQuery
metadata reads that share no code with the loader predicate:

  * ``INFORMATION_SCHEMA.PARTITIONS.total_rows`` for the superpartitioned ``vet_%`` / ``ref_ranges_%``
    families (partition metadata is a flat ~10 MB regardless of callset size; ``partition_id`` *is*
    the ``sample_id`` because these tables are integer-range partitioned on ``sample_id`` with step
    1), and
  * a per-sample ``COUNT(*)`` on the regular ploidy table (scans a single column, ~100 MB on a
    500k-sample callset).

The crux of the independence is ``total_rows`` vs. ``total_logical_bytes``: a partition that is
present but partial or duplicated has a wrong row count, which these checks catch and the loader
predicate cannot.

Scope notes (VS-1989):
  * The exact footer-vs-BigQuery per-sample comparison (the gold-standard loss+duplication detector
    from the ticket description) is intentionally deferred: it needs a new heavy image dependency
    (pyarrow) and a footer read per file, which is prohibitive at AoU scale.
  * ``ref_ranges`` has no cheap per-sample duplication signal -- its row count tracks GQ-band
    transitions, not genome length, and legitimate samples reach many times the median -- so only
    whole-sample presence is verified for it and the gap is recorded rather than papered over. The
    only detector that does work is exact: ``COUNT(*)`` vs. ``COUNT(DISTINCT packed_ref_data)`` per
    sample, which costs a full per-sample scan rather than a partition-metadata read and so is out
    of scope here. Comparing a sample's row count against its counterpart in the parent callset was
    also considered and rejected: a child callset is seeded by copying its parent, so a defect
    already present in the parent is copied forward identically and reads as agreement, not
    disagreement. Cross-generation comparison can only catch a change introduced at or after the
    copy -- never a defect the parent already had -- so it is not a substitute for a real detector
    and is not implemented. Do not re-attempt either approach without re-reading this note.
  * The ploidy *duplication* reading holds only where a sample's ploidy rows were written at ingest;
    a later backfill writes a full per-sample set regardless of what the reference tables contained,
    masking a doubled sample. See ``PLOIDY_BACKFILL_CAVEAT``.
"""

import logging
import statistics
from collections import Counter, defaultdict

try:
    from google.cloud import bigquery
except ImportError:
    bigquery = None  # type: ignore  # will fail at runtime if BigQuery calls are made without the package installed

# Reuse only the SQL-injection guard from the loader module. This is an input validator, not the
# loader's "what is already loaded?" oracle, so importing it does not compromise the independence
# these checks exist to provide.
from parse_and_group_files import _validate_table_prefixes


log = logging.getLogger(__name__)

# Families for which a per-sample row-count-to-median duplication screen produces a usable signal.
# Only ``vet`` qualifies (tight distribution); ``ref_ranges`` is deliberately excluded (see module
# docstring). Callers may override.
DEFAULT_DUPLICATION_SCREEN_FAMILIES = ["vet"]

# Default ratio-to-median above which a vet sample is flagged as a possible duplicate. Miguel's
# Foxtrot calibration found vet tight enough that 1.2x is safe against false positives; 1.6x is a
# conservative default that still catches a doubled sample (~2.0x).
DEFAULT_VET_DUPLICATION_THRESHOLD = 1.6

PLOIDY_BACKFILL_CAVEAT = (
    "Ploidy duplication is only detectable here where a sample's ploidy rows were written at ingest. "
    "A later backfill writes a full per-sample set regardless, so a doubled sample would still show "
    "the modal row count. Do not read a green ploidy result as a duplication guarantee for callsets "
    "or sample ranges whose ploidy was backfilled (provenance is not captured today -- see VS-1989)."
)


def family_for_table(table_name, superpartitioned_table_prefixes):
    """
    Return the superpartitioned family prefix a concrete table name belongs to, or None.

    A superpartitioned table name is ``<prefix>_<digits>`` (e.g. ``vet_001`` -> ``vet``,
    ``ref_ranges_042`` -> ``ref_ranges``). The full-string match on ``^<prefix>_[0-9]+$`` avoids a
    shorter prefix swallowing a longer table name.
    """
    for prefix in superpartitioned_table_prefixes:
        remainder = table_name[len(prefix) + 1:]
        if table_name.startswith(prefix + "_") and remainder.isdigit():
            return prefix
    return None


def get_partition_row_counts(project_id, dataset_name, superpartitioned_table_prefixes):
    """
    Return per-partition row counts for the superpartitioned families as a list of
    ``(table_name, sample_id, total_rows)`` tuples.

    Reads ``INFORMATION_SCHEMA.PARTITIONS.total_rows`` directly. Unlike the loader predicate this
    does NOT filter on ``total_logical_bytes > 0`` -- a present-but-empty partition (row count 0) is
    exactly the partial-load signal we want to surface, so it must be included.
    """
    _validate_table_prefixes(superpartitioned_table_prefixes, [])

    if not superpartitioned_table_prefixes:
        return []

    superpartitioned_regex = "|".join(
        f"^{prefix}_[0-9]+$" for prefix in superpartitioned_table_prefixes
    )
    query = f"""
        SELECT table_name AS table_name,
               SAFE_CAST(partition_id AS INT64) AS sample_id,
               total_rows AS total_rows
        FROM `{project_id}.{dataset_name}.INFORMATION_SCHEMA.PARTITIONS`
        WHERE
            REGEXP_CONTAINS(table_name, '{superpartitioned_regex}') AND
            NOT STARTS_WITH(partition_id, '__') AND
            SAFE_CAST(partition_id AS INT64) IS NOT NULL
        ORDER BY table_name, sample_id
    """

    try:
        client = bigquery.Client(project=project_id)
        results = client.query(query)
        rows = [(row.table_name, row.sample_id, row.total_rows) for row in results]
        log.info(f"Read {len(rows)} partition row counts from INFORMATION_SCHEMA.PARTITIONS")
        return rows
    except Exception as e:
        log.error(f"ERROR: Could not read partition row counts from BigQuery: {e}")
        raise


def get_ploidy_row_counts(project_id, dataset_name, ploidy_table):
    """
    Return ``{sample_id: row_count}`` for a regular (non-superpartitioned) per-sample table by
    grouping on ``sample_id``. Samples with zero rows simply do not appear.
    """
    _validate_table_prefixes([], [ploidy_table])

    query = f"""
        SELECT sample_id AS sample_id, COUNT(*) AS n
        FROM `{project_id}.{dataset_name}.{ploidy_table}`
        GROUP BY sample_id
    """

    try:
        client = bigquery.Client(project=project_id)
        results = client.query(query)
        counts = {row.sample_id: row.n for row in results}
        log.info(f"Read per-sample row counts for {len(counts)} samples from {ploidy_table}")
        return counts
    except Exception as e:
        log.error(f"ERROR: Could not read per-sample row counts for {ploidy_table}: {e}")
        raise


def assess_family_completeness(partition_rows, regular_counts, expected_by_family,
                               superpartitioned_table_prefixes, regular_table_prefixes):
    """
    Check that every sample GCS says should exist is present with a non-empty row count in each of
    its families (the superpartitioned ``vet_%`` / ``ref_ranges_%`` and the regular ploidy table).

    A sample expected in a superpartitioned family whose partition has ``total_rows == 0`` is
    reported separately as ``empty`` -- present to the loader predicate (it may have bytes) but empty
    in fact, i.e. a partial load the loader-shared check cannot see.

    Returns a dict with an overall ``ok`` and a per-family breakdown.
    """
    present_by_family = defaultdict(set)
    empty_by_family = defaultdict(set)

    for table_name, sample_id, total_rows in partition_rows:
        fam = family_for_table(table_name, superpartitioned_table_prefixes)
        if fam is None:
            continue
        if total_rows and total_rows > 0:
            present_by_family[fam].add(sample_id)
        else:
            empty_by_family[fam].add(sample_id)

    for prefix in regular_table_prefixes:
        counts = regular_counts.get(prefix, {})
        present_by_family[prefix] = {sid for sid, n in counts.items() if n and n > 0}

    per_family = {}
    overall_ok = True
    for fam, expected in sorted(expected_by_family.items()):
        present = present_by_family.get(fam, set())
        # An "empty" partition only matters if that sample was expected for this run.
        empty = sorted(empty_by_family.get(fam, set()) & set(expected))
        missing = sorted(set(expected) - present)
        fam_ok = not missing and not empty
        overall_ok = overall_ok and fam_ok
        per_family[fam] = {
            "ok": fam_ok,
            "expected": len(expected),
            "present": len(present & set(expected)),
            "missing_samples": missing,
            "empty_partition_samples": empty,
        }

    return {"ok": overall_ok, "per_family": per_family}


def assess_cardinality(counts, expected_samples, expected_count=None):
    """
    Check that every expected sample has the *same* per-sample row count and that none is missing.

    By default the reference is the callset's own modal count, deliberately *not* a fixed number: for
    ploidy the intuitive "24" is emergent (the count of distinct non-PAR contigs, which chrM or
    non-WGS ingest can change), so hardcoding it would false-positive on exome/BGE/chrM callsets. A
    sample below the reference indicates a partial load; a sample at a multiple of it indicates
    duplication (subject to PLOIDY_BACKFILL_CAVEAT).

    When ``expected_count`` is supplied the reference becomes that exact constant instead of the mode
    (e.g. pass 24 for a WGS ploidy table). This recovers the ticket's exact-cardinality guarantee and
    closes the one gap mode-based detection has -- if a majority of samples share the *wrong* count,
    that wrong count becomes the mode and the correct minority would otherwise be flagged. The
    observed mode is always reported alongside, so a mismatch between an override and reality is
    visible rather than silently overriding it.

    ``counts`` is ``{sample_id: row_count}``; only ``expected_samples`` are assessed (extra samples
    already in the table from prior loads are ignored).
    """
    expected = set(expected_samples)
    present = {sid: counts[sid] for sid in expected if sid in counts and counts[sid] > 0}
    missing = sorted(expected - set(present))

    reference_source = "override" if expected_count is not None else "mode"

    if not present:
        return {
            "ok": not expected,
            "mode": None,
            "reference_count": expected_count,
            "reference_source": reference_source,
            "min": None,
            "max": None,
            "distinct_samples": 0,
            "total_rows": 0,
            "missing_samples": missing,
            "deviating_samples": [],
        }

    values = list(present.values())
    observed_mode = Counter(values).most_common(1)[0][0]
    reference = expected_count if expected_count is not None else observed_mode
    deviating = sorted(
        ({"sample_id": sid, "count": n} for sid, n in present.items() if n != reference),
        key=lambda d: d["sample_id"],
    )

    return {
        "ok": not missing and not deviating,
        "mode": observed_mode,
        "reference_count": reference,
        "reference_source": reference_source,
        "min": min(values),
        "max": max(values),
        "distinct_samples": len(present),
        "total_rows": sum(values),
        "missing_samples": missing,
        "deviating_samples": deviating,
    }


def assess_duplication_screen(partition_rows, family, expected_samples,
                              superpartitioned_table_prefixes, threshold):
    """
    Flag samples in ``family`` whose ``total_rows`` is at least ``threshold`` times the callset
    median -- a heuristic screen for duplication. Returns the flagged samples; the caller decides
    whether flags gate (strict) or merely warn (default).

    Only meaningful for families with a tight per-sample distribution (``vet``); callers must not
    apply it to ``ref_ranges``.
    """
    expected = set(expected_samples)
    rows_by_sample = {
        sample_id: total_rows
        for table_name, sample_id, total_rows in partition_rows
        if family_for_table(table_name, superpartitioned_table_prefixes) == family
        and sample_id in expected
        and total_rows and total_rows > 0
    }

    if not rows_by_sample:
        return {"family": family, "threshold": threshold, "median": None,
                "samples_screened": 0, "outliers": []}

    median = statistics.median(rows_by_sample.values())
    outliers = []
    if median > 0:
        for sid, rows in rows_by_sample.items():
            ratio = rows / median
            if ratio >= threshold:
                outliers.append({"sample_id": sid, "rows": rows, "ratio": round(ratio, 3)})
    outliers.sort(key=lambda d: d["ratio"], reverse=True)

    return {
        "family": family,
        "threshold": threshold,
        "median": median,
        "samples_screened": len(rows_by_sample),
        "outliers": outliers,
    }


def run_structural_checks(project_id, dataset_name, expected_by_family,
                          superpartitioned_table_prefixes=None, regular_table_prefixes=None,
                          vet_duplication_threshold=DEFAULT_VET_DUPLICATION_THRESHOLD,
                          strict_vet_screen=False,
                          expected_ploidy_rows_per_sample=None,
                          duplication_screen_families=None):
    """
    Run the independent structural checks and return a result dict.

    ``expected_by_family`` maps a family name (a superpartitioned prefix such as ``vet`` /
    ``ref_ranges``, or a regular table name such as ``sample_chromosome_ploidy``) to the set of
    ``sample_id`` values GCS says should be loaded for it.

    ``expected_ploidy_rows_per_sample``, when set, is the exact per-sample row count each regular
    per-sample table (ploidy) is validated against instead of the callset mode; leave it unset to
    infer the reference from the data. (There is one regular table -- ploidy -- today; if others with
    differing exact counts are ever added this would need to become a per-table mapping.)

    The returned dict carries flat booleans the WDL reads shallowly (``completeness_ok``,
    ``cardinality_ok``, ``duplication_flagged``) plus a nested ``details`` block for humans and logs.
    ``completeness_ok`` and ``cardinality_ok`` are the hard-gate signals; the duplication screen only
    contributes to the gate when ``strict_vet_screen`` is set.
    """
    if superpartitioned_table_prefixes is None:
        superpartitioned_table_prefixes = ["vet", "ref_ranges"]
    if regular_table_prefixes is None:
        regular_table_prefixes = ["sample_chromosome_ploidy"]
    if duplication_screen_families is None:
        duplication_screen_families = DEFAULT_DUPLICATION_SCREEN_FAMILIES

    partition_rows = get_partition_row_counts(
        project_id, dataset_name, superpartitioned_table_prefixes
    )
    regular_counts = {
        prefix: get_ploidy_row_counts(project_id, dataset_name, prefix)
        for prefix in regular_table_prefixes
    }

    completeness = assess_family_completeness(
        partition_rows, regular_counts, expected_by_family,
        superpartitioned_table_prefixes, regular_table_prefixes,
    )

    # Per-sample cardinality consistency, applied to each regular per-sample table (ploidy is the
    # exemplar): every sample should carry the callset's modal row count.
    cardinality = {}
    cardinality_ok = True
    for prefix in regular_table_prefixes:
        result = assess_cardinality(
            regular_counts.get(prefix, {}),
            expected_by_family.get(prefix, set()),
            expected_count=expected_ploidy_rows_per_sample,
        )
        result["backfill_caveat"] = PLOIDY_BACKFILL_CAVEAT
        cardinality[prefix] = result
        cardinality_ok = cardinality_ok and result["ok"]

    # Duplication screen on the usable superpartitioned families only (vet). ref_ranges is recorded
    # as intentionally unchecked.
    duplication = {}
    for family in duplication_screen_families:
        if family in superpartitioned_table_prefixes:
            duplication[family] = assess_duplication_screen(
                partition_rows, family, expected_by_family.get(family, set()),
                superpartitioned_table_prefixes, vet_duplication_threshold,
            )
    unscreened = {
        "checked": False,
        "reason": ("No cheap per-sample duplication detector exists for these families (row count "
                   "tracks GQ-band transitions, not genome length). The exact detector -- COUNT(*) "
                   "vs. COUNT(DISTINCT packed_ref_data) per sample -- needs a full per-sample scan "
                   "and is out of scope. Comparing against the parent callset's row count was also "
                   "considered and rejected: a child callset copies its parent, so a pre-existing "
                   "parent-side defect is copied forward identically and would read as agreement. "
                   "See the module docstring and VS-1989."),
        "families": sorted(f for f in superpartitioned_table_prefixes if f not in duplication),
    }

    duplication_flagged = any(d["outliers"] for d in duplication.values())

    return {
        "completeness_ok": completeness["ok"],
        "cardinality_ok": cardinality_ok,
        "duplication_flagged": duplication_flagged,
        "strict_vet_screen": strict_vet_screen,
        "details": {
            "family_completeness": completeness,
            "cardinality": cardinality,
            "duplication_screen": duplication,
            "duplication_unscreened": unscreened,
        },
    }
