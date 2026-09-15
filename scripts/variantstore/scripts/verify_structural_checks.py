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
    (pyarrow) and a footer read per file, which is prohibitive at AoU scale. It would also catch
    little that the checks here do not. A partial *load* cannot occur: BigQuery load jobs are atomic
    (``load_parquet_to_bq.py`` loads each batch via ``load_table_from_uri`` with ``WRITE_APPEND``, and
    a failed job commits nothing), so a present partition is either complete, empty (caught by
    completeness), or duplicated (caught by the duplication screen). The one residual truncation source
    is a Parquet file generated upstream with too few rows -- and there the footer count and the
    BigQuery count agree, so footer-vs-BigQuery would pass it too. That case is instead surfaced
    cheaply by ``assess_truncation_screen``, a below-median heuristic (the low-side mirror of the
    duplication screen); catching it exactly would need a gVCF-level variant count, which is out of
    scope. Do not re-attempt the footer comparison without re-reading this note.
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
import math
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

# Regular (non-superpartitioned) tables held to cardinality consistency by the cardinality check.
# When an explicit expected count is configured (e.g. expected_ploidy_rows_per_sample=24 for WGS),
# all expected samples must match that exact count. When unset, the callset mode is inferred:
# legitimate contig variations (such as 23 or 25 when mode is 24, since RefRangesCreator only records
# ploidy for contigs with usable reference blocks) pass, while duplicated loads (>= 1.5x baseline, such
# as 48 vs 24) are flagged as deviating and fail cardinality. Only ploidy qualifies for cardinality
# checking: tables whose per-sample count legitimately varies across samples (such as
# vcf_header_lines_scratch) are excluded entirely.
DEFAULT_CARDINALITY_TABLE_PREFIXES = ["sample_chromosome_ploidy"]

# Families that Parquet generation produces TOGETHER, per sample, so their sample sets must be
# identical within a run (the premise the cross-family consistency check rests on): a sample's vet,
# ref_ranges and ploidy files are all written in the same pass (CreateVariantIngestFiles constructs all
# three creators unconditionally per sample). ``vcf_header_lines_scratch`` is deliberately NOT a member
# -- it is produced only on the header-loading path, on its own per-sample cadence, and a supported
# headers-only ingest produces it while producing none of the data families. The cross-family check
# activates only when at least one member is present, so a headers-only run (no data family present) is
# not required to have them and passes; see assess_cross_family_consistency.
DEFAULT_CO_PRODUCED_FAMILIES = ["vet", "ref_ranges", "sample_chromosome_ploidy"]

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


def _extract_count_info(val):
    """
    Extract (row_count, distinct_chromosomes) from a count value that may be:
    - a raw int: e.g. 24 -> (24, None)
    - a tuple: e.g. (24, 24) -> (24, 24)
    - a dict: e.g. {"count": 24, "distinct": 24} -> (24, 24)
    """
    if isinstance(val, tuple):
        return val[0], (val[1] if len(val) > 1 else None)
    if isinstance(val, dict):
        return val.get("count", val.get("total", 0)), val.get("distinct")
    return val, None


def get_ploidy_row_counts(project_id, dataset_name, ploidy_table):
    """
    Return per-sample row counts for a regular (non-superpartitioned) table. For tables containing
    a ``chromosome`` column (such as ``sample_chromosome_ploidy``), returns
    ``{sample_id: (row_count, distinct_chromosome_count)}`` so duplicate row detection can compare
    COUNT(*) against COUNT(DISTINCT chromosome). For other regular tables, returns ``{sample_id: row_count}``.
    Samples with zero rows do not appear.
    """
    _validate_table_prefixes([], [ploidy_table])

    if ploidy_table == "sample_chromosome_ploidy":
        query = f"""
            SELECT sample_id AS sample_id, COUNT(*) AS n, COUNT(DISTINCT chromosome) AS n_distinct
            FROM `{project_id}.{dataset_name}.{ploidy_table}`
            GROUP BY sample_id
        """
    else:
        query = f"""
            SELECT sample_id AS sample_id, COUNT(*) AS n
            FROM `{project_id}.{dataset_name}.{ploidy_table}`
            GROUP BY sample_id
        """

    try:
        client = bigquery.Client(project=project_id)
        results = client.query(query)
        if ploidy_table == "sample_chromosome_ploidy":
            counts = {row.sample_id: (row.n, getattr(row, "n_distinct", None)) for row in results}
        else:
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
    in fact, i.e. a partial load the loader-shared check cannot see. It appears under
    ``empty_partition_samples`` only and is never also counted under ``missing_samples`` -- "missing"
    means the sample has no partition at all, so the two lists are disjoint.

    Returns a dict with an overall ``ok`` and a per-family breakdown.
    """
    present_by_family = defaultdict(set)
    empty_by_family = defaultdict(set)

    for table_name, sample_id, total_rows in partition_rows:
        fam = family_for_table(table_name, superpartitioned_table_prefixes)
        if fam is None:
            continue
        # Nonzero partition row count check. Per VS-1989 design (see module docstring Scope notes),
        # exact per-file Parquet footer reads are intentionally deferred due to pyarrow dependency and
        # scale limits at 400k+ samples. BigQuery load jobs commit atomically; partial/truncated data is
        # screened downstream via assess_truncation_screen on vet and assess_cardinality on ploidy.
        if total_rows and total_rows > 0:
            present_by_family[fam].add(sample_id)
        else:
            empty_by_family[fam].add(sample_id)

    for prefix in regular_table_prefixes:
        counts = regular_counts.get(prefix, {})
        present_by_family[prefix] = {
            sid for sid, val in counts.items()
            if _extract_count_info(val)[0] and _extract_count_info(val)[0] > 0
        }

    per_family = {}
    overall_ok = True
    for fam, expected in sorted(expected_by_family.items()):
        expected_set = set(expected)
        present = present_by_family.get(fam, set())
        # A present-but-empty partition (total_rows == 0) is a partial load, reported on its own; an
        # empty only matters if that sample was expected for this run.
        empty_set = empty_by_family.get(fam, set()) & expected_set
        empty = sorted(empty_set)
        # "Missing" means no partition at all -- neither present nor empty. Excluding the empties here
        # keeps the two lists disjoint so a present-but-empty sample is not double-counted as missing.
        missing = sorted(expected_set - present - empty_set)
        fam_ok = not missing and not empty
        overall_ok = overall_ok and fam_ok
        per_family[fam] = {
            "ok": fam_ok,
            "expected": len(expected_set),
            "present": len(present & expected_set),
            "missing_samples": missing,
            "empty_partition_samples": empty,
        }

    return {"ok": overall_ok, "per_family": per_family}


def assess_cross_family_consistency(expected_by_family, co_produced_families):
    """
    Cross-check the co-produced families' expected-sample sets against each other (VS-1989; narrows VS-2016).

    assess_family_completeness derives each family's expected set from that family's OWN GCS output
    files, so a sample for which one family's file was never produced is simply not "expected" in that
    family and its absence there passes vacuously -- meanwhile its files in the OTHER families load and
    verify cleanly, and deletion of a source that is in fact incomplete gets authorized. On the Parquet
    ingest path every sample produces a vet, a ref_ranges and a ploidy file together
    (CreateVariantIngestFiles constructs all three creators unconditionally per sample, and a zero-row
    file is still written on close), so within a single run these families' sample sets must be
    identical. Any sample present in some but absent from another is therefore a partial upload -- a real
    incompleteness -- and is caught here from the GCS listing alone, with no external sample source.

    ``co_produced_families`` is the group known to be generated together per sample (vet / ref_ranges /
    ploidy; see DEFAULT_CO_PRODUCED_FAMILIES). The check is **presence-activated**: it fires only when at
    least one member actually has files this run. This is what makes it correct across the phased ingest
    without a phase flag. A member present-but-absent-from-another is a per-sample partial upload; a
    member entirely absent while a sibling is present is a whole-family partial upload -- both fail. But
    a run that produced NONE of the group (a supported headers-only ingest, whose only family is
    ``vcf_header_lines_scratch`` -- deliberately not a member of this group -- while the data prefixes are
    still configured, so their files are simply absent) has no member present, so the check is dormant
    and passes. Membership must be limited to families with a single shared per-sample cadence; a family
    on a different production cadence (headers) must not be listed, or a legitimate phase that produces
    one but not the other would false-positive.

    This still does NOT close the whole-sample gap tracked in VS-2016: a sample absent from EVERY family
    was never in the GCS listing at all, so it is in no family's set and no union, and stays invisible
    here; catching it needs a non-GCS expected-sample source (this run's ingest FOFN). Unlike a whole
    family, whose identity we know from the co-produced group, a never-produced sample's identity is
    unknowable without that source.

    ``expected_by_family`` is ``{family: set(sample_id)}``. Returns an overall ``ok`` plus, per member
    family, the sample_ids in the cross-family union that this family is missing. When no member is
    present the returned ``per_family`` is empty and ``ok`` is True.
    """
    families = sorted(co_produced_families)
    # Presence-activation: unless at least one co-produced family carries files this run, there is
    # nothing to cross-check -- treating an absent group as a gap would fail a legitimate headers-only
    # ingest (which produces none of these families) even though the data prefixes are configured.
    if not any(expected_by_family.get(fam) for fam in families):
        return {"ok": True, "union_size": 0, "per_family": {}}

    union = set()
    for fam in families:
        union |= expected_by_family.get(fam, set())

    per_family = {}
    overall_ok = True
    for fam in families:
        present = expected_by_family.get(fam, set())
        missing = sorted(union - present)
        fam_ok = not missing
        overall_ok = overall_ok and fam_ok
        per_family[fam] = {"ok": fam_ok, "missing_samples": missing}

    return {"ok": overall_ok, "union_size": len(union), "per_family": per_family}


def assess_cardinality(counts, expected_samples, expected_count=None):
    """
    Check that every expected sample has the expected per-sample row count and that none is missing
    or duplicated.

    When ``expected_count`` is supplied, cardinality is strictly enforced against that exact constant
    (e.g. pass 24 for a WGS ploidy table). Any sample whose row count differs from ``expected_count``
    is flagged as deviating, failing ``ok``.

    When ``expected_count`` is None, duplication and partial loads are detected via:
    1. Exact chromosome collision: SamplePloidyCreator writes at most one row per chromosome, so
       comparing COUNT(*) against COUNT(DISTINCT chromosome) identifies duplicated loads exactly.
    2. Modal ratio: row count >= 1.5x baseline flags duplicated loads (such as 48 vs 24).
    3. Modal floor: row count < mode - 2 (for mode >= 10) flags partial loads missing autosomes (such
       as 20 vs 24) while allowing legitimate karyotype variations (such as 23 vs 24 or 25 with chrM).

    ``counts`` can map sample_id to raw row counts (int), tuples ``(total_rows, distinct_chromosomes)``,
    or dicts ``{"count": int, "distinct": int}``. Only ``expected_samples`` are assessed.
    """
    expected = set(expected_samples)
    present = {
        sid: counts[sid] for sid in expected
        if sid in counts and (_extract_count_info(counts[sid])[0] or 0) > 0
    }
    missing = sorted(expected - set(present))

    reference_source = "override" if expected_count is not None else "mode"

    if not present:
        return {
            "ok": not expected,
            "mode": None,
            "reference_count": expected_count,
            "reference_source": reference_source,
            "baseline": None,
            "min": None,
            "max": None,
            "distinct_samples": 0,
            "total_rows": 0,
            "missing_samples": missing,
            "deviating_samples": [],
        }

    values = [_extract_count_info(val)[0] for val in present.values()]
    freq = Counter(values)
    max_freq = max(freq.values())
    top_modes = sorted([val for val, count in freq.items() if count == max_freq])

    if len(top_modes) == 1:
        observed_mode = top_modes[0]
    else:
        # Non-unique modes (e.g. {20, 20, 24, 24} or {23, 23, 24, 24}):
        # Break ties deterministically independent of dictionary iteration order. If tied modes
        # include a gross duplicate (ratio >= 1.5), choose the lower count as the un-duplicated
        # baseline (e.g. 24 over 48). Otherwise, choose the higher count as the conservative
        # complete-genome baseline (e.g. 24 over 20), ensuring truncated samples missing autosomes
        # (< mode - 2) are reliably flagged.
        min_top = top_modes[0]
        max_top = top_modes[-1]
        if min_top > 0 and max_top / min_top >= 1.5:
            observed_mode = min_top
        else:
            observed_mode = max_top

    if expected_count is not None:
        reference = expected_count
        baseline = expected_count
        deviating = []
        for sid, val in sorted(present.items(), key=lambda d: d[0]):
            n, n_distinct = _extract_count_info(val)
            if n != expected_count or (n_distinct is not None and n > n_distinct):
                deviating.append({"sample_id": sid, "count": n})
    else:
        reference = observed_mode
        deviating = []
        if len(values) == 2:
            lower_val = min(values)
            upper_val = max(values)
            baseline = lower_val
            if lower_val > 0 and upper_val / lower_val >= 1.5:
                # Duplication for N=2: flag the high sample (and any chromosome collision)
                for sid, val in sorted(present.items(), key=lambda d: d[0]):
                    n, n_distinct = _extract_count_info(val)
                    if (n_distinct is not None and n > n_distinct) or n == upper_val:
                        deviating.append({"sample_id": sid, "count": n})
            elif upper_val >= 10 and lower_val < upper_val - 2:
                # Truncation for N=2: flag the low sample (and any chromosome collision)
                for sid, val in sorted(present.items(), key=lambda d: d[0]):
                    n, n_distinct = _extract_count_info(val)
                    if (n_distinct is not None and n > n_distinct) or n == lower_val:
                        deviating.append({"sample_id": sid, "count": n})
            else:
                # Permitted minor variation (e.g. 23 vs 24): flag only chromosome collision
                for sid, val in sorted(present.items(), key=lambda d: d[0]):
                    n, n_distinct = _extract_count_info(val)
                    if n_distinct is not None and n > n_distinct:
                        deviating.append({"sample_id": sid, "count": n})
        else:
            baseline = observed_mode
            # Legitimate karyotype variations in human callsets include female samples lacking chrY
            # (mode - 1), or samples lacking chrM when mode includes it (mode - 2). A count below
            # mode - 2 indicates missing autosomes (a partial ploidy load), which fails cardinality.
            min_allowed = max(1, observed_mode - 2) if observed_mode >= 10 else max(1, observed_mode - 1)
            for sid, val in sorted(present.items(), key=lambda d: d[0]):
                n, n_distinct = _extract_count_info(val)
                if ((n_distinct is not None and n > n_distinct)
                        or (baseline > 0 and n / baseline >= 1.5)
                        or (n < min_allowed)):
                    deviating.append({"sample_id": sid, "count": n})

    return {
        "ok": not missing and not deviating,
        "mode": observed_mode,
        "reference_count": reference,
        "reference_source": reference_source,
        "baseline": baseline,
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
    whether a flag blocks Parquet deletion (the default) or is waived (``allow_flagged_vet_loads``).

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
        return {"family": family, "threshold": threshold, "median": None, "baseline": None,
                "samples_screened": 0, "outliers": [], "singleton_flagged": False}

    median = statistics.median(rows_by_sample.values())
    # In a two-sample load, median averages the two values, diluting a 2x duplicate to 1.333x
    # (which misses the default 1.6x threshold). Use the lower value as the duplication baseline.
    baseline = min(rows_by_sample.values()) if len(rows_by_sample) == 2 else median

    outliers = []
    if baseline > 0:
        for sid, rows in rows_by_sample.items():
            ratio = rows / baseline
            if ratio >= threshold:
                outliers.append({"sample_id": sid, "rows": rows, "ratio": round(ratio, 3)})
    outliers.sort(key=lambda d: d["ratio"], reverse=True)

    singleton_flagged = len(rows_by_sample) == 1

    return {
        "family": family,
        "threshold": threshold,
        "median": median,
        "baseline": baseline,
        "samples_screened": len(rows_by_sample),
        "outliers": outliers,
        "singleton_flagged": singleton_flagged,
    }


def assess_truncation_screen(partition_rows, family, expected_samples,
                             superpartitioned_table_prefixes, threshold):
    """
    Flag samples in ``family`` whose ``total_rows`` is at most ``median / threshold`` -- the low-side
    mirror of ``assess_duplication_screen`` and a heuristic screen for a grossly truncated partition.
    The same ``threshold`` governs both sides: a sample reads as a possible duplicate above
    `threshold * median`` and as a possible truncation below ``median / threshold``. Returns the
    flagged samples; the caller decides whether a flag blocks Parquet deletion (the default) or is
    waived (``allow_flagged_vet_loads``).

    This is the only cheap detector for the one truncation source the rest of the row-count gate
    misses: a Parquet file generated upstream with far fewer rows than its peers. Truncation cannot
    arise from the load itself -- BigQuery load jobs are atomic (``load_parquet_to_bq.py`` loads each
    batch via ``load_table_from_uri`` with ``WRITE_APPEND``, and a failed job commits nothing), so a
    present partition is either complete, empty (zero rows, caught by completeness), or duplicated
    (caught by the duplication screen), never partially committed. A load-vs-source row comparison
    (the deferred footer-vs-BigQuery check) would therefore add cost -- a footer read per file, a new
    pyarrow dependency -- without catching anything this does not: where a source Parquet is itself
    truncated, its footer count and the BigQuery count agree, so that comparison would pass it too.
    Catching the upstream case exactly would need a gVCF-level count, which is out of scope; a
    below-peers heuristic is what remains, and that is this.

    Only meaningful for families with a tight per-sample distribution (``vet``); callers must not apply
    it to ``ref_ranges``, whose row count varies by orders of magnitude between legitimate samples.
    Zero-row partitions are excluded here -- that is the completeness check's empty-partition case --
    so this screen judges only partitions that are present and non-empty.
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
        return {"family": family, "threshold": threshold, "median": None, "baseline": None,
                "samples_screened": 0, "outliers": [], "singleton_flagged": False}

    median = statistics.median(rows_by_sample.values())
    # In a two-sample load, median averages the two values, diluting a 0.5x truncation to 0.667x
    # (which misses the default 1/1.6 = 0.625 floor). Use the upper value as the truncation baseline.
    baseline = max(rows_by_sample.values()) if len(rows_by_sample) == 2 else median

    outliers = []
    if baseline > 0:
        floor = baseline / threshold
        for sid, rows in rows_by_sample.items():
            if rows <= floor:
                outliers.append({"sample_id": sid, "rows": rows, "ratio": round(rows / baseline, 3)})
    outliers.sort(key=lambda d: d["ratio"])

    singleton_flagged = len(rows_by_sample) == 1

    return {
        "family": family,
        "threshold": threshold,
        "median": median,
        "baseline": baseline,
        "samples_screened": len(rows_by_sample),
        "outliers": outliers,
        "singleton_flagged": singleton_flagged,
    }


def run_structural_checks(project_id, dataset_name, expected_by_family,
                          superpartitioned_table_prefixes=None, regular_table_prefixes=None,
                          vet_duplication_threshold=DEFAULT_VET_DUPLICATION_THRESHOLD,
                          allow_flagged_vet_loads=False,
                          expected_ploidy_rows_per_sample=None,
                          duplication_screen_families=None,
                          cardinality_table_prefixes=None,
                          co_produced_families=None):
    """
    Run the independent structural checks and return a result dict.

    ``expected_by_family`` maps a family name (a superpartitioned prefix such as ``vet`` /
    ``ref_ranges``, or a regular table name such as ``sample_chromosome_ploidy``) to the set of
    ``sample_id`` values GCS says should be loaded for it.

    ``regular_table_prefixes`` are the non-superpartitioned per-sample tables checked for completeness
    (every expected sample present with > 0 rows). ``cardinality_table_prefixes`` is the subset of
    those additionally held to a UNIFORM per-sample row count; it defaults to ploidy alone. A regular
    table whose per-sample row count legitimately varies -- ``vcf_header_lines_scratch`` writes one row
    per header line and samples differ -- must NOT be listed here, or every sample would read as
    "deviating" and fail a valid ingest; completeness alone covers it.

    ``expected_ploidy_rows_per_sample``, when set, is the exact per-sample row count each
    cardinality-checked table is validated against (e.g. 24 for WGS). When unset, the callset's
    modal count is used as the reference baseline: legitimate contig variations (e.g. 23 or 25 when
    mode is 24, as RefRangesCreator only records ploidy for contigs with usable reference blocks) pass,
    while duplicated loads (>= 1.5x baseline, such as 48 vs 24) are flagged as deviating and fail
    cardinality. (Today ploidy is the only cardinality-checked table; if others with differing exact
    counts are ever added this single override would need to become a per-table mapping.)

    The returned dict carries flat booleans read shallowly downstream (``completeness_ok``,
    ``cardinality_ok``, ``cross_family_ok``, ``duplication_flagged``, ``truncation_flagged``) plus a
    nested ``details`` block for humans and logs. ``completeness_ok``, ``cardinality_ok`` and
    ``cross_family_ok`` are the exact signals that gate ``all_loaded`` (and so the fail-loud abort).
    The duplication and truncation screens never
    affect ``all_loaded``; they gate only the separate ``safe_to_delete_parquet`` predicate, and there
    only when ``allow_flagged_vet_loads`` is false (its default). ``allow_flagged_vet_loads`` itself is
    carried through unchanged, purely so the log summary can say whether a flag will block deletion.
    """
    # Reject a nonsensical ratio before any BigQuery read. 0 divides by zero in the truncation screen
    # (median / threshold); any value <= 1 makes ordinary samples satisfy an outlier condition on both
    # sides (>= threshold * median above, <= median / threshold below), so the screens would flag half
    # the callset. Fail fast with an actionable message rather than crash mid-verification or silently
    # flag everything.
    try:
        threshold_ok = math.isfinite(vet_duplication_threshold) and vet_duplication_threshold > 1
    except TypeError:
        threshold_ok = False
    if not threshold_ok:
        raise ValueError(
            f"vet_duplication_threshold must be a finite number > 1 (got {vet_duplication_threshold!r}); "
            "a ratio <= 1 flags ordinary samples on both screens and 0 divides by zero in the "
            "truncation screen."
        )

    if superpartitioned_table_prefixes is None:
        superpartitioned_table_prefixes = ["vet", "ref_ranges"]
    if regular_table_prefixes is None:
        regular_table_prefixes = ["sample_chromosome_ploidy"]
    if duplication_screen_families is None:
        duplication_screen_families = DEFAULT_DUPLICATION_SCREEN_FAMILIES
    if cardinality_table_prefixes is None:
        cardinality_table_prefixes = DEFAULT_CARDINALITY_TABLE_PREFIXES
    if co_produced_families is None:
        co_produced_families = DEFAULT_CO_PRODUCED_FAMILIES

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

    # Cross-check the co-produced families' sample sets against one another. Completeness above judges
    # each family only against its own GCS listing, so a family whose file was never produced -- for a
    # sample, or for the whole run -- passes vacuously there; this catches both from the GCS listing
    # alone. The check ranges over the co-produced-data group (vet / ref_ranges / ploidy), which the Java
    # ingest emits together per sample, NOT over every configured prefix: vcf_header_lines_scratch is on a
    # different cadence and a headers-only phase produces only it. It is presence-activated -- dormant
    # unless a group member has files this run -- so that supported headers-only ingest passes while a
    # per-sample or whole-family gap within the data group still fails. The group is intersected with the
    # configured prefixes so a run that does not load one of them cannot false-positive on its absence.
    # Exact, so it gates all_loaded alongside completeness and cardinality (see
    # assess_cross_family_consistency).
    configured_families = set(superpartitioned_table_prefixes) | set(regular_table_prefixes)
    active_co_produced_families = [f for f in co_produced_families if f in configured_families]
    # If the header family is configured and at least one data family is present in this run,
    # headers and data are being loaded together (e.g. GvsQuickstartIntegration with load_vcf_headers=true).
    # Include headers in the cross-family check so a missing header Parquet cannot authorize deleting
    # an incomplete source set. Keep it excluded for a true headers-only run (no data families present).
    header_family = "vcf_header_lines_scratch"
    if header_family in configured_families and any(expected_by_family.get(f) for f in active_co_produced_families):
        active_co_produced_families.append(header_family)
    cross_family = assess_cross_family_consistency(expected_by_family, active_co_produced_families)

    # Per-sample cardinality consistency, applied only to regular tables that carry a UNIFORM
    # per-sample row count (ploidy). Tables whose per-sample count legitimately varies -- notably
    # vcf_header_lines_scratch, one row per header line with counts differing across samples -- are
    # excluded here (completeness above still covers them); enforcing the mode on them would
    # false-positive every sample as "deviating" and fail a valid ingest.
    cardinality = {}
    cardinality_ok = True
    for prefix in regular_table_prefixes:
        if prefix not in cardinality_table_prefixes:
            continue
        result = assess_cardinality(
            regular_counts.get(prefix, {}),
            expected_by_family.get(prefix, set()),
            expected_count=expected_ploidy_rows_per_sample,
        )
        result["backfill_caveat"] = PLOIDY_BACKFILL_CAVEAT
        cardinality[prefix] = result
        cardinality_ok = cardinality_ok and result["ok"]

    # Duplication and truncation screens on the usable superpartitioned families only (vet): same
    # families, same threshold, opposite sides of the median. ref_ranges is recorded as intentionally
    # unchecked. The truncation screen is the low-side mirror -- it catches a grossly under-rowed
    # partition, the one truncation source a row-count gate can produce cheaply (see
    # assess_truncation_screen for why an exact load-vs-source comparison would add nothing here).
    duplication = {}
    truncation = {}
    for family in duplication_screen_families:
        if family in superpartitioned_table_prefixes:
            duplication[family] = assess_duplication_screen(
                partition_rows, family, expected_by_family.get(family, set()),
                superpartitioned_table_prefixes, vet_duplication_threshold,
            )
            truncation[family] = assess_truncation_screen(
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

    duplication_flagged = any(d["outliers"] or d.get("singleton_flagged", False) for d in duplication.values())
    truncation_flagged = any(t["outliers"] or t.get("singleton_flagged", False) for t in truncation.values())

    return {
        "completeness_ok": completeness["ok"],
        "cardinality_ok": cardinality_ok,
        "cross_family_ok": cross_family["ok"],
        "duplication_flagged": duplication_flagged,
        "truncation_flagged": truncation_flagged,
        "allow_flagged_vet_loads": allow_flagged_vet_loads,
        "details": {
            "family_completeness": completeness,
            "cross_family_consistency": cross_family,
            "cardinality": cardinality,
            "duplication_screen": duplication,
            "truncation_screen": truncation,
            "duplication_unscreened": unscreened,
        },
    }
