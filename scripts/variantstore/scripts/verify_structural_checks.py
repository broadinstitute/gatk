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

The functions here ask a different question -- "how many rows?" -- through cheap BigQuery reads
that share no code with the loader predicate:

  * ``INFORMATION_SCHEMA.PARTITIONS.total_rows`` for the superpartitioned ``vet_%`` / ``ref_ranges_%``
    families (partition metadata is a flat ~10 MB regardless of callset size; ``partition_id`` *is*
    the ``sample_id`` because these tables are integer-range partitioned on ``sample_id`` with step
    1),
  * a per-sample ``COUNT(*)`` and ``COUNT(DISTINCT chromosome)`` on the regular ploidy table (scans
    two columns, ~100 MB on a 500k-sample callset), and
  * the row count each source Parquet file declares in its own footer, read by BigQuery through an
    external table over the files (VS-2032; see ``get_parquet_footer_row_counts``).

The crux of the independence is ``total_rows`` vs. ``total_logical_bytes``: a partition that is
present but partial or duplicated has a wrong row count, which these checks catch and the loader
predicate cannot.

Scope notes:
  * The footer comparison (``assess_footer_row_counts``, VS-2032) is the exact per-sample loss and
    duplication detector. For each ``(table, sample_id)`` it compares the sum of ``num_rows`` across
    the sample's Parquet footers against what BigQuery holds for it, and any disagreement, in either
    direction, flags the sample with no threshold involved. It covers ``vet``, ``ref_ranges`` and
    ploidy alike, so it replaces the ratio-to-median ``vet`` duplication screen VS-1989 shipped and
    gives ``ref_ranges`` the duplication detector it otherwise has none of: that family admits no
    count-ratio heuristic, since its row count tracks GQ-band transitions rather than genome length
    and legitimate samples reach many times the median.
  * It is free. BigQuery answers ``COUNT(*)`` grouped by ``_FILE_NAME`` from footer metadata, so no
    data column is read and nothing is billed however many files there are; it needs no Parquet
    reader, no image change and no client-side read path. Do not repeat the two older rationales for
    deferring it: that pyarrow cannot be had on musl (it has published ``musllinux_1_2_x86_64`` wheels
    since 20.0.0, and ``build_base.Dockerfile`` is stale on this point), or that the read I/O is
    prohibitive at AoU scale (a footer fetch is one or two ranged GETs, cents for ~1M objects, and
    moot anyway since BigQuery does the reading).
  * Its limit: ``num_rows`` is the file's own self-report, so the comparison is exact with respect to
    *loading* and not with respect to *generation*. The same file loaded twice reads footer N against
    partition 2N and is caught; a file written upstream with doubled rows reads footer 2N against
    partition 2N and passes as agreement. Truncation has the same shape: a file generated with too
    few rows agrees with its partition. That upstream case is surfaced, heuristically, by
    ``assess_truncation_screen``, a below-median ``vet`` screen; catching it exactly would need a
    gVCF-level count, which is out of scope. A partial *load* cannot occur: BigQuery load jobs are
    atomic (``load_parquet_to_bq.py`` loads each batch via ``load_table_from_uri`` with
    ``WRITE_APPEND``, and a failed job commits nothing).
  * Two other detectors were considered for ``ref_ranges`` and not adopted. ``COUNT(*)`` vs.
    ``COUNT(DISTINCT packed_ref_data)`` per sample is exact, but a full scan of the dataset's largest
    table rather than a metadata read. Comparing a sample's row count against its counterpart in the
    parent callset cannot work: a child callset is seeded by copying its parent, so a defect already
    present in the parent is copied forward identically and reads as agreement. Do not re-attempt
    either without re-reading this note.
  * Ploidy cardinality is not a ``ref_ranges`` detector, and must not be read as one. What it sees is
    a doubled *ploidy* file or load; because a sample's vet, ref_ranges and ploidy files are written
    in one pass (see ``DEFAULT_CO_PRODUCED_FAMILIES``), a flag is evidence about that pass rather than
    a readout of what landed in ``ref_ranges``. It is kept alongside the footer comparison because it
    sees the one thing footers cannot: a ploidy file written upstream with a repeated chromosome. And
    the reading holds only where a sample's ploidy rows were written at ingest; see
    ``PLOIDY_BACKFILL_CAVEAT``.
"""

import logging
import math
import statistics
from collections import defaultdict

try:
    from google.cloud import bigquery
except ImportError:
    bigquery = None  # type: ignore  # will fail at runtime if BigQuery calls are made without the package installed

# Reuse only the SQL-injection guard from the loader module. This is an input validator, not the
# loader's "what is already loaded?" oracle, so importing it does not compromise the independence
# these checks exist to provide.
from parse_and_group_files import _validate_table_prefixes


log = logging.getLogger(__name__)

# Families for which a per-sample row-count-to-median truncation screen produces a usable signal.
# Only ``vet`` qualifies (tight distribution); ``ref_ranges`` is deliberately excluded (see module
# docstring). Callers may override.
DEFAULT_TRUNCATION_SCREEN_FAMILIES = ["vet"]

# Families whose source Parquet footers are compared, per sample, against the rows BigQuery holds
# (VS-2032). Each is a directory directly under the Parquet output dir, which GenerateParquetFiles
# uploads to as ``<family>/[<table_number>/]<file>.parquet``. vcf_header_lines_scratch is not a member:
# its table is a scratch table drained by a later step, so its row count is not a stable target.
DEFAULT_FOOTER_CHECKED_FAMILIES = ["vet", "ref_ranges", "sample_chromosome_ploidy"]

# Name the external table over the source Parquet goes by within the footer query. Local to that one
# query job; nothing is created in any dataset.
FOOTER_EXTERNAL_TABLE = "source_parquet"

# Regular (non-superpartitioned) tables held to cardinality consistency by the cardinality check.
# The check is exact and per-sample: COUNT(*) must equal COUNT(DISTINCT chromosome), which a
# duplicated load cannot satisfy. Only ploidy qualifies, because only ploidy has a column whose values
# are unique per sample by construction; tables whose per-sample row count legitimately varies with
# nothing to compare it against (such as vcf_header_lines_scratch) are excluded entirely.
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

# Default ratio-to-median at or below whose reciprocal a vet sample is flagged as possibly truncated:
# the floor is baseline / threshold, so 1.6 here means 0.625x. It mirrors the 1.6x ceiling of the
# retired duplication screen, which the Foxtrot calibration measured on the high side only (2 samples
# of 540,545); nothing has been measured below the median, and a variant-count distribution has no
# reason to be symmetric in ratio space -- the upper tail is bounded by biology while the lower tail
# absorbs low-coverage samples, smaller callable fractions, and more aggressive GQ dropping. 1.6 is
# therefore a provisional default, not a derived one. Set to 0 to disable the truncation screen
# entirely (see run_structural_checks). VS-1989.
DEFAULT_VET_TRUNCATION_THRESHOLD = 1.6

# Sentinel value of the truncation threshold meaning "do not run the truncation screen at all".
TRUNCATION_SCREEN_DISABLED = 0

PLOIDY_BACKFILL_CAVEAT = (
    "This check covers the ploidy table's own rows and nothing else -- it is not a vet or ref_ranges "
    "duplication detector; the footer row count comparison is (VS-2032). Even for ploidy it is only "
    "meaningful where a sample's rows were written at ingest: a later backfill writes one row per "
    "chromosome by construction, so a backfilled sample satisfies "
    "COUNT(*) == COUNT(DISTINCT chromosome) whatever happened at load time. Do not read a green "
    "ploidy result as a duplication guarantee for callsets or sample ranges whose ploidy was "
    "backfilled (provenance is not captured today -- see VS-1989)."
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


def footer_source_uri(output_gcs_dir, family):
    """
    Return the external table source URI covering one family's Parquet under ``output_gcs_dir``.

    A trailing wildcard matches across ``/``, so ``<dir>/vet/*.parquet`` spans every ``vet/<NNN>/``
    subdirectory. Both halves of the scoping matter. The ``<family>/`` component keeps the URI out of
    ``<dir>/quarantine/`` and any sibling directory; a callset-wide wildcard has been observed to match
    a sibling callset and double-count. The ``.parquet`` suffix excludes ``.parquet.quarantined``
    objects and the verification diagnostics JSON.
    """
    if not output_gcs_dir.startswith("gs://"):
        raise ValueError(f"output_gcs_dir must be a gs:// URI (got {output_gcs_dir!r})")
    return f"{output_gcs_dir.rstrip('/')}/{family}/*.parquet"


def get_parquet_footer_row_counts(project_id, output_gcs_dir, families):
    """
    Return ``{gcs_uri: num_rows}`` for every Parquet file under each family's directory, as each
    file's own footer declares it.

    BigQuery reads the footers itself: ``COUNT(*)`` over an external table, grouped by the
    ``_FILE_NAME`` pseudo-column, is answered from footer metadata without touching a data column, so
    the query bills zero bytes however many files it covers. One query per family, because the
    families have different schemas and an external table takes one.

    Callers must pass only families that have files: BigQuery fails a query whose wildcard URI matches
    nothing.
    """
    if not families:
        return {}

    client = bigquery.Client(project=project_id)
    query = (f"SELECT _FILE_NAME AS file_name, COUNT(*) AS num_rows "
             f"FROM {FOOTER_EXTERNAL_TABLE} GROUP BY file_name")
    file_rows = {}
    for family in families:
        external = bigquery.ExternalConfig("PARQUET")
        external.source_uris = [footer_source_uri(output_gcs_dir, family)]
        job_config = bigquery.QueryJobConfig(table_definitions={FOOTER_EXTERNAL_TABLE: external})
        try:
            job = client.query(query, job_config=job_config)
            rows = list(job.result())
        except Exception as e:
            log.error(f"ERROR: Could not read Parquet footer row counts for {family} from "
                      f"{external.source_uris[0]}: {e}")
            raise
        for row in rows:
            file_rows[row.file_name] = row.num_rows
        # Logged because the zero cost is a measured property, not a guarantee: a change to how
        # BigQuery plans this query would show up here first.
        log.info(f"Read footer row counts for {len(rows)} {family} Parquet file(s) "
                 f"(bytes billed: {job.total_bytes_billed}, slot ms: {job.slot_millis})")
    return file_rows


def assess_footer_row_counts(file_rows, gcs_pairs_to_files, partition_rows, regular_counts,
                             footer_families, superpartitioned_table_prefixes,
                             regular_table_prefixes):
    """
    Compare, per ``(table, sample_id)``, the rows the sample's source Parquet footers declare against
    the rows BigQuery holds for it (VS-2032). Disagreement in either direction is exact: more rows in
    BigQuery is a duplicated load, fewer is a lost one. No threshold is involved.

    ``file_rows`` is ``{gcs_uri: num_rows}`` from ``get_parquet_footer_row_counts``. Only files in the
    GCS listing (``gcs_pairs_to_files``) are counted, so a file written under the output dir after
    discovery cannot inflate a sample's sum; such files are reported as ``unlisted_files`` and
    otherwise ignored. A listed file the footer query returned no row for -- unreadable, or a valid
    file with zero rows, which ``GROUP BY _FILE_NAME`` cannot tell apart -- cannot be counted, so its
    sample cannot be cleared: it is reported under ``unread`` and flagged alongside the mismatches.

    Ingest writes exactly one file per family per sample, so more than one listed file for a pair is
    flagged on its own under ``multiple_files``. Summing would otherwise pass a sample loaded once
    from each of two files, e.g. after a re-run at a different scatter width renamed its file.

    A sample with no BigQuery rows at all compares against 0. That is also a completeness failure,
    which fails the load outright; it is reported here too rather than special-cased.

    Returns an overall ``ok``, ``flagged_samples`` (every sample with a mismatch, unread file or
    multiple files in any family), and a per-family breakdown.
    """
    bigquery_rows = {(table_name, sample_id): total_rows or 0
                     for table_name, sample_id, total_rows in partition_rows}
    for prefix, counts in regular_counts.items():
        for sample_id, val in counts.items():
            bigquery_rows[(prefix, sample_id)] = _extract_count_info(val)[0] or 0

    per_family = {fam: {"ok": True, "pairs_checked": 0, "mismatched": [], "unread": [],
                        "multiple_files": []}
                  for fam in footer_families}
    listed = set()
    for (table_name, sample_id), files in sorted(gcs_pairs_to_files.items()):
        fam = family_for_table(table_name, superpartitioned_table_prefixes)
        if fam is None and table_name in regular_table_prefixes:
            fam = table_name
        if fam not in per_family:
            continue
        listed.update(files)
        result = per_family[fam]
        result["pairs_checked"] += 1
        if len(files) > 1:
            result["multiple_files"].append({"table": table_name, "sample_id": sample_id,
                                             "files": sorted(files)})
        unread = sorted(f for f in files if f not in file_rows)
        if unread:
            result["unread"].append({"table": table_name, "sample_id": sample_id, "files": unread})
            continue
        footer = sum(file_rows[f] for f in files)
        in_bigquery = bigquery_rows.get((table_name, sample_id), 0)
        if footer != in_bigquery:
            result["mismatched"].append({"table": table_name, "sample_id": sample_id,
                                         "footer_rows": footer, "bigquery_rows": in_bigquery})

    flagged = set()
    for result in per_family.values():
        result["ok"] = not result["mismatched"] and not result["unread"] and not result["multiple_files"]
        flagged.update(m["sample_id"] for m in result["mismatched"])
        flagged.update(u["sample_id"] for u in result["unread"])
        flagged.update(m["sample_id"] for m in result["multiple_files"])

    return {
        "ok": not flagged,
        "flagged_samples": sorted(flagged),
        "unlisted_files": sorted(set(file_rows) - listed),
        "per_family": per_family,
    }


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
        # Presence only. Whether a present partition holds the right number of rows is the footer
        # comparison's question (assess_footer_row_counts).
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


def assess_cardinality(counts, expected_samples):
    """
    Check that every expected sample is present in a cardinality-checked regular table, and that no
    sample holds more rows than it holds distinct chromosomes.

    ``SamplePloidyCreator`` writes at most one row per chromosome per sample, so
    ``COUNT(*) == COUNT(DISTINCT chromosome)`` is an exact duplication detector: a doubled load
    repeats chromosomes, and nothing legitimate does. That is the whole check, and it needs no
    reference count to compare against -- it holds for a single sample as readily as for a callset.

    Nothing here compares one sample's row count against another's. Those counts legitimately vary:
    a chromosome gets a ploidy row only where a non-PAR reference block on it survived ``drop_state``
    (``RefRangesCreator.apply``), so a sample's contig count depends on its data, not just its sex --
    every gVCF ingested for Foxtrot carries chrY, female participants included -- chrM is
    recorded or not, and an exome or BGE callset covers fewer contigs. Every screen built on that
    variation was either redundant with the exact comparison
    -- the 1.5x modal ceiling, which fired only where ``COUNT(*) > COUNT(DISTINCT chromosome)``
    already had -- or unsound: a floor below the callset mode reads a legitimately smaller contig set
    as a partial load, and on a small batch there is no mode worth measuring against in the first
    place.

    The gap that leaves: ploidy rows lost for *some* of a sample's chromosomes, with no repeats, is
    not detectable from row counts at all. The count is merely smaller, which is also what a
    legitimate karyotype looks like, so no threshold separates the two. Catching it needs the expected
    contig set for the sample's assay, which GVS does not record today. Whole-sample loss is still
    caught, by completeness. VS-1989.

    ``counts`` maps sample_id to ``(row_count, distinct_chromosomes)``, to
    ``{"count": int, "distinct": int}``, or to a raw row count. A value carrying no
    distinct-chromosome count degenerates to a presence check, the exact comparison being
    unavailable. Only ``expected_samples`` are assessed.
    """
    expected = set(expected_samples)
    present = {
        sid: counts[sid] for sid in expected
        if sid in counts and (_extract_count_info(counts[sid])[0] or 0) > 0
    }
    missing = sorted(expected - set(present))

    deviating = []
    for sid, val in sorted(present.items(), key=lambda d: d[0]):
        n, n_distinct = _extract_count_info(val)
        if n_distinct is not None and n > n_distinct:
            deviating.append({"sample_id": sid, "count": n, "distinct_chromosomes": n_distinct})

    values = [_extract_count_info(val)[0] for val in present.values()]

    return {
        "ok": not missing and not deviating,
        "min": min(values) if values else None,
        "max": max(values) if values else None,
        "distinct_samples": len(present),
        "total_rows": sum(values),
        "missing_samples": missing,
        "deviating_samples": deviating,
    }


def assess_truncation_screen(partition_rows, family, expected_samples,
                             superpartitioned_table_prefixes, threshold):
    """
    Flag samples in ``family`` whose ``total_rows`` is at most ``median / threshold`` -- a heuristic
    screen for a grossly truncated partition. Returns the flagged samples in ``flagged_samples``; the
    caller decides whether those samples' Parquet is quarantined from deletion (the default) or deleted
    with everything else (``allow_flagged_vet_loads``). A flag is scoped to the samples named here and
    never to the run.

    The default threshold is provisional. It mirrors the 1.6x ceiling of the retired duplication
    screen, which the Foxtrot calibration measured on the high side only, so the 0.625x floor is an
    assumption of symmetry rather than a measurement -- and a variant-count distribution has no reason
    to be symmetric in ratio space: the upper tail is bounded by biology while the lower tail absorbs
    low-coverage samples, smaller callable fractions, and more aggressive GQ dropping. Pass
    ``TRUNCATION_SCREEN_DISABLED`` to turn it off. See VS-1989.

    This is the only detector for the one truncation source the footer comparison misses: a Parquet
    file generated upstream with far fewer rows than its peers. Its footer count and its BigQuery
    count agree, so ``assess_footer_row_counts`` passes it. Truncation cannot arise from the load
    itself -- BigQuery load jobs are atomic (``load_parquet_to_bq.py`` loads each batch via
    ``load_table_from_uri`` with ``WRITE_APPEND``, and a failed job commits nothing). Catching the
    upstream case exactly would need a gVCF-level count, which is out of scope; a below-peers
    heuristic is what remains, and that is this.

    Only meaningful for families with a tight per-sample distribution (``vet``); callers must not apply
    it to ``ref_ranges``, whose row count varies by orders of magnitude between legitimate samples.
    Zero-row partitions are excluded here -- that is the completeness check's empty-partition case --
    so this screen judges only partitions that are present and non-empty.
    """
    # Disabled is reported rather than silently skipped: a caller reading details["truncation_screen"]
    # must be able to tell "nothing fell below the floor" apart from "nothing was looked at".
    if threshold == TRUNCATION_SCREEN_DISABLED:
        return {"family": family, "threshold": threshold, "disabled": True, "median": None,
                "baseline": None, "samples_screened": 0, "outliers": [], "singleton_flagged": False,
                "flagged_samples": []}

    # The baseline population is this run's expected samples, NOT every sample in the dataset, so a
    # sample loaded by an earlier run contributes no cohort consensus here. That is a property of the
    # batch rather than of the callset, and it is why a small incremental top-up degenerates: N=1 has
    # no cohort at all and flags unconditionally below, and N=2 takes its baseline from one of the two
    # values it saw. ``partition_rows`` is already unscoped (see get_partition_row_counts), so drawing
    # the baseline from the whole callset would cost no additional BigQuery read -- it is left scoped
    # pending a decision on mixed-assay datasets. VS-2027; operator workarounds are in
    # scripts/variantstore/parquet_loading/README.md.
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
                "samples_screened": 0, "outliers": [], "singleton_flagged": False,
                "flagged_samples": []}

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
        # The samples this screen objects to, and so the samples whose Parquet is quarantined rather
        # than deleted. A singleton load has no peers to compare against, so the screen cannot clear
        # the one sample it saw and names it instead of an outlier set.
        "flagged_samples": (sorted(rows_by_sample) if singleton_flagged
                            else sorted(o["sample_id"] for o in outliers)),
    }


def run_structural_checks(project_id, dataset_name, expected_by_family, gcs_pairs_to_files,
                          output_gcs_dir,
                          superpartitioned_table_prefixes=None, regular_table_prefixes=None,
                          vet_truncation_threshold=DEFAULT_VET_TRUNCATION_THRESHOLD,
                          allow_flagged_vet_loads=False,
                          truncation_screen_families=None,
                          cardinality_table_prefixes=None,
                          co_produced_families=None,
                          footer_checked_families=None):
    """
    Run the independent structural checks and return a result dict.

    ``expected_by_family`` maps a family name (a superpartitioned prefix such as ``vet`` /
    ``ref_ranges``, or a regular table name such as ``sample_chromosome_ploidy``) to the set of
    ``sample_id`` values GCS says should be loaded for it. ``gcs_pairs_to_files`` maps each
    ``(table_name, sample_id)`` in the GCS listing to its Parquet URIs, and ``output_gcs_dir`` is the
    directory those files were discovered under; the footer comparison needs both.

    ``regular_table_prefixes`` are the non-superpartitioned per-sample tables checked for completeness
    (every expected sample present with > 0 rows). ``cardinality_table_prefixes`` is the subset of
    those additionally held to a UNIFORM per-sample row count; it defaults to ploidy alone. A regular
    table whose per-sample row count legitimately varies -- ``vcf_header_lines_scratch`` writes one row
    per header line and samples differ -- must NOT be listed here, or every sample would read as
    "deviating" and fail a valid ingest; completeness alone covers it.

    Cardinality is assessed exactly, per sample, with no configurable expected count: COUNT(*) must
    equal COUNT(DISTINCT chromosome). See assess_cardinality for why no count-to-count comparison is
    made and what that leaves undetected.

    The returned dict carries flat booleans read shallowly downstream (``completeness_ok``,
    ``cardinality_ok``, ``cross_family_ok``, ``row_count_mismatch_flagged``, ``truncation_flagged``)
    plus a nested ``details`` block for humans and logs. ``completeness_ok``, ``cardinality_ok`` and
    ``cross_family_ok`` are the exact signals that gate ``all_loaded`` (and so the fail-loud abort).

    The footer comparison and the truncation screen never affect ``all_loaded``, and do not gate
    deletion for the run either: they name samples, whose Parquet the caller quarantines out of the
    bulk delete's reach so it survives for inspection, and the rest of the callset is deleted normally.
    ``flagged_samples`` is the union of both. ``row_count_mismatch_samples`` is the footer comparison's
    share on its own, because the two are not alike: the footer comparison is exact and cannot be
    waived, while the truncation screen is a heuristic that ``allow_flagged_vet_loads`` waives. The
    waiver is applied by the caller; it is carried through here unchanged so the log summary can say
    so. Pass ``TRUNCATION_SCREEN_DISABLED`` as ``vet_truncation_threshold`` to switch the truncation
    screen off.
    """
    # Reject a nonsensical ratio before any BigQuery read. A value <= 1 makes ordinary samples satisfy
    # the outlier condition (<= median / threshold), so the screen would flag half the callset. Fail
    # fast with an actionable message rather than silently flag everything. The disable sentinel lets
    # an operator who finds the uncalibrated screen flagging ordinary samples switch it off. VS-1989.
    try:
        threshold_ok = math.isfinite(vet_truncation_threshold) and vet_truncation_threshold > 1
    except TypeError:
        threshold_ok = False
    if not threshold_ok and vet_truncation_threshold != TRUNCATION_SCREEN_DISABLED:
        raise ValueError(
            f"vet_truncation_threshold must be a finite number > 1, or exactly "
            f"{TRUNCATION_SCREEN_DISABLED} to disable the screen (got {vet_truncation_threshold!r}); "
            "a ratio <= 1 flags ordinary samples."
        )

    if superpartitioned_table_prefixes is None:
        superpartitioned_table_prefixes = ["vet", "ref_ranges"]
    if regular_table_prefixes is None:
        regular_table_prefixes = ["sample_chromosome_ploidy"]
    if truncation_screen_families is None:
        truncation_screen_families = DEFAULT_TRUNCATION_SCREEN_FAMILIES
    if cardinality_table_prefixes is None:
        cardinality_table_prefixes = DEFAULT_CARDINALITY_TABLE_PREFIXES
    if co_produced_families is None:
        co_produced_families = DEFAULT_CO_PRODUCED_FAMILIES
    if footer_checked_families is None:
        footer_checked_families = DEFAULT_FOOTER_CHECKED_FAMILIES

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

    # Per-sample cardinality consistency, applied only to regular tables with a per-chromosome unique
    # row (ploidy). Tables whose rows carry no such uniqueness -- notably vcf_header_lines_scratch, one
    # row per header line -- are excluded here, there being nothing exact to compare their counts
    # against; completeness above still covers them.
    cardinality = {}
    cardinality_ok = True
    for prefix in regular_table_prefixes:
        if prefix not in cardinality_table_prefixes:
            continue
        result = assess_cardinality(
            regular_counts.get(prefix, {}),
            expected_by_family.get(prefix, set()),
        )
        result["backfill_caveat"] = PLOIDY_BACKFILL_CAVEAT
        cardinality[prefix] = result
        cardinality_ok = cardinality_ok and result["ok"]

    # Footer comparison (VS-2032), exact and per sample. Limited to the configured families, and of
    # those to the ones with files this run: a wildcard URI that matches nothing fails the query.
    configured_footer_families = [f for f in footer_checked_families if f in configured_families]
    footer_families_with_files = [f for f in configured_footer_families if expected_by_family.get(f)]
    file_rows = get_parquet_footer_row_counts(project_id, output_gcs_dir, footer_families_with_files)
    footer = assess_footer_row_counts(
        file_rows, gcs_pairs_to_files, partition_rows, regular_counts,
        configured_footer_families, superpartitioned_table_prefixes, regular_table_prefixes,
    )

    # Truncation screen on the usable superpartitioned families only (vet). ref_ranges is excluded: its
    # row count varies by orders of magnitude between legitimate samples (see the module docstring).
    truncation = {}
    for family in truncation_screen_families:
        if family in superpartitioned_table_prefixes:
            truncation[family] = assess_truncation_screen(
                partition_rows, family, expected_by_family.get(family, set()),
                superpartitioned_table_prefixes, vet_truncation_threshold,
            )
    truncation_flagged_samples = sorted(
        {sid for result in truncation.values() for sid in result["flagged_samples"]}
    )

    return {
        "completeness_ok": completeness["ok"],
        "cardinality_ok": cardinality_ok,
        "cross_family_ok": cross_family["ok"],
        "row_count_mismatch_flagged": not footer["ok"],
        "truncation_flagged": bool(truncation_flagged_samples),
        "flagged_samples": sorted(set(footer["flagged_samples"]) | set(truncation_flagged_samples)),
        "row_count_mismatch_samples": footer["flagged_samples"],
        "truncation_flagged_samples": truncation_flagged_samples,
        "allow_flagged_vet_loads": allow_flagged_vet_loads,
        "details": {
            "family_completeness": completeness,
            "cross_family_consistency": cross_family,
            "cardinality": cardinality,
            "footer_row_counts": footer,
            "truncation_screen": truncation,
        },
    }
