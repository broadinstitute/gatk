#!/usr/bin/env python3
"""
Verify that all Parquet files in GCS have been loaded to BigQuery.

Two layers of verification run here, and only together do they gate deletion of the source Parquet:

1. Shared-predicate presence check. Parses every GCS path into a (table_name, sample_id) pair and
   compares against get_already_loaded_tables_and_sample_ids -- the same predicate the loader uses to
   decide what to skip. This is cheap and catches whole-sample absence, but because it asks the
   loader's own question it can never contradict the loader: it is blind to a partition that is
   present-but-partial or present-but-duplicated (VS-1989).

2. Independent structural checks (verify_structural_checks.py). These ask "how many rows?" via
   INFORMATION_SCHEMA.PARTITIONS.total_rows and a per-sample COUNT(*) on the ploidy table -- signals
   the loader predicate never consults -- so they can catch partial loads and duplication. Family
   completeness and ploidy cardinality are exact checks: they gate all_loaded, the factual "is the
   load complete?" signal that fail-loud aborts on. The vet duplication and truncation screens are
   heuristics -- they never affect all_loaded, and instead gate the separate safe_to_delete_parquet
   predicate, which is what authorizes deleting the source Parquet. A screen flag blocks that deletion
   by default and is waived only with --allow-flagged-vet-loads, so a flagged-but-complete load
   succeeds and simply retains its Parquet rather than failing the import.

Neither layer requires or consults a parquet_load_status tracking table.
"""

import argparse
import json
import logging
import os
import sys
from collections import defaultdict

from parse_and_group_files import (
    parse_table_and_sample_id_from_file_path,
    get_already_loaded_tables_and_sample_ids,
)
from verify_structural_checks import (
    family_for_table,
    run_structural_checks,
    DEFAULT_VET_DUPLICATION_THRESHOLD,
)

log = logging.getLogger(__name__)

# Cap on the per-sample ID lists embedded in the results JSON. On a large-callset failure these lists
# (missing / empty / deviating samples, duplication outliers) could otherwise hold hundreds of
# thousands of entries, bloating verification_results.json -- which the WDL re-parses in full on every
# one of its read_json calls -- and the Cromwell logs. The gating booleans and the human-readable
# counts logged by _log_structural_summary are computed from the full, uncapped structural result;
# only the copy written to disk is bounded.
STRUCTURAL_DETAIL_LIST_CAP = 1000


def _cap_structural_detail_lists(obj, cap=STRUCTURAL_DETAIL_LIST_CAP):
    """
    Return a copy of the structural-checks detail block with any per-sample list longer than ``cap``
    truncated to its first ``cap`` entries. Where a list is truncated a sibling ``<key>_total`` key
    records its true length, so the count survives even though the full enumeration does not (the
    authoritative set lives in BigQuery). Recurses through the nested dicts (per_family, cardinality,
    duplication_screen); list elements themselves (ints or small ``{sample_id, ...}`` dicts) are left
    untouched.
    """
    if isinstance(obj, dict):
        capped = {}
        for key, value in obj.items():
            if isinstance(value, list) and len(value) > cap:
                capped[key] = value[:cap]
                capped[f"{key}_total"] = len(value)
            else:
                capped[key] = _cap_structural_detail_lists(value, cap)
        return capped
    return obj


def _log_structural_summary(structural):
    """Log a human-readable summary of the independent structural checks."""
    details = structural["details"]

    for family, fam in sorted(details["family_completeness"]["per_family"].items()):
        if fam["ok"]:
            log.info(f"  [completeness] {family}: {fam['present']}/{fam['expected']} expected samples present")
        else:
            log.error(
                f"  [completeness] {family}: {len(fam['missing_samples'])} missing, "
                f"{len(fam['empty_partition_samples'])} present-but-empty "
                f"(expected {fam['expected']})"
            )

    cross = details.get("cross_family_consistency")
    if cross is not None:
        if cross["ok"]:
            log.info(f"  [cross-family] all families share the same {cross['union_size']} expected sample(s)")
        else:
            for family, fam in sorted(cross["per_family"].items()):
                if fam["missing_samples"]:
                    log.error(
                        f"  [cross-family] {family}: {len(fam['missing_samples'])} sample(s) present in "
                        f"other families but absent here (of {cross['union_size']} across all families)"
                    )

    for table, card in sorted(details["cardinality"].items()):
        ref_desc = f"{card.get('reference_count')} rows/sample ({card.get('reference_source', 'mode')})"
        if card["ok"]:
            log.info(
                f"  [cardinality] {table}: {card['distinct_samples']} samples, all with {ref_desc}"
            )
        else:
            log.error(
                f"  [cardinality] {table}: expected {ref_desc}, observed mode={card['mode']} "
                f"min={card['min']} max={card['max']}; "
                f"{len(card['missing_samples'])} missing, {len(card['deviating_samples'])} off-reference"
            )

    for family, screen in sorted(details["duplication_screen"].items()):
        outliers = screen["outliers"]
        level = "  [duplication]"
        if not outliers:
            log.info(f"{level} {family}: no samples >= {screen['threshold']}x median ({screen['median']})")
        elif structural["allow_flagged_vet_loads"]:
            log.warning(f"{level} {family}: {len(outliers)} sample(s) >= {screen['threshold']}x median (warning only; --allow-flagged-vet-loads set)")
        else:
            log.error(f"{level} {family}: {len(outliers)} sample(s) >= {screen['threshold']}x median (blocks Parquet deletion)")

    for family, screen in sorted(details.get("truncation_screen", {}).items()):
        outliers = screen["outliers"]
        level = "  [truncation]"
        if not outliers:
            log.info(f"{level} {family}: no samples <= median/{screen['threshold']} ({screen['median']})")
        elif structural["allow_flagged_vet_loads"]:
            log.warning(f"{level} {family}: {len(outliers)} sample(s) <= median/{screen['threshold']} (warning only; --allow-flagged-vet-loads set)")
        else:
            log.error(f"{level} {family}: {len(outliers)} sample(s) <= median/{screen['threshold']} (blocks Parquet deletion)")

    unscreened = details["duplication_unscreened"]
    if unscreened["families"]:
        log.info(f"  [duplication] not screened for {unscreened['families']}: {unscreened['reason']}")


def compute_structural_checks_ok(structural):
    """
    Reduce a run_structural_checks result to the exact structural boolean that feeds all_loaded.

    Only the exact checks contribute: family completeness (every expected sample present and
    non-empty), per-sample ploidy cardinality, and cross-family consistency (a sample present in some
    families must be present in all -- catches the cross-family gap completeness judges vacuously). The
    duplication and truncation screens are heuristics and deliberately do NOT gate all_loaded --
    all_loaded is the factual "is the load complete?" signal that fail-loud aborts on, and a heuristic
    flag does not make a complete load incomplete. The screens instead gate the separate
    safe_to_delete_parquet predicate (see compute_safe_to_delete_parquet).
    """
    return (
        structural["completeness_ok"]
        and structural["cardinality_ok"]
        and structural["cross_family_ok"]
    )


def compute_all_loaded(missing_pairs, unmatched_files, structural_checks_ok):
    """
    The factual "is the load complete?" gate, and the signal fail-loud aborts on. Returns True only
    when no (table, sample_id) pair the loader should have produced is missing from BigQuery, no GCS
    file was left unmatched, and the exact structural checks pass (see compute_structural_checks_ok).

    all_loaded is deliberately NOT the deletion gate: the heuristic vet screens are excluded here so a
    flagged-but-complete load still reads as loaded and its task succeeds (retaining its Parquet)
    rather than aborting the import. Whether it is safe to delete the source Parquet is the separate
    compute_safe_to_delete_parquet predicate, which layers the screen policy on top of this. Kept a
    pure function of its inputs so both predicates stay trivially testable in isolation.
    """
    return (
        len(missing_pairs) == 0
        and len(unmatched_files) == 0
        and structural_checks_ok
    )


def compute_safe_to_delete_parquet(all_loaded, structural, allow_flagged_vet_loads):
    """
    The deletion gate -- the single predicate DeleteParquetFiles is downstream of. Returns True only
    when the load is factually complete (all_loaded) AND no heuristic vet screen objects, unless the
    operator has explicitly waived the screens with allow_flagged_vet_loads.

    Gating the irreversible delete on a screen flag by default is the safe direction: a screen false
    positive here only over-retains Parquet (cheap and reversible), whereas trusting a flagged load and
    deleting its source is not. Waiving the screens is therefore an opt-out (--allow-flagged-vet-loads),
    so the default never deletes anything a screen objected to. Kept a pure function so the
    deletion-authorizing logic is unit-tested in isolation -- a regression here authorizes an
    irreversible delete.
    """
    if not all_loaded:
        return False
    if allow_flagged_vet_loads:
        return True
    return not (structural["duplication_flagged"] or structural["truncation_flagged"])


def describe_incomplete_reasons(results):
    """
    Human-readable reasons a verification came back not-all-loaded, for the fail-loud error line.

    Kept a pure function of the results dict so the operator-facing diagnosis is unit-tested: this line
    is what an operator reads when a run has aborted. Only the exact checks can make all_loaded false
    (missing/unmatched files, family completeness, ploidy cardinality, cross-family consistency); the
    vet duplication and truncation screens never fail all_loaded -- they gate safe_to_delete_parquet
    instead -- so they are deliberately absent here.
    """
    reasons = []
    missing_count = results.get("missing_files", 0) or 0
    unmatched_count = results.get("unmatched_files", 0) or 0
    if missing_count:
        reasons.append(f"{missing_count} file(s) not yet loaded")
    if unmatched_count:
        reasons.append(f"{unmatched_count} file(s) could not be parsed or matched to a table/sample_id")
    if not results.get("family_completeness_ok", True):
        reasons.append("family completeness check failed (missing or empty partitions)")
    if not results.get("ploidy_cardinality_ok", True):
        reasons.append("ploidy cardinality check failed (missing or off-reference samples)")
    if not results.get("cross_family_consistency_ok", True):
        reasons.append("cross-family consistency check failed (a sample present in some families is absent from another)")
    return reasons


def verify_all_loaded(project_id, dataset_name, gcs_files_list, output_dir,
                      superpartitioned_table_prefixes=None, regular_table_prefixes=None,
                      vet_duplication_threshold=DEFAULT_VET_DUPLICATION_THRESHOLD,
                      allow_flagged_vet_loads=False,
                      expected_ploidy_rows_per_sample=None):
    """
    Compare GCS-derived (table_name, sample_id) pairs against what is actually
    present in BigQuery to find loads that are missing, then run independent
    row-count-based structural checks that the loader-shared predicate cannot make.

    Args:
        project_id: BigQuery project ID
        dataset_name: BigQuery dataset name
        gcs_files_list: Path to file listing all GCS Parquet URIs
        output_dir: Directory to write results
        superpartitioned_table_prefixes: Prefixes for superpartitioned tables (default: ["vet", "ref_ranges"])
        regular_table_prefixes: Prefixes for regular tables (default: ["sample_chromosome_ploidy"])
        vet_duplication_threshold: Ratio-to-median above which a vet sample is flagged (default 1.6)
        allow_flagged_vet_loads: If False (default), a vet duplication- or truncation-screen flag
            blocks deletion of the source Parquet (the load still succeeds -- all_loaded stays factual
            -- and the Parquet is retained). If True, the screens are waived and deletion may proceed
            despite a flag.
        expected_ploidy_rows_per_sample: If set, the exact per-sample ploidy row count to validate
            against (e.g. 24 for WGS) instead of the callset mode; leave unset to infer from the data.

    Returns:
        Dictionary with verification results
    """
    if superpartitioned_table_prefixes is None:
        superpartitioned_table_prefixes = ["vet", "ref_ranges"]
    if regular_table_prefixes is None:
        regular_table_prefixes = ["sample_chromosome_ploidy"]

    os.makedirs(output_dir, exist_ok=True)

    # Read GCS file list
    with open(gcs_files_list) as f:
        gcs_files = [line.strip() for line in f if line.strip()]

    log.info(f"Found {len(gcs_files)} files in GCS")

    # Parse each GCS file path to determine its (table_name, sample_id).
    # Files whose path cannot be parsed are counted as unmatched and reported.
    gcs_pairs_to_files = {}   # (table_name, sample_id) -> list of file paths
    unmatched_files = []

    for file_path in gcs_files:
        table_name, sample_id = parse_table_and_sample_id_from_file_path(
            file_path, superpartitioned_table_prefixes, regular_table_prefixes
        )
        if table_name is None:
            unmatched_files.append(file_path)
            continue
        key = (table_name, sample_id)
        gcs_pairs_to_files.setdefault(key, []).append(file_path)

    if unmatched_files:
        log.warning(f"Could not parse {len(unmatched_files)} file path(s); they will be excluded from verification:")
        for p in unmatched_files[:20]:
            log.warning(f"  {p}")
        if len(unmatched_files) > 20:
            log.warning(f"  ... and {len(unmatched_files) - 20} more")

    total_files = len(gcs_files)
    all_gcs_pairs = set(gcs_pairs_to_files.keys())
    log.info(f"Parsed {len(all_gcs_pairs)} unique (table_name, sample_id) pairs from GCS file list")

    # Query BigQuery for all already-loaded (table_name, sample_id) pairs.
    log.info("Querying BigQuery for loaded data...")
    loaded_pairs = get_already_loaded_tables_and_sample_ids(
        project_id, dataset_name,
        superpartitioned_table_prefixes=superpartitioned_table_prefixes,
        regular_table_prefixes=regular_table_prefixes,
    )

    # Determine which GCS (table_name, sample_id) pairs are not yet in BigQuery.
    missing_pairs = all_gcs_pairs - loaded_pairs

    # Count files as loaded/missing based on their pair's presence in BigQuery.
    loaded_files_count = sum(
        len(files) for pair, files in gcs_pairs_to_files.items() if pair not in missing_pairs
    )
    missing_files_count = total_files - loaded_files_count - len(unmatched_files)

    log.info("Results:")
    log.info(f"  Total files in GCS:       {total_files}")
    log.info(f"  Unmatched (unparseable):  {len(unmatched_files)}")
    log.info(f"  Files with loaded pairs:  {loaded_files_count}")
    log.info(f"  Files with missing pairs: {missing_files_count}")
    log.info(f"  Missing (table, sample_id) pairs: {len(missing_pairs)}")

    # --- Independent structural checks (VS-1989) --------------------------------------------------
    # Derive, per family, the sample_ids GCS says should be loaded, then run row-count-based checks
    # that share no code with the loader's get_already_loaded_tables_and_sample_ids predicate.
    #
    # This independence has a boundary: expected_by_family is sourced from all_gcs_pairs, which is the
    # same GCS file listing (gcs_files_list, produced upstream by DiscoverParquetFiles) that the loader
    # itself scatters over. Within that listing the checks are now cross-family consistent: a sample
    # present in some co-produced data families but absent from another, AND an entirely absent
    # co-produced family, are both caught by assess_cross_family_consistency. It ranges over the group the
    # Java ingest emits together per sample (vet / ref_ranges / ploidy) rather than over every configured
    # prefix, so an absent group member is compared as an empty set -- either gap is a partial upload
    # rather than passing vacuously in the family that lacks the data -- while the check stays presence-
    # activated: it is dormant unless a group member has files this run, so a supported headers-only ingest
    # (which produces only vcf_header_lines_scratch, deliberately not a group member) is not falsely
    # reported incomplete. What remains outside the listing is the whole-sample gap: a sample for which
    # Parquet generation never produced output for ANY family was never in the listing at all, so it is in
    # no family's expected set and no cross-family union, and stays invisible to every check here -- as
    # opposed to a sample whose files were produced but never loaded into BigQuery, which the checks do
    # catch. Unlike a whole family, whose identity is known from the co-produced group, a never-produced
    # sample's identity is unknowable without a non-GCS source. Closing that residual gap needs an expected-sample source that
    # does not derive from GCS: specifically this run's input sample set (the ingest FOFN). Note that
    # sample_info cannot serve this wholesale -- it accumulates every sample ever ingested, including
    # ones since withdrawn or deleted, so comparing against it in bulk would false-positive samples that
    # are legitimately absent from this run; the expected set has to be scoped to the FOFN for this run.
    # Tracked as a follow-up in VS-2016 (family-completeness independence); see also VS-1989.
    expected_by_family = defaultdict(set)
    for table_name, sample_id in all_gcs_pairs:
        family = family_for_table(table_name, superpartitioned_table_prefixes)
        if family is not None:
            expected_by_family[family].add(sample_id)
        elif table_name in regular_table_prefixes:
            expected_by_family[table_name].add(sample_id)

    log.info("Running independent structural checks...")
    structural = run_structural_checks(
        project_id, dataset_name, expected_by_family,
        superpartitioned_table_prefixes=superpartitioned_table_prefixes,
        regular_table_prefixes=regular_table_prefixes,
        vet_duplication_threshold=vet_duplication_threshold,
        allow_flagged_vet_loads=allow_flagged_vet_loads,
        expected_ploidy_rows_per_sample=expected_ploidy_rows_per_sample,
    )
    _log_structural_summary(structural)

    # Two separate predicates, both pure helpers so the logic is unit-tested in isolation:
    #  * all_loaded -- factual "is the load complete?", from the exact checks only; fail-loud aborts on
    #    it, so a heuristic screen flag must NOT enter here or a complete load would wrongly abort.
    #  * safe_to_delete_parquet -- the deletion gate, which layers the vet screen policy on top of
    #    all_loaded: a flag blocks deletion unless allow_flagged_vet_loads waives the screens.
    structural_checks_ok = compute_structural_checks_ok(structural)
    all_loaded = compute_all_loaded(missing_pairs, unmatched_files, structural_checks_ok)
    safe_to_delete_parquet = compute_safe_to_delete_parquet(
        all_loaded, structural, allow_flagged_vet_loads
    )

    # Write list of missing file paths if there are any
    missing_files_list_path = None
    if missing_pairs:
        missing_files_list_path = f"{output_dir}/missing_files.txt"
        missing_file_paths = sorted(
            path
            for pair in sorted(missing_pairs)
            for path in gcs_pairs_to_files.get(pair, [])
        )
        with open(missing_files_list_path, 'w') as f:
            for path in missing_file_paths:
                f.write(f"{path}\n")
        log.error(f"Wrote {len(missing_file_paths)} missing file path(s) to {missing_files_list_path}")
        log.error(f"Missing (table_name, sample_id) pairs:")
        for pair in sorted(missing_pairs)[:20]:
            log.error(f"  {pair}")
        if len(missing_pairs) > 20:
            log.error(f"  ... and {len(missing_pairs) - 20} more")

    results_dict = {
        "all_loaded": all_loaded,
        # The deletion gate DeleteParquetFiles is downstream of: all_loaded AND no unwaived vet-screen
        # flag. Distinct from all_loaded so a flagged-but-complete load succeeds yet retains Parquet.
        "safe_to_delete_parquet": safe_to_delete_parquet,
        "total_files": total_files,
        "loaded_files": loaded_files_count,
        "missing_files": missing_files_count,
        "missing_files_list": missing_files_list_path,
        "unmatched_files": len(unmatched_files),
        # Flat structural-check booleans, read shallowly by the WDL (VS-1989).
        "structural_checks_ok": structural_checks_ok,
        "family_completeness_ok": structural["completeness_ok"],
        "ploidy_cardinality_ok": structural["cardinality_ok"],
        "cross_family_consistency_ok": structural["cross_family_ok"],
        "vet_duplication_flagged": structural["duplication_flagged"],
        "vet_truncation_flagged": structural["truncation_flagged"],
        # Full per-check detail for humans and logs, with per-sample lists bounded so a large-callset
        # failure cannot bloat this file (which the WDL re-parses on every read_json call).
        "structural_checks": _cap_structural_detail_lists(structural["details"]),
    }

    results_file = f"{output_dir}/verification_results.json"
    with open(results_file, 'w') as f:
        json.dump(results_dict, f, indent=2)

    log.info(f"Results written to {results_file}")
    return results_dict


def main():
    parser = argparse.ArgumentParser(
        description="Verify all Parquet files have been loaded to BigQuery"
    )
    parser.add_argument("--project-id", required=True, help="BigQuery project ID")
    parser.add_argument("--dataset-name", required=True, help="BigQuery dataset name")
    parser.add_argument(
        "--gcs-files-list",
        required=True,
        help="File containing list of all GCS Parquet URIs"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write verification results"
    )
    parser.add_argument(
        "--superpartitioned-table-prefixes",
        nargs="+",
        default=["vet", "ref_ranges"],
        help="Table prefixes for superpartitioned tables (default: vet ref_ranges)"
    )
    parser.add_argument(
        "--regular-table-prefixes",
        nargs="+",
        default=["sample_chromosome_ploidy"],
        help="Table prefixes for regular (non-superpartitioned) tables (default: sample_chromosome_ploidy)"
    )
    parser.add_argument(
        "--vet-duplication-threshold",
        type=float,
        default=DEFAULT_VET_DUPLICATION_THRESHOLD,
        help=(
            "Ratio-to-median above which a vet sample is flagged as a possible duplicate "
            f"(default: {DEFAULT_VET_DUPLICATION_THRESHOLD})"
        )
    )
    parser.add_argument(
        "--allow-flagged-vet-loads",
        action="store_true",
        help=(
            "Permit deleting the source Parquet even when the vet duplication or truncation screen "
            "flags a sample. By default a screen flag blocks deletion (the load still succeeds and its "
            "Parquet is retained); pass this to waive the screens and allow deletion anyway."
        )
    )
    parser.add_argument(
        "--expected-ploidy-rows-per-sample",
        type=int,
        default=None,
        help=(
            "Exact per-sample ploidy row count to validate against (e.g. 24 for WGS) instead of the "
            "callset mode. Leave unset to infer the reference from the data."
        )
    )

    args = parser.parse_args()

    results = verify_all_loaded(
        project_id=args.project_id,
        dataset_name=args.dataset_name,
        gcs_files_list=args.gcs_files_list,
        output_dir=args.output_dir,
        superpartitioned_table_prefixes=args.superpartitioned_table_prefixes,
        regular_table_prefixes=args.regular_table_prefixes,
        vet_duplication_threshold=args.vet_duplication_threshold,
        allow_flagged_vet_loads=args.allow_flagged_vet_loads,
        expected_ploidy_rows_per_sample=args.expected_ploidy_rows_per_sample,
    )

    if results["all_loaded"]:
        log.info("✓ SUCCESS: All files have been loaded!")
    else:
        reasons = describe_incomplete_reasons(results)
        if reasons:
            log.error("✗ INCOMPLETE: " + "; ".join(reasons))
        else:
            # Fallback in case all_loaded is False but no counts are available.
            log.error("✗ INCOMPLETE: Verification failed for unknown reasons")

        if results.get("missing_files_list"):
            log.error(f"  See missing files list: {results['missing_files_list']}")
        if results.get("unmatched_files_list"):
            log.error(f"  See unmatched files list: {results['unmatched_files_list']}")

    if not results["all_loaded"]:
        sys.exit(1)

    return 0


if __name__ == "__main__":
    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format='%(levelname)s - %(message)s')
    sys.exit(main())
