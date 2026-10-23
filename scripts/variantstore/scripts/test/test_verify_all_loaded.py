"""
Unit tests for verify_all_loaded.py, focused on the VS-1989 composition change, its Part-2 split of
the deletion gate, and the VS-2032 replacement of the heuristic vet duplication screen by an exact
footer row count comparison.

The exact structural checks (family completeness, ploidy cardinality, cross-family consistency) block
all_loaded -- the factual "is the load complete?" signal fail-loud aborts on -- even when the
shared-predicate presence check is fully satisfied. Neither the footer comparison nor the vet
truncation screen touches all_loaded or blocks deletion for the run; both instead name the samples
whose Parquet is quarantined out of the bulk delete's reach. They differ in waivability:
--allow-flagged-vet-loads waives the truncation screen (heuristic) and does NOT waive the footer
comparison (exact). The results JSON must carry both predicates, the quarantine work list, the
structural booleans, and the pre-existing keys the WDL already reads.

The BigQuery-touching helpers (get_already_loaded_tables_and_sample_ids, run_structural_checks) are
patched; the real GCS-path parser runs against realistic fixture paths.
"""

import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import verify_all_loaded

# Six fixture files: two samples (1, 2) across vet, ref_ranges, and sample_chromosome_ploidy.
FIXTURE_FILES = [
    "gs://b/vet/vet_001_1_input_vcf_0_S1.vcf.gz.parquet",
    "gs://b/vet/vet_001_2_input_vcf_0_S2.vcf.gz.parquet",
    "gs://b/ref_ranges/ref_ranges_001_1_input_vcf_0_S1.vcf.gz.parquet",
    "gs://b/ref_ranges/ref_ranges_001_2_input_vcf_0_S2.vcf.gz.parquet",
    "gs://b/sample_chromosome_ploidy/sample_chromosome_ploidy_1_S1.parquet",
    "gs://b/sample_chromosome_ploidy/sample_chromosome_ploidy_2_S2.parquet",
]

# The directory those fixture paths live under -- what the footer comparison derives its external
# table URIs from, and what the WDL passes as --output-gcs-dir.
OUTPUT_GCS_DIR = "gs://b"

ALL_PAIRS = {
    ("vet_001", 1), ("vet_001", 2),
    ("ref_ranges_001", 1), ("ref_ranges_001", 2),
    ("sample_chromosome_ploidy", 1), ("sample_chromosome_ploidy", 2),
}


def _structural(completeness_ok=True, cardinality_ok=True, cross_family_ok=True,
                mismatch_samples=(), truncation_samples=(), allow_flagged_vet_loads=False):
    """A run_structural_checks return value with the keys verify_all_loaded consumes.

    The two flag booleans are derived from the sample lists rather than passed separately, because in
    a real result they come from the same pass over the same data: a raised flag always names at least
    one sample, and a named sample always raises the flag.
    """
    mismatch_samples = list(mismatch_samples)
    truncation_samples = list(truncation_samples)
    return {
        "completeness_ok": completeness_ok,
        "cardinality_ok": cardinality_ok,
        "cross_family_ok": cross_family_ok,
        "flagged_samples": sorted(set(mismatch_samples) | set(truncation_samples)),
        "row_count_mismatch_samples": mismatch_samples,
        "truncation_flagged_samples": truncation_samples,
        "row_count_mismatch_flagged": bool(mismatch_samples),
        "truncation_flagged": bool(truncation_samples),
        "allow_flagged_vet_loads": allow_flagged_vet_loads,
        "details": {
            "family_completeness": {"ok": completeness_ok, "per_family": {}},
            "cross_family_consistency": {"ok": cross_family_ok, "union_size": 0, "per_family": {}},
            "cardinality": {},
            "footer_row_counts": {
                "ok": not mismatch_samples,
                "flagged_samples": mismatch_samples,
                "unlisted_files": [],
                "per_family": {},
            },
            "truncation_screen": {},
        },
    }


class VerifyAllLoadedTestBase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gcs_list = os.path.join(self.tmp, "gcs_files.txt")
        with open(self.gcs_list, "w") as f:
            f.write("\n".join(FIXTURE_FILES) + "\n")
        self.out_dir = os.path.join(self.tmp, "out")

    def _run(self, loaded_pairs, structural, allow_flagged_vet_loads=False, **kwargs):
        kwargs.setdefault("output_gcs_dir", OUTPUT_GCS_DIR)
        with patch("verify_all_loaded.get_already_loaded_tables_and_sample_ids",
                   return_value=loaded_pairs), \
             patch("verify_all_loaded.run_structural_checks",
                   return_value=structural) as mock_struct:
            self.mock_struct = mock_struct
            return verify_all_loaded.verify_all_loaded(
                project_id="proj",
                dataset_name="ds",
                gcs_files_list=self.gcs_list,
                output_dir=self.out_dir,
                allow_flagged_vet_loads=allow_flagged_vet_loads,
                **kwargs,
            )

    def _written_json(self):
        with open(os.path.join(self.out_dir, "verification_results.json")) as f:
            return json.load(f)

    def _quarantine_list(self):
        """The quarantine work list as written to disk -- what the WDL move step actually reads."""
        with open(os.path.join(self.out_dir, "quarantine_files.txt")) as f:
            return [line.strip() for line in f if line.strip()]


class TestHappyPath(VerifyAllLoadedTestBase):
    def test_all_loaded_and_json_keys(self):
        r = self._run(set(ALL_PAIRS), _structural())

        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["safe_to_delete_parquet"])
        # Pre-existing keys the WDL already reads must be preserved.
        self.assertEqual(r["total_files"], 6)
        self.assertEqual(r["loaded_files"], 6)
        self.assertEqual(r["missing_files"], 0)
        self.assertEqual(r["unmatched_files"], 0)
        self.assertIsNone(r["missing_files_list"])
        # New structural booleans (read shallowly by the WDL).
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(r["family_completeness_ok"])
        self.assertTrue(r["ploidy_cardinality_ok"])
        self.assertTrue(r["cross_family_consistency_ok"])
        self.assertFalse(r["row_count_mismatch_flagged"])
        self.assertFalse(r["vet_truncation_flagged"])
        self.assertIn("structural_checks", r)
        # The JSON on disk matches the returned dict.
        self.assertEqual(self._written_json(), r)

    def test_expected_by_family_passed_to_structural_checks(self):
        self._run(set(ALL_PAIRS), _structural())
        expected_by_family = self.mock_struct.call_args[0][2]
        self.assertEqual(set(expected_by_family["vet"]), {1, 2})
        self.assertEqual(set(expected_by_family["ref_ranges"]), {1, 2})
        self.assertEqual(set(expected_by_family["sample_chromosome_ploidy"]), {1, 2})

    def test_output_gcs_dir_is_passed_to_structural_checks(self):
        # The footer comparison addresses <output_gcs_dir>/<family>/*.parquet, so this argument is the
        # whole basis of what gets compared; dropping it would silently compare nothing.
        self._run(set(ALL_PAIRS), _structural())
        self.assertEqual(self.mock_struct.call_args[0][4], OUTPUT_GCS_DIR)


class TestExactChecksGateAllLoaded(VerifyAllLoadedTestBase):
    """The exact checks fail all_loaded (and so safe_to_delete_parquet) even when the shared predicate
    sees every pair present."""

    def test_cardinality_failure_blocks_all_loaded(self):
        r = self._run(set(ALL_PAIRS), _structural(cardinality_ok=False))
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["safe_to_delete_parquet"])
        self.assertFalse(r["structural_checks_ok"])
        self.assertFalse(r["ploidy_cardinality_ok"])
        self.assertEqual(r["missing_files"], 0)

    def test_completeness_failure_blocks_all_loaded(self):
        r = self._run(set(ALL_PAIRS), _structural(completeness_ok=False))
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["safe_to_delete_parquet"])
        self.assertFalse(r["family_completeness_ok"])

    def test_cross_family_failure_blocks_all_loaded(self):
        # A sample present in some families but absent from another (a partial upload) fails all_loaded
        # even though every pair the listing names is present in BigQuery.
        r = self._run(set(ALL_PAIRS), _structural(cross_family_ok=False))
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["safe_to_delete_parquet"])
        self.assertFalse(r["structural_checks_ok"])
        self.assertFalse(r["cross_family_consistency_ok"])
        self.assertEqual(r["missing_files"], 0)


class TestFlaggedSamplesQuarantineRatherThanBlockDeletion(VerifyAllLoadedTestBase):
    """A footer mismatch or a truncation flag is orthogonal to load completeness and scoped to the
    samples it names: all_loaded stays True (so the task succeeds), the bulk delete stays authorized
    for the rest of the callset, and the flagged samples' Parquet is listed for quarantine instead."""

    def test_footer_mismatch_quarantines_without_blocking_deletion(self):
        r = self._run(set(ALL_PAIRS), _structural(mismatch_samples=[1]))
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(r["row_count_mismatch_flagged"])
        # The delete proceeds for the rest of the callset -- the whole point of quarantining.
        self.assertTrue(r["safe_to_delete_parquet"])
        self.assertEqual(r["quarantine_samples"], [1])
        self.assertEqual(r["quarantine_sample_count"], 1)
        # The count GvsImportGenomes' abort condition reads.
        self.assertEqual(r["quarantine_row_count_mismatch_sample_count"], 1)
        self.assertFalse(r["quarantine_waived"])

    def test_quarantine_covers_the_flagged_sample_across_every_family(self):
        # A flag says this sample's ingest is unexplained, so all of its Parquet is held back, not just
        # the flagging family's -- a sample is only re-ingestable from a complete set of its files.
        r = self._run(set(ALL_PAIRS), _structural(mismatch_samples=[1]))
        self.assertEqual(r["quarantine_files"], 3)
        self.assertEqual(sorted(self._quarantine_list()), sorted([
            "gs://b/vet/vet_001_1_input_vcf_0_S1.vcf.gz.parquet",
            "gs://b/ref_ranges/ref_ranges_001_1_input_vcf_0_S1.vcf.gz.parquet",
            "gs://b/sample_chromosome_ploidy/sample_chromosome_ploidy_1_S1.parquet",
        ]))

    def test_unflagged_samples_are_not_quarantined(self):
        self._run(set(ALL_PAIRS), _structural(mismatch_samples=[1]))
        self.assertNotIn("S2", "".join(self._quarantine_list()))

    def test_truncation_flag_quarantines_without_blocking_deletion(self):
        r = self._run(set(ALL_PAIRS), _structural(truncation_samples=[2]))
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(r["vet_truncation_flagged"])
        self.assertTrue(r["safe_to_delete_parquet"])
        self.assertEqual(r["quarantine_samples"], [2])
        self.assertEqual(r["quarantine_files"], 3)

    def test_both_checks_flagging_different_samples_quarantines_both(self):
        r = self._run(set(ALL_PAIRS),
                      _structural(mismatch_samples=[1], truncation_samples=[2]))
        self.assertTrue(r["safe_to_delete_parquet"])
        self.assertEqual(r["quarantine_samples"], [1, 2])
        self.assertEqual(r["quarantine_files"], 6)

    def test_clean_run_writes_an_empty_quarantine_list(self):
        # The list is always written so the WDL move step can be unconditional.
        r = self._run(set(ALL_PAIRS), _structural())
        self.assertEqual(r["quarantine_files"], 0)
        self.assertEqual(r["quarantine_samples"], [])
        self.assertEqual(self._quarantine_list(), [])
        self.assertTrue(os.path.exists(os.path.join(self.out_dir, "quarantine_files.txt")))

    def test_unquarantinable_flagged_sample_blocks_deletion(self):
        # A flagged sample_id with no GCS file behind it: the quarantine would protect nothing, so the
        # bulk delete must not run, or it destroys exactly the Parquet the checks asked to keep.
        r = self._run(set(ALL_PAIRS), _structural(mismatch_samples=[999]))
        self.assertTrue(r["all_loaded"])
        self.assertFalse(r["safe_to_delete_parquet"])
        self.assertEqual(r["quarantine_uncovered_samples"], [999])


class TestWaiverCoversTheTruncationScreenOnly(VerifyAllLoadedTestBase):
    """VS-2032: --allow-flagged-vet-loads waives the uncalibrated truncation screen and nothing else.
    The footer comparison is exact -- the sample's rows in BigQuery are not the rows its Parquet holds
    -- so there is no sound reading of that finding under which its Parquet should be deleted."""

    def test_truncation_flag_waived_quarantines_nothing(self):
        r = self._run(
            set(ALL_PAIRS),
            _structural(truncation_samples=[1], allow_flagged_vet_loads=True),
            allow_flagged_vet_loads=True,
        )
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["vet_truncation_flagged"])
        self.assertTrue(r["safe_to_delete_parquet"])
        self.assertTrue(r["quarantine_waived"])
        self.assertEqual(r["quarantine_files"], 0)
        self.assertEqual(self._quarantine_list(), [])

    def test_footer_mismatch_is_quarantined_even_when_waived(self):
        r = self._run(
            set(ALL_PAIRS),
            _structural(mismatch_samples=[1], allow_flagged_vet_loads=True),
            allow_flagged_vet_loads=True,
        )
        self.assertTrue(r["row_count_mismatch_flagged"])
        self.assertEqual(r["quarantine_samples"], [1])
        self.assertEqual(r["quarantine_files"], 3)
        self.assertEqual(len(self._quarantine_list()), 3)
        # The abort count the WDL reads is unaffected by the waiver, so a waived run with a mismatch
        # still stops under parquet_fail_on_quarantine.
        self.assertEqual(r["quarantine_row_count_mismatch_sample_count"], 1)

    def test_waiver_separates_the_two_findings_in_one_run(self):
        # Sample 1 mismatched, sample 2 truncated, waiver on: only sample 1's Parquet is held back.
        r = self._run(
            set(ALL_PAIRS),
            _structural(mismatch_samples=[1], truncation_samples=[2],
                        allow_flagged_vet_loads=True),
            allow_flagged_vet_loads=True,
        )
        self.assertEqual(r["quarantine_samples"], [1])
        self.assertEqual(r["quarantine_files"], 3)
        self.assertNotIn("S2", "".join(self._quarantine_list()))


class TestStructuralDetailCapped(VerifyAllLoadedTestBase):
    def test_large_lists_capped_in_json(self):
        # A failure that produces a huge per-sample list must not bloat the results JSON: the embedded
        # list is capped and a sibling *_total records the true length.
        cap = verify_all_loaded.STRUCTURAL_DETAIL_LIST_CAP
        huge = list(range(cap + 500))
        structural = _structural(completeness_ok=False)
        structural["details"]["family_completeness"]["per_family"]["vet"] = {
            "ok": False, "expected": len(huge), "present": 0,
            "missing_samples": huge, "empty_partition_samples": [],
        }

        r = self._run(set(ALL_PAIRS), structural)

        vet = r["structural_checks"]["family_completeness"]["per_family"]["vet"]
        self.assertEqual(len(vet["missing_samples"]), cap)
        self.assertEqual(vet["missing_samples_total"], len(huge))
        # The on-disk JSON carries the same bounded copy.
        self.assertEqual(self._written_json(), r)

    def test_small_lists_untouched(self):
        structural = _structural(completeness_ok=False)
        structural["details"]["family_completeness"]["per_family"]["vet"] = {
            "ok": False, "expected": 3, "present": 1,
            "missing_samples": [4], "empty_partition_samples": [3],
        }

        r = self._run(set(ALL_PAIRS), structural)

        vet = r["structural_checks"]["family_completeness"]["per_family"]["vet"]
        self.assertEqual(vet["missing_samples"], [4])
        self.assertNotIn("missing_samples_total", vet)


class TestComputeStructuralChecksOk(unittest.TestCase):
    """The exact structural signal that feeds all_loaded: completeness, cardinality and cross-family
    consistency only. The footer comparison and the truncation screen are deliberately excluded --
    they gate safe_to_delete_parquet, not all_loaded."""

    def _ok(self, structural):
        return verify_all_loaded.compute_structural_checks_ok(structural)

    def test_all_ok_passes(self):
        self.assertTrue(self._ok(_structural()))

    def test_completeness_failure_gates(self):
        self.assertFalse(self._ok(_structural(completeness_ok=False)))

    def test_cardinality_failure_gates(self):
        self.assertFalse(self._ok(_structural(cardinality_ok=False)))

    def test_cross_family_failure_gates(self):
        self.assertFalse(self._ok(_structural(cross_family_ok=False)))

    def test_quarantining_findings_do_not_affect_exact_signal(self):
        # Including the footer mismatch: exact though it is, it names samples rather than condemning
        # the run, and a failing all_loaded would abort before anything could be quarantined.
        self.assertTrue(self._ok(_structural(mismatch_samples=[1])))
        self.assertTrue(self._ok(_structural(truncation_samples=[1])))


class TestComputeAllLoaded(unittest.TestCase):
    """The factual load-complete gate (what fail-loud aborts on): True only when nothing is
    missing/unmatched and the exact structural checks pass. The quarantining findings never enter."""

    def test_true_when_everything_clean(self):
        self.assertTrue(verify_all_loaded.compute_all_loaded(set(), [], True))

    def test_missing_pair_blocks(self):
        self.assertFalse(verify_all_loaded.compute_all_loaded({("vet_001", 2)}, [], True))

    def test_unmatched_file_blocks(self):
        self.assertFalse(verify_all_loaded.compute_all_loaded(set(), ["gs://b/weird.parquet"], True))

    def test_structural_failure_blocks_even_when_files_all_present(self):
        self.assertFalse(verify_all_loaded.compute_all_loaded(set(), [], False))


def _quarantine(waived=False, samples=(), paths=(), uncovered_samples=(), mismatch_samples=None):
    """A compute_quarantine return value.

    mismatch_samples defaults to all of samples, matching the common fixture where the footer
    comparison is what flagged them; pass it explicitly for a truncation-only flag.
    """
    return {"waived": waived, "samples": list(samples), "paths": list(paths),
            "uncovered_samples": list(uncovered_samples),
            "row_count_mismatch_samples":
                list(samples if mismatch_samples is None else mismatch_samples)}


class TestComputeQuarantine(unittest.TestCase):
    """Resolving flagged sample_ids to the Parquet paths held back from the bulk delete."""

    PAIRS_TO_FILES = {
        ("vet_001", 1): ["gs://b/vet/vet_001_1_a.parquet"],
        ("ref_ranges_001", 1): ["gs://b/ref_ranges/ref_ranges_001_1_a.parquet"],
        ("sample_chromosome_ploidy", 1): ["gs://b/sample_chromosome_ploidy/p_1.parquet"],
        ("vet_001", 2): ["gs://b/vet/vet_001_2_a.parquet"],
        ("ref_ranges_001", 2): ["gs://b/ref_ranges/ref_ranges_001_2_a.parquet"],
        ("sample_chromosome_ploidy", 2): ["gs://b/sample_chromosome_ploidy/p_2.parquet"],
    }

    def _q(self, structural, allow=False):
        return verify_all_loaded.compute_quarantine(structural, self.PAIRS_TO_FILES, allow)

    def test_no_flags_quarantines_nothing(self):
        q = self._q(_structural())
        self.assertEqual(q, _quarantine())

    def test_flagged_sample_pulls_all_its_families(self):
        q = self._q(_structural(mismatch_samples=[1]))
        self.assertEqual(q["samples"], [1])
        self.assertEqual(q["paths"], sorted([
            "gs://b/vet/vet_001_1_a.parquet",
            "gs://b/ref_ranges/ref_ranges_001_1_a.parquet",
            "gs://b/sample_chromosome_ploidy/p_1.parquet",
        ]))
        self.assertEqual(q["uncovered_samples"], [])

    def test_paths_are_sorted_and_deduplicated_by_pair(self):
        q = self._q(_structural(mismatch_samples=[1, 2]))
        self.assertEqual(q["paths"], sorted(q["paths"]))
        self.assertEqual(len(q["paths"]), 6)

    def test_multiple_files_per_pair_all_quarantined(self):
        # A sample ingested from several input VCFs has several Parquet files per table.
        pairs = {("vet_001", 1): ["gs://b/vet/a.parquet", "gs://b/vet/b.parquet"]}
        q = verify_all_loaded.compute_quarantine(
            _structural(mismatch_samples=[1]), pairs, False)
        self.assertEqual(q["paths"], ["gs://b/vet/a.parquet", "gs://b/vet/b.parquet"])

    def test_waiver_drops_the_truncation_flag_before_resolving_paths(self):
        q = self._q(_structural(truncation_samples=[1]), allow=True)
        self.assertEqual(q, _quarantine(waived=True))

    def test_waiver_does_not_drop_a_footer_mismatch(self):
        q = self._q(_structural(mismatch_samples=[1]), allow=True)
        self.assertEqual(q["samples"], [1])
        self.assertEqual(len(q["paths"]), 3)
        self.assertEqual(q["row_count_mismatch_samples"], [1])
        # waived still reports the operator's setting, even though it changed nothing here.
        self.assertTrue(q["waived"])

    def test_flagged_sample_with_no_files_is_reported_uncovered(self):
        q = self._q(_structural(mismatch_samples=[1, 999]))
        self.assertEqual(q["samples"], [1, 999])
        self.assertEqual(q["uncovered_samples"], [999])
        self.assertEqual(len(q["paths"]), 3)

    def test_missing_flagged_samples_keys_are_tolerated(self):
        # Defensive: a structural result predating the keys must not crash the verifier.
        q = verify_all_loaded.compute_quarantine({}, self.PAIRS_TO_FILES, False)
        self.assertEqual(q, _quarantine())


class TestQuarantineMismatchCountInResults(VerifyAllLoadedTestBase):
    """The results JSON key the WDL abort condition reads, end to end through verify_all_loaded."""

    def test_truncation_only_run_quarantines_but_reports_no_mismatch(self):
        r = self._run(set(ALL_PAIRS), _structural(truncation_samples=[1]))
        self.assertTrue(r["all_loaded"])
        self.assertEqual(r["quarantine_sample_count"], 1)
        self.assertEqual(r["quarantine_row_count_mismatch_sample_count"], 0)
        self.assertEqual(len(self._quarantine_list()), 3)

    def test_clean_run_reports_zero(self):
        r = self._run(set(ALL_PAIRS), _structural())
        self.assertEqual(r["quarantine_row_count_mismatch_sample_count"], 0)

    def test_waived_truncation_run_reports_zero(self):
        r = self._run(set(ALL_PAIRS), _structural(truncation_samples=[1],
                                                  allow_flagged_vet_loads=True),
                      allow_flagged_vet_loads=True)
        self.assertEqual(r["quarantine_row_count_mismatch_sample_count"], 0)

    def test_waived_run_with_a_mismatch_still_reports_it(self):
        r = self._run(set(ALL_PAIRS), _structural(mismatch_samples=[1],
                                                  allow_flagged_vet_loads=True),
                      allow_flagged_vet_loads=True)
        self.assertEqual(r["quarantine_row_count_mismatch_sample_count"], 1)


class TestQuarantineMismatchAttribution(unittest.TestCase):
    """Only a footer mismatch aborts the run (GvsImportGenomes' parquet_fail_on_quarantine), so the
    mismatch subset of the quarantined samples has to be reported separately and accurately -- the
    truncation threshold is not calibrated, and must not be able to fail a completed ingest."""

    PAIRS_TO_FILES = {
        ("vet_001", 1): ["gs://b/vet/vet_001_1_a.parquet"],
        ("vet_001", 2): ["gs://b/vet/vet_001_2_a.parquet"],
    }

    def _q(self, structural, allow=False):
        return verify_all_loaded.compute_quarantine(structural, self.PAIRS_TO_FILES, allow)

    def test_clean_run_attributes_nothing(self):
        self.assertEqual(self._q(_structural())["row_count_mismatch_samples"], [])

    def test_footer_mismatch_is_attributed(self):
        q = self._q(_structural(mismatch_samples=[1]))
        self.assertEqual(q["samples"], [1])
        self.assertEqual(q["row_count_mismatch_samples"], [1])

    def test_truncation_only_flag_is_quarantined_but_not_attributed(self):
        # The sample is still held back from the delete; it just does not abort the run.
        q = self._q(_structural(truncation_samples=[1]))
        self.assertEqual(q["samples"], [1])
        self.assertEqual(q["row_count_mismatch_samples"], [])

    def test_mixed_findings_attribute_only_the_mismatched_sample(self):
        q = self._q(_structural(mismatch_samples=[2], truncation_samples=[1]))
        self.assertEqual(q["samples"], [1, 2])
        self.assertEqual(q["row_count_mismatch_samples"], [2])

    def test_waiver_attributes_the_mismatch_it_could_not_waive(self):
        q = self._q(_structural(mismatch_samples=[1], truncation_samples=[2]), allow=True)
        self.assertEqual(q["samples"], [1])
        self.assertEqual(q["row_count_mismatch_samples"], [1])

    def test_attribution_never_names_an_unquarantined_sample(self):
        # The abort must not be able to fire for a sample whose Parquet was not actually held back.
        for structural, allow in [
            (_structural(mismatch_samples=[1], truncation_samples=[2]), False),
            (_structural(mismatch_samples=[1], truncation_samples=[2]), True),
            (_structural(truncation_samples=[1, 2]), True),
        ]:
            q = self._q(structural, allow)
            self.assertLessEqual(set(q["row_count_mismatch_samples"]), set(q["samples"]))

    def test_missing_keys_attribute_nothing(self):
        self.assertEqual(
            verify_all_loaded.compute_quarantine(
                {}, self.PAIRS_TO_FILES, False)["row_count_mismatch_samples"],
            [],
        )


class TestComputeSafeToDeleteParquet(unittest.TestCase):
    """The bulk deletion gate: all_loaded, and every flagged sample quarantinable."""

    def _safe(self, all_loaded, quarantine):
        return verify_all_loaded.compute_safe_to_delete_parquet(all_loaded, quarantine)

    def test_requires_all_loaded(self):
        # An incomplete load is never safe to delete, whatever the quarantine says.
        self.assertFalse(self._safe(False, _quarantine()))
        self.assertFalse(self._safe(False, _quarantine(waived=True)))

    def test_clean_load_is_safe(self):
        self.assertTrue(self._safe(True, _quarantine()))

    def test_flagged_but_quarantinable_load_is_still_safe_to_delete(self):
        # The change of policy: a flag no longer blocks the run's delete, because the flagged samples'
        # Parquet is moved out of its reach first.
        self.assertTrue(self._safe(True, _quarantine(samples=[1], paths=["gs://b/vet/a.parquet"])))

    def test_waived_flags_are_safe_to_delete(self):
        self.assertTrue(self._safe(True, _quarantine(waived=True)))

    def test_unquarantinable_flagged_sample_blocks(self):
        self.assertFalse(self._safe(True, _quarantine(samples=[999], uncovered_samples=[999])))


class TestSharedPredicateStillGates(VerifyAllLoadedTestBase):
    def test_missing_pair_blocks_even_when_structural_ok(self):
        loaded = set(ALL_PAIRS) - {("vet_001", 2)}
        r = self._run(loaded, _structural())
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["safe_to_delete_parquet"])
        self.assertEqual(r["missing_files"], 1)
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(os.path.exists(r["missing_files_list"]))


class TestDescribeIncompleteReasons(unittest.TestCase):
    """The fail-loud operator message names every not-all-loaded cause. Only the exact gating checks
    appear; neither quarantining finding fails all_loaded, so neither must surface here."""

    def _reasons(self, **over):
        base = {
            "missing_files": 0,
            "unmatched_files": 0,
            "family_completeness_ok": True,
            "ploidy_cardinality_ok": True,
            "cross_family_consistency_ok": True,
        }
        base.update(over)
        return verify_all_loaded.describe_incomplete_reasons(base)

    def test_clean_results_yield_no_reasons(self):
        self.assertEqual(self._reasons(), [])

    def test_missing_and_unmatched_named(self):
        reasons = self._reasons(missing_files=3, unmatched_files=2)
        self.assertIn("3 file(s) not yet loaded", reasons)
        self.assertTrue(any("2 file(s) could not be parsed" in r for r in reasons))

    def test_completeness_and_cardinality_named(self):
        reasons = self._reasons(family_completeness_ok=False, ploidy_cardinality_ok=False)
        self.assertTrue(any("family completeness" in r for r in reasons))
        self.assertTrue(any("ploidy cardinality" in r for r in reasons))

    def test_cross_family_named(self):
        reasons = self._reasons(cross_family_consistency_ok=False)
        self.assertTrue(any("cross-family consistency" in r for r in reasons))

    def test_quarantining_findings_do_not_appear(self):
        # These gate safe_to_delete_parquet, never all_loaded, so they must not appear as a fail-loud
        # reason even when flagged. The footer mismatch belongs here too: it is exact, but it aborts
        # via the WDL's quarantine gate, not by claiming the load is incomplete.
        self.assertEqual(
            self._reasons(row_count_mismatch_flagged=True, vet_truncation_flagged=True), [])


def _log_details(**over):
    """The details dict _log_structural_summary walks, with every section it reads present."""
    details = {
        "family_completeness": {"ok": True, "per_family": {}},
        "cross_family_consistency": {"ok": True, "union_size": 2, "per_family": {}},
        "cardinality": {},
        "footer_row_counts": {"ok": True, "unlisted_files": [], "per_family": {}},
        "truncation_screen": {},
    }
    details.update(over)
    return details


class TestLogStructuralSummary(unittest.TestCase):
    def test_truncation_outlier_logs_baseline(self):
        structural = {
            "allow_flagged_vet_loads": False,
            "details": _log_details(
                cardinality={
                    "sample_chromosome_ploidy": {
                        "ok": True,
                        "min": 23,
                        "max": 24,
                        "distinct_samples": 2,
                    }
                },
                truncation_screen={
                    "vet": {
                        "threshold": 1.6,
                        "median": 150,
                        "baseline": 200,
                        "samples_screened": 2,
                        "outliers": [{"sample_id": 1, "rows": 100, "ratio": 0.5}],
                        "singleton_flagged": False,
                    }
                },
            ),
        }
        with self.assertLogs("verify_all_loaded", level="INFO") as cm:
            verify_all_loaded._log_structural_summary(structural)

        logs = "\n".join(cm.output)
        self.assertIn("1 sample(s) <= baseline/1.6 (200) (quarantining their Parquet)", logs)
        self.assertIn("one row per chromosome for every one of them (min=23, max=24 rows/sample)", logs)

    def test_disabled_truncation_screen_is_logged_as_disabled_not_clean(self):
        # A switched-off screen must not read like a screen that looked and found nothing -- that is
        # the whole hazard of giving the uncalibrated screen an off switch.
        structural = {
            "allow_flagged_vet_loads": False,
            "details": _log_details(
                truncation_screen={
                    "vet": {
                        "threshold": 0,
                        "disabled": True,
                        "median": None,
                        "baseline": None,
                        "samples_screened": 0,
                        "outliers": [],
                        "singleton_flagged": False,
                    }
                },
            ),
        }
        with self.assertLogs("verify_all_loaded", level="INFO") as cm:
            verify_all_loaded._log_structural_summary(structural)

        logs = "\n".join(cm.output)
        self.assertIn("screen disabled", logs)
        self.assertIn("no low-side check was performed", logs)
        self.assertNotIn("no samples <=", logs)

    def test_clean_footer_comparison_reports_what_it_compared(self):
        # "agree" is only meaningful alongside the number of pairs behind it: a footer query that
        # compared nothing would otherwise log exactly like one that compared the whole callset.
        structural = {
            "allow_flagged_vet_loads": False,
            "details": _log_details(
                footer_row_counts={
                    "ok": True,
                    "unlisted_files": [],
                    "per_family": {
                        "ref_ranges": {"ok": True, "pairs_checked": 2, "mismatched": [], "unread": []},
                    },
                },
            ),
        }
        with self.assertLogs("verify_all_loaded", level="INFO") as cm:
            verify_all_loaded._log_structural_summary(structural)

        logs = "\n".join(cm.output)
        self.assertIn("[footer] ref_ranges: footer and BigQuery row counts agree for 2", logs)

    def test_footer_mismatch_is_logged_as_unwaivable(self):
        structural = {
            # Even with the waiver set, the footer line must not read as waived.
            "allow_flagged_vet_loads": True,
            "details": _log_details(
                footer_row_counts={
                    "ok": False,
                    "unlisted_files": ["gs://b/ref_ranges/001/stray.parquet"],
                    "per_family": {
                        "ref_ranges": {
                            "ok": False,
                            "pairs_checked": 2,
                            "mismatched": [{"table": "ref_ranges_001", "sample_id": 1,
                                            "footer_rows": 100, "bigquery_rows": 200}],
                            "unread": [],
                        },
                    },
                },
            ),
        }
        with self.assertLogs("verify_all_loaded", level="INFO") as cm:
            verify_all_loaded._log_structural_summary(structural)

        logs = "\n".join(cm.output)
        self.assertIn("not waivable", logs)
        self.assertIn("1 Parquet file(s) under the output dir", logs)
        self.assertIn("stray.parquet", logs)


class TestTruncationThresholdPlumbing(VerifyAllLoadedTestBase):
    """VS-1989: the truncation ratio reaches run_structural_checks as its own argument. VS-2032
    removed the duplication threshold it used to share a default with, so these now also stand as the
    check that the surviving threshold is still plumbed at all."""

    def test_default_is_forwarded(self):
        self._run(set(ALL_PAIRS), _structural())
        kwargs = self.mock_struct.call_args.kwargs
        self.assertEqual(kwargs["vet_truncation_threshold"],
                         verify_all_loaded.DEFAULT_VET_TRUNCATION_THRESHOLD)

    def test_explicit_threshold_is_forwarded(self):
        self._run(set(ALL_PAIRS), _structural(), vet_truncation_threshold=1.25)
        self.assertEqual(self.mock_struct.call_args.kwargs["vet_truncation_threshold"], 1.25)

    def test_disable_sentinel_is_forwarded(self):
        self._run(set(ALL_PAIRS), _structural(),
                  vet_truncation_threshold=verify_all_loaded.TRUNCATION_SCREEN_DISABLED)
        self.assertEqual(self.mock_struct.call_args.kwargs["vet_truncation_threshold"],
                         verify_all_loaded.TRUNCATION_SCREEN_DISABLED)

    def test_duplication_threshold_is_gone(self):
        # VS-2032 retired the screen; a resurrected argument would silently do nothing.
        self.assertFalse(hasattr(verify_all_loaded, "DEFAULT_VET_DUPLICATION_THRESHOLD"))
        self._run(set(ALL_PAIRS), _structural())
        self.assertNotIn("vet_duplication_threshold", self.mock_struct.call_args.kwargs)


class TestCommandLineInterface(VerifyAllLoadedTestBase):
    def _main_kwargs(self, extra_argv):
        """Drive main() through the real argument parser and report what it passed downstream."""
        argv = ["verify_all_loaded.py", "--project-id", "p", "--dataset-name", "d",
                "--gcs-files-list", self.gcs_list, "--output-dir", self.out_dir] + extra_argv
        stub = {"all_loaded": True, "safe_to_delete_parquet": True, "quarantine_files": 0}
        with patch.object(sys, "argv", argv), \
             patch("verify_all_loaded.verify_all_loaded", return_value=stub) as mock_verify:
            verify_all_loaded.main()
        return mock_verify.call_args.kwargs

    def test_output_gcs_dir_is_forwarded(self):
        kwargs = self._main_kwargs(["--output-gcs-dir", "gs://b/data"])
        self.assertEqual(kwargs["output_gcs_dir"], "gs://b/data")

    def test_output_gcs_dir_is_required(self):
        # Required rather than defaulted: there is no value that could be guessed from the listing
        # without risking a footer comparison that silently addresses the wrong prefix.
        with self.assertRaises(SystemExit) as cm:
            self._main_kwargs([])
        self.assertEqual(cm.exception.code, 2)

    def test_truncation_threshold_default(self):
        kwargs = self._main_kwargs(["--output-gcs-dir", "gs://b"])
        self.assertEqual(kwargs["vet_truncation_threshold"],
                         verify_all_loaded.DEFAULT_VET_TRUNCATION_THRESHOLD)

    def test_truncation_threshold_accepts_the_disable_sentinel(self):
        kwargs = self._main_kwargs(["--output-gcs-dir", "gs://b",
                                    "--vet-truncation-threshold", "0"])
        self.assertEqual(kwargs["vet_truncation_threshold"], 0)

    def test_waiver_flag_is_forwarded(self):
        kwargs = self._main_kwargs(["--output-gcs-dir", "gs://b", "--allow-flagged-vet-loads"])
        self.assertTrue(kwargs["allow_flagged_vet_loads"])

    def test_duplication_threshold_flag_is_rejected(self):
        # The WDL passed this until VS-2032; an image/WDL version skew must fail loudly at argparse
        # rather than quietly ignoring a threshold the operator thinks they set.
        with self.assertRaises(SystemExit) as cm:
            self._main_kwargs(["--output-gcs-dir", "gs://b", "--vet-duplication-threshold", "2.0"])
        self.assertEqual(cm.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
