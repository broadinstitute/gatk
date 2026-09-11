"""
Unit tests for verify_all_loaded.py, focused on the VS-1989 composition change and its Part-2 split of
the deletion gate: the exact structural checks (family completeness, ploidy cardinality) block
all_loaded -- the factual "is the load complete?" signal fail-loud aborts on -- even when the
shared-predicate presence check is fully satisfied; while the heuristic vet duplication/truncation
screens never touch all_loaded and instead gate the separate safe_to_delete_parquet predicate, blocking
deletion by default and waived only under allow_flagged_vet_loads. The results JSON must carry both
predicates plus the structural booleans and the pre-existing keys the WDL already reads.

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

ALL_PAIRS = {
    ("vet_001", 1), ("vet_001", 2),
    ("ref_ranges_001", 1), ("ref_ranges_001", 2),
    ("sample_chromosome_ploidy", 1), ("sample_chromosome_ploidy", 2),
}


def _structural(completeness_ok=True, cardinality_ok=True, cross_family_ok=True,
                duplication_flagged=False, truncation_flagged=False, allow_flagged_vet_loads=False):
    """A run_structural_checks return value with the keys verify_all_loaded consumes."""
    return {
        "completeness_ok": completeness_ok,
        "cardinality_ok": cardinality_ok,
        "cross_family_ok": cross_family_ok,
        "duplication_flagged": duplication_flagged,
        "truncation_flagged": truncation_flagged,
        "allow_flagged_vet_loads": allow_flagged_vet_loads,
        "details": {
            "family_completeness": {"ok": completeness_ok, "per_family": {}},
            "cross_family_consistency": {"ok": cross_family_ok, "union_size": 0, "per_family": {}},
            "cardinality": {},
            "duplication_screen": {},
            "truncation_screen": {},
            "duplication_unscreened": {"families": ["ref_ranges"], "reason": "not screened"},
        },
    }


class VerifyAllLoadedTestBase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.gcs_list = os.path.join(self.tmp, "gcs_files.txt")
        with open(self.gcs_list, "w") as f:
            f.write("\n".join(FIXTURE_FILES) + "\n")
        self.out_dir = os.path.join(self.tmp, "out")

    def _run(self, loaded_pairs, structural, allow_flagged_vet_loads=False,
             expected_ploidy_rows_per_sample=None):
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
                expected_ploidy_rows_per_sample=expected_ploidy_rows_per_sample,
            )

    def _written_json(self):
        with open(os.path.join(self.out_dir, "verification_results.json")) as f:
            return json.load(f)


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
        self.assertFalse(r["vet_duplication_flagged"])
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

    def test_expected_ploidy_override_threaded_to_structural_checks(self):
        self._run(set(ALL_PAIRS), _structural(), expected_ploidy_rows_per_sample=24)
        self.assertEqual(self.mock_struct.call_args.kwargs["expected_ploidy_rows_per_sample"], 24)

    def test_expected_ploidy_override_defaults_to_none(self):
        self._run(set(ALL_PAIRS), _structural())
        self.assertIsNone(self.mock_struct.call_args.kwargs["expected_ploidy_rows_per_sample"])


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
        # even though every listed pair is present in BigQuery -- the cross-family gap completeness
        # judges vacuously.
        r = self._run(set(ALL_PAIRS), _structural(cross_family_ok=False))
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["safe_to_delete_parquet"])
        self.assertFalse(r["structural_checks_ok"])
        self.assertFalse(r["cross_family_consistency_ok"])
        self.assertEqual(r["missing_files"], 0)


class TestVetScreensGateDeletionNotAllLoaded(VerifyAllLoadedTestBase):
    """A vet-screen flag is orthogonal to load completeness: all_loaded stays True (so the task
    succeeds), but the flag blocks Parquet deletion by default and is waived only under
    allow_flagged_vet_loads."""

    def test_duplication_flag_blocks_deletion_but_not_all_loaded(self):
        r = self._run(set(ALL_PAIRS), _structural(duplication_flagged=True))
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(r["vet_duplication_flagged"])
        self.assertFalse(r["safe_to_delete_parquet"])

    def test_duplication_flag_waived_allows_deletion(self):
        r = self._run(
            set(ALL_PAIRS),
            _structural(duplication_flagged=True, allow_flagged_vet_loads=True),
            allow_flagged_vet_loads=True,
        )
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["vet_duplication_flagged"])
        self.assertTrue(r["safe_to_delete_parquet"])

    def test_truncation_flag_blocks_deletion_but_not_all_loaded(self):
        r = self._run(set(ALL_PAIRS), _structural(truncation_flagged=True))
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(r["vet_truncation_flagged"])
        self.assertFalse(r["safe_to_delete_parquet"])

    def test_truncation_flag_waived_allows_deletion(self):
        r = self._run(
            set(ALL_PAIRS),
            _structural(truncation_flagged=True, allow_flagged_vet_loads=True),
            allow_flagged_vet_loads=True,
        )
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["vet_truncation_flagged"])
        self.assertTrue(r["safe_to_delete_parquet"])


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
    consistency only. The vet screens are deliberately excluded -- they gate safe_to_delete_parquet,
    not all_loaded."""

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

    def test_screen_flags_do_not_affect_exact_signal(self):
        self.assertTrue(self._ok(_structural(duplication_flagged=True)))
        self.assertTrue(self._ok(_structural(truncation_flagged=True)))


class TestComputeAllLoaded(unittest.TestCase):
    """The factual load-complete gate (what fail-loud aborts on): True only when nothing is
    missing/unmatched and the exact structural checks pass. The vet screens never enter here."""

    def test_true_when_everything_clean(self):
        self.assertTrue(verify_all_loaded.compute_all_loaded(set(), [], True))

    def test_missing_pair_blocks(self):
        self.assertFalse(verify_all_loaded.compute_all_loaded({("vet_001", 2)}, [], True))

    def test_unmatched_file_blocks(self):
        self.assertFalse(verify_all_loaded.compute_all_loaded(set(), ["gs://b/weird.parquet"], True))

    def test_structural_failure_blocks_even_when_files_all_present(self):
        self.assertFalse(verify_all_loaded.compute_all_loaded(set(), [], False))


class TestComputeSafeToDeleteParquet(unittest.TestCase):
    """The deletion gate proper: all_loaded AND no unwaived vet-screen flag."""

    def _safe(self, all_loaded, structural, allow):
        return verify_all_loaded.compute_safe_to_delete_parquet(all_loaded, structural, allow)

    def test_requires_all_loaded(self):
        # Even with no flags and the screens waived, an incomplete load is never safe to delete.
        self.assertFalse(self._safe(False, _structural(), False))
        self.assertFalse(self._safe(False, _structural(), True))

    def test_clean_load_is_safe(self):
        self.assertTrue(self._safe(True, _structural(), False))

    def test_duplication_flag_blocks_by_default(self):
        self.assertFalse(self._safe(True, _structural(duplication_flagged=True), False))

    def test_truncation_flag_blocks_by_default(self):
        self.assertFalse(self._safe(True, _structural(truncation_flagged=True), False))

    def test_flags_waived_when_allowed(self):
        self.assertTrue(self._safe(True, _structural(duplication_flagged=True), True))
        self.assertTrue(self._safe(True, _structural(truncation_flagged=True), True))
        self.assertTrue(self._safe(
            True, _structural(duplication_flagged=True, truncation_flagged=True), True))


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
    """The fail-loud operator message names every not-all-loaded cause. Only the exact checks appear;
    a vet-screen flag never fails all_loaded, so it must not surface here."""

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

    def test_screen_flags_do_not_appear(self):
        # Screens gate safe_to_delete_parquet, never all_loaded, so they must not appear as a
        # fail-loud reason even when flagged.
        self.assertEqual(self._reasons(vet_duplication_flagged=True, vet_truncation_flagged=True), [])


if __name__ == "__main__":
    unittest.main()
