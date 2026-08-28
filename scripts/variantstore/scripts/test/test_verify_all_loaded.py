"""
Unit tests for verify_all_loaded.py, focused on the VS-1989 composition change: the independent
structural checks must be able to block Parquet deletion (all_loaded=False) even when the
shared-predicate presence check is fully satisfied, and the results JSON must carry both the new
structural booleans and the pre-existing keys the WDL already reads.

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


def _structural(completeness_ok=True, cardinality_ok=True, duplication_flagged=False,
                strict_vet_screen=False):
    """A run_structural_checks return value with the keys verify_all_loaded consumes."""
    return {
        "completeness_ok": completeness_ok,
        "cardinality_ok": cardinality_ok,
        "duplication_flagged": duplication_flagged,
        "strict_vet_screen": strict_vet_screen,
        "details": {
            "family_completeness": {"ok": completeness_ok, "per_family": {}},
            "cardinality": {},
            "duplication_screen": {},
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

    def _run(self, loaded_pairs, structural, strict_vet_screen=False):
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
                strict_vet_screen=strict_vet_screen,
            )

    def _written_json(self):
        with open(os.path.join(self.out_dir, "verification_results.json")) as f:
            return json.load(f)


class TestHappyPath(VerifyAllLoadedTestBase):
    def test_all_loaded_and_json_keys(self):
        r = self._run(set(ALL_PAIRS), _structural())

        self.assertTrue(r["all_loaded"])
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
        self.assertFalse(r["vet_duplication_flagged"])
        self.assertIn("structural_checks", r)
        # The JSON on disk matches the returned dict.
        self.assertEqual(self._written_json(), r)

    def test_expected_by_family_passed_to_structural_checks(self):
        self._run(set(ALL_PAIRS), _structural())
        expected_by_family = self.mock_struct.call_args[0][2]
        self.assertEqual(set(expected_by_family["vet"]), {1, 2})
        self.assertEqual(set(expected_by_family["ref_ranges"]), {1, 2})
        self.assertEqual(set(expected_by_family["sample_chromosome_ploidy"]), {1, 2})


class TestStructuralChecksGate(VerifyAllLoadedTestBase):
    def test_cardinality_failure_blocks_deletion(self):
        # Every pair is present to the shared predicate, but a partial load fails cardinality.
        r = self._run(set(ALL_PAIRS), _structural(cardinality_ok=False))
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["structural_checks_ok"])
        self.assertFalse(r["ploidy_cardinality_ok"])
        self.assertEqual(r["missing_files"], 0)

    def test_completeness_failure_blocks_deletion(self):
        r = self._run(set(ALL_PAIRS), _structural(completeness_ok=False))
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["family_completeness_ok"])

    def test_vet_duplication_warns_by_default(self):
        r = self._run(set(ALL_PAIRS), _structural(duplication_flagged=True))
        self.assertTrue(r["all_loaded"])
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(r["vet_duplication_flagged"])

    def test_vet_duplication_gates_under_strict(self):
        r = self._run(
            set(ALL_PAIRS),
            _structural(duplication_flagged=True, strict_vet_screen=True),
            strict_vet_screen=True,
        )
        self.assertFalse(r["all_loaded"])
        self.assertFalse(r["structural_checks_ok"])
        self.assertTrue(r["vet_duplication_flagged"])


class TestSharedPredicateStillGates(VerifyAllLoadedTestBase):
    def test_missing_pair_blocks_even_when_structural_ok(self):
        loaded = set(ALL_PAIRS) - {("vet_001", 2)}
        r = self._run(loaded, _structural())
        self.assertFalse(r["all_loaded"])
        self.assertEqual(r["missing_files"], 1)
        self.assertTrue(r["structural_checks_ok"])
        self.assertTrue(os.path.exists(r["missing_files_list"]))


if __name__ == "__main__":
    unittest.main()
