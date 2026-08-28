"""
Unit tests for verify_structural_checks.py -- the independent, row-count-based checks added for
VS-1989.

The BigQuery-reading functions are tested by patching verify_structural_checks.bigquery and asserting
on the generated SQL (mirroring the style in test_parquet_loading.py). The assessment functions are
pure and are tested directly. run_structural_checks is tested with the two query helpers patched so no
BigQuery is touched.
"""

import os
import sys
import unittest
from unittest.mock import MagicMock, patch

# Add parent directory to path for imports (matches test_parquet_loading.py).
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import verify_structural_checks as vsc
from verify_structural_checks import (
    family_for_table,
    get_partition_row_counts,
    get_ploidy_row_counts,
    assess_family_completeness,
    assess_cardinality,
    assess_duplication_screen,
    run_structural_checks,
)


def _row(**kwargs):
    """A stand-in BigQuery Row: attribute access returns the given values."""
    row = MagicMock()
    for k, val in kwargs.items():
        setattr(row, k, val)
    return row


class TestFamilyForTable(unittest.TestCase):
    def test_matches_superpartitioned_names(self):
        self.assertEqual(family_for_table("vet_001", ["vet", "ref_ranges"]), "vet")
        self.assertEqual(family_for_table("ref_ranges_042", ["vet", "ref_ranges"]), "ref_ranges")

    def test_requires_digits_after_prefix(self):
        self.assertIsNone(family_for_table("vet_x", ["vet"]))
        self.assertIsNone(family_for_table("vet_", ["vet"]))

    def test_does_not_match_longer_name_with_shared_start(self):
        # "vet" must not swallow a differently-named table that merely starts with the letters.
        self.assertIsNone(family_for_table("vetting_1", ["vet"]))

    def test_returns_none_for_unknown(self):
        self.assertIsNone(family_for_table("sample_chromosome_ploidy", ["vet", "ref_ranges"]))


class TestGetPartitionRowCounts(unittest.TestCase):
    @patch("verify_structural_checks.bigquery")
    def test_query_reads_total_rows_and_not_bytes(self, mock_bq):
        mock_client = MagicMock()
        mock_bq.Client.return_value = mock_client
        mock_client.query.return_value = []

        get_partition_row_counts("proj", "ds", ["vet", "ref_ranges"])

        query = mock_client.query.call_args[0][0]
        self.assertIn("INFORMATION_SCHEMA.PARTITIONS", query)
        self.assertIn("total_rows", query)
        self.assertIn("^vet_[0-9]+$", query)
        self.assertIn("^ref_ranges_[0-9]+$", query)
        # Independence from the loader predicate: it must NOT gate on bytes.
        self.assertNotIn("total_logical_bytes", query)

    @patch("verify_structural_checks.bigquery")
    def test_returns_tuples(self, mock_bq):
        mock_client = MagicMock()
        mock_bq.Client.return_value = mock_client
        mock_client.query.return_value = [
            _row(table_name="vet_001", sample_id=1, total_rows=100),
            _row(table_name="ref_ranges_001", sample_id=1, total_rows=50),
        ]

        result = get_partition_row_counts("proj", "ds", ["vet", "ref_ranges"])

        self.assertEqual(result, [("vet_001", 1, 100), ("ref_ranges_001", 1, 50)])

    @patch("verify_structural_checks.bigquery")
    def test_empty_prefixes_skips_query(self, mock_bq):
        mock_client = MagicMock()
        mock_bq.Client.return_value = mock_client

        self.assertEqual(get_partition_row_counts("proj", "ds", []), [])
        mock_client.query.assert_not_called()

    def test_invalid_prefix_raises(self):
        with self.assertRaises(ValueError):
            get_partition_row_counts("proj", "ds", ["vet; DROP TABLE x--"])


class TestGetPloidyRowCounts(unittest.TestCase):
    @patch("verify_structural_checks.bigquery")
    def test_query_groups_by_sample(self, mock_bq):
        mock_client = MagicMock()
        mock_bq.Client.return_value = mock_client
        mock_client.query.return_value = [_row(sample_id=1, n=24), _row(sample_id=2, n=24)]

        result = get_ploidy_row_counts("proj", "ds", "sample_chromosome_ploidy")

        query = mock_client.query.call_args[0][0]
        self.assertIn("sample_chromosome_ploidy", query)
        self.assertIn("GROUP BY sample_id", query)
        self.assertIn("COUNT(*)", query)
        self.assertEqual(result, {1: 24, 2: 24})

    def test_invalid_table_raises(self):
        with self.assertRaises(ValueError):
            get_ploidy_row_counts("proj", "ds", "bad table")


class TestAssessCardinality(unittest.TestCase):
    def test_uniform_counts_pass(self):
        r = assess_cardinality({i: 24 for i in range(1, 6)}, set(range(1, 6)))
        self.assertTrue(r["ok"])
        self.assertEqual(r["mode"], 24)
        self.assertEqual(r["distinct_samples"], 5)
        self.assertEqual(r["deviating_samples"], [])

    def test_does_not_assume_24(self):
        # A callset whose modal count is 25 (e.g. chrM ingested) still passes when uniform.
        r = assess_cardinality({i: 25 for i in range(1, 4)}, set(range(1, 4)))
        self.assertTrue(r["ok"])
        self.assertEqual(r["mode"], 25)

    def test_missing_sample_fails(self):
        r = assess_cardinality({1: 24, 2: 24}, {1, 2, 3})
        self.assertFalse(r["ok"])
        self.assertEqual(r["missing_samples"], [3])

    def test_partial_and_duplicate_flagged(self):
        # sample 3 short (partial load), sample 4 doubled (duplication).
        r = assess_cardinality({1: 24, 2: 24, 3: 20, 4: 48, 5: 24}, set(range(1, 6)))
        self.assertFalse(r["ok"])
        self.assertEqual({d["sample_id"] for d in r["deviating_samples"]}, {3, 4})
        self.assertEqual(r["min"], 20)
        self.assertEqual(r["max"], 48)

    def test_extra_samples_in_table_ignored(self):
        # A sample already in the table from a prior load (99) is not expected this run: ignored.
        r = assess_cardinality({1: 24, 2: 24, 99: 24}, {1, 2})
        self.assertTrue(r["ok"])
        self.assertEqual(r["distinct_samples"], 2)

    def test_zero_count_treated_as_missing(self):
        r = assess_cardinality({1: 24, 2: 0}, {1, 2})
        self.assertFalse(r["ok"])
        self.assertEqual(r["missing_samples"], [2])

    def test_no_expected_is_ok(self):
        r = assess_cardinality({}, set())
        self.assertTrue(r["ok"])
        self.assertIsNone(r["mode"])


class TestAssessFamilyCompleteness(unittest.TestCase):
    def _run(self, part, reg, exp):
        return assess_family_completeness(
            part, reg, exp, ["vet", "ref_ranges"], ["sample_chromosome_ploidy"]
        )

    def test_all_present_passes(self):
        part = [("vet_001", 1, 100), ("vet_001", 2, 100),
                ("ref_ranges_001", 1, 50), ("ref_ranges_001", 2, 50)]
        reg = {"sample_chromosome_ploidy": {1: 24, 2: 24}}
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = self._run(part, reg, exp)
        self.assertTrue(r["ok"])
        self.assertEqual(r["per_family"]["vet"]["present"], 2)

    def test_empty_partition_is_partial_load(self):
        # sample 2 has a vet partition with 0 rows -> present to a bytes-based check, empty in fact.
        part = [("vet_001", 1, 100), ("vet_001", 2, 0),
                ("ref_ranges_001", 1, 50), ("ref_ranges_001", 2, 50)]
        reg = {"sample_chromosome_ploidy": {1: 24, 2: 24}}
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = self._run(part, reg, exp)
        self.assertFalse(r["ok"])
        self.assertEqual(r["per_family"]["vet"]["empty_partition_samples"], [2])

    def test_missing_sample_in_one_family(self):
        part = [("vet_001", 1, 100), ("ref_ranges_001", 1, 50)]  # sample 2 absent everywhere
        reg = {"sample_chromosome_ploidy": {1: 24}}
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = self._run(part, reg, exp)
        self.assertFalse(r["ok"])
        self.assertEqual(r["per_family"]["vet"]["missing_samples"], [2])
        self.assertEqual(r["per_family"]["sample_chromosome_ploidy"]["missing_samples"], [2])

    def test_empty_partition_for_unexpected_sample_ignored(self):
        # A 0-row partition for a sample not expected this run must not fail completeness.
        part = [("vet_001", 1, 100), ("vet_001", 9, 0), ("ref_ranges_001", 1, 50)]
        reg = {"sample_chromosome_ploidy": {1: 24}}
        exp = {"vet": {1}, "ref_ranges": {1}, "sample_chromosome_ploidy": {1}}
        r = self._run(part, reg, exp)
        self.assertTrue(r["ok"])


class TestAssessDuplicationScreen(unittest.TestCase):
    def test_flags_high_outlier(self):
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 220)]
        r = assess_duplication_screen(part, "vet", set(range(1, 10)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["median"], 100)
        self.assertEqual([o["sample_id"] for o in r["outliers"]], [9])
        self.assertAlmostEqual(r["outliers"][0]["ratio"], 2.2)

    def test_uniform_no_outliers(self):
        part = [("vet_001", i, 100) for i in range(1, 6)]
        r = assess_duplication_screen(part, "vet", set(range(1, 6)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["outliers"], [])

    def test_other_family_not_screened(self):
        part = [("vet_001", i, 100) for i in range(1, 6)]
        r = assess_duplication_screen(part, "ref_ranges", set(range(1, 6)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["samples_screened"], 0)
        self.assertEqual(r["outliers"], [])

    def test_only_expected_samples_screened(self):
        part = [("vet_001", i, 100) for i in range(1, 6)] + [("vet_001", 99, 1000)]
        r = assess_duplication_screen(part, "vet", set(range(1, 6)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["samples_screened"], 5)
        self.assertEqual(r["outliers"], [])


class TestRunStructuralChecks(unittest.TestCase):
    """run_structural_checks with the two BigQuery helpers patched (no BQ touched)."""

    def _patch(self, partition_rows, ploidy_counts):
        p1 = patch("verify_structural_checks.get_partition_row_counts", return_value=partition_rows)
        p2 = patch("verify_structural_checks.get_ploidy_row_counts", return_value=ploidy_counts)
        p1.start()
        p2.start()
        self.addCleanup(p1.stop)
        self.addCleanup(p2.stop)

    def test_all_good(self):
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: 24, 2: 24})
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = run_structural_checks("proj", "ds", exp)
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])
        self.assertFalse(r["duplication_flagged"])
        self.assertIn("ref_ranges", r["details"]["duplication_unscreened"]["families"])
        self.assertNotIn("vet", r["details"]["duplication_unscreened"]["families"])
        self.assertIn("backfill_caveat", r["details"]["cardinality"]["sample_chromosome_ploidy"])

    def test_partial_ploidy_load_fails_cardinality(self):
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: 24, 2: 20})  # sample 2 short
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = run_structural_checks("proj", "ds", exp)
        self.assertTrue(r["completeness_ok"])
        self.assertFalse(r["cardinality_ok"])

    def test_vet_duplication_flag_warns_not_gates_by_default(self):
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 300)]
        part += [("ref_ranges_001", i, 50) for i in range(1, 10)]
        self._patch(part, {i: 24 for i in range(1, 10)})
        exp = {"vet": set(range(1, 10)), "ref_ranges": set(range(1, 10)),
               "sample_chromosome_ploidy": set(range(1, 10))}

        r = run_structural_checks("proj", "ds", exp, strict_vet_screen=False)
        self.assertTrue(r["duplication_flagged"])
        # Completeness and cardinality (the hard gates) are unaffected by a warn-only screen.
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])
        self.assertFalse(r["strict_vet_screen"])


if __name__ == "__main__":
    unittest.main()
