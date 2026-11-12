import unittest
from collections import namedtuple
from unittest.mock import patch

import process_vrs_alleles as pva

PROJECT = "test-project"
DATASET = "test_dataset"
SAMPLES = [12, 3, 7]

SampleRow = namedtuple("SampleRow", ["sample_id"])
IdRow = namedtuple("IdRow", ["vrs_allele_id"])
PairRow = namedtuple("PairRow", ["sample_id", "vrs_allele_id"])


def _norm(sql):
    """Collapse all runs of whitespace to single spaces so assertions are formatting-independent."""
    return " ".join(sql.split())


class TestVrsAlleleSql(unittest.TestCase):
    """Intent-level guards on the SQL builders (no BigQuery client required)."""

    def test_sample_ids_render_as_sorted_literals(self):
        # Literals, so BigQuery can prune the sample_id-partitioned vet tables.
        sql = _norm(pva.missing_vet_ids_sql(PROJECT, DATASET, SAMPLES))
        self.assertIn("v.sample_id IN (3, 7, 12)", sql)

    def test_empty_sample_ids_rejected(self):
        # An empty IN () is a syntax error; worse, a builder that dropped the filter would touch every sample.
        with self.assertRaises(ValueError):
            pva.clean_up_scratch_sql(PROJECT, DATASET, [])

    def test_insert_targets_vrs_allele_and_skips_existing_ids(self):
        sql = _norm(pva.vrs_allele_insert_sql(PROJECT, DATASET, SAMPLES))
        self.assertIn(f"INSERT INTO `{PROJECT}.{DATASET}.vrs_allele`", sql)
        self.assertIn("t.vrs_allele_id IS NULL", sql)
        self.assertIn("GROUP BY s.vrs_allele_id", sql)
        self.assertIn("s.sample_id IN (3, 7, 12)", sql)

    def test_insert_writes_every_vrs_allele_column(self):
        sql = _norm(pva.vrs_allele_insert_sql(PROJECT, DATASET, SAMPLES))
        for column in pva.VRS_ALLELE_COLUMNS:
            self.assertIn(column, sql)
        self.assertNotIn("sample_id)", sql.split("SELECT")[0])  # not in the INSERT column list

    def test_conflict_check_compares_scratch_with_existing_rows(self):
        sql = _norm(pva.conflicting_ids_sql(PROJECT, DATASET, SAMPLES))
        self.assertIn(f"`{PROJECT}.{DATASET}.vrs_allele_scratch`", sql)
        self.assertIn(f"`{PROJECT}.{DATASET}.vrs_allele`", sql)
        self.assertIn("COUNT(DISTINCT TO_JSON_STRING(STRUCT(", sql)

    def test_missing_check_reads_numbered_vet_tables_only(self):
        sql = _norm(pva.missing_vet_ids_sql(PROJECT, DATASET, SAMPLES))
        self.assertIn(f"`{PROJECT}.{DATASET}.vet_*`", sql)
        self.assertIn("REGEXP_CONTAINS(_TABLE_SUFFIX, r'^[0-9]+$')", sql)
        self.assertIn("a.vrs_allele_id IS NULL", sql)

    def test_missing_check_ignores_spanning_deletion_placeholder(self):
        # vet.vrs_allele_ids holds '' for a `*` ALT; it has no vrs_allele row and must not be reported.
        sql = _norm(pva.missing_vet_ids_sql(PROJECT, DATASET, SAMPLES))
        self.assertIn("vet_id != ''", sql)

    def test_clean_up_is_scoped_delete(self):
        sql = _norm(pva.clean_up_scratch_sql(PROJECT, DATASET, SAMPLES))
        self.assertTrue(sql.startswith(f"DELETE FROM `{PROJECT}.{DATASET}.vrs_allele_scratch`"), sql)
        self.assertIn("WHERE sample_id IN (3, 7, 12)", sql)


class FakeBigQuery:
    """Stands in for _execute: answers each query by its label and records the order they ran in."""

    def __init__(self, samples=SAMPLES, conflicting=(), missing=(), duplicates=()):
        self.answers = {
            "vrs_allele_scratch samples": [SampleRow(s) for s in samples],
            "vrs_allele conflicting ids": [IdRow(i) for i in conflicting],
            "vet vrs_allele_ids missing from vrs_allele": [PairRow(s, i) for s, i in missing],
            "vrs_allele duplicate ids": [IdRow(i) for i in duplicates],
        }
        self.labels = []

    def __call__(self, client, label, sql):
        self.labels.append(label)
        return {"results": self.answers.get(label, [])}


class TestProcessVrsAlleles(unittest.TestCase):
    """The order of operations: every check runs before the scratch rows are deleted, and a failed
    check leaves them in place."""

    def run_with(self, fake):
        with patch.object(pva, "_execute", fake):
            pva.process_vrs_alleles(PROJECT, DATASET, client=object())

    def test_clean_run_merges_checks_then_deletes(self):
        fake = FakeBigQuery()
        self.run_with(fake)
        self.assertEqual(fake.labels, [
            "vrs_allele_scratch samples",
            "vrs_allele conflicting ids",
            "vrs_allele",
            "vet vrs_allele_ids missing from vrs_allele",
            "vrs_allele duplicate ids",
            "vrs_allele_scratch",
        ])

    def test_empty_scratch_does_nothing(self):
        fake = FakeBigQuery(samples=[])
        self.run_with(fake)
        self.assertEqual(fake.labels, ["vrs_allele_scratch samples"])

    def test_conflicting_ids_abort_before_merge(self):
        fake = FakeBigQuery(conflicting=["ga4gh:VA.clash"])
        with self.assertRaisesRegex(RuntimeError, "ga4gh:VA.clash"):
            self.run_with(fake)
        self.assertNotIn("vrs_allele", fake.labels)
        self.assertNotIn("vrs_allele_scratch", fake.labels)

    def test_missing_vet_ids_abort_before_delete(self):
        fake = FakeBigQuery(missing=[(7, "ga4gh:VA.dangling")])
        with self.assertRaisesRegex(RuntimeError, "ga4gh:VA.dangling"):
            self.run_with(fake)
        self.assertNotIn("vrs_allele_scratch", fake.labels)

    def test_duplicate_ids_abort_before_delete(self):
        fake = FakeBigQuery(duplicates=["ga4gh:VA.twice"])
        with self.assertRaisesRegex(RuntimeError, "ga4gh:VA.twice"):
            self.run_with(fake)
        self.assertNotIn("vrs_allele_scratch", fake.labels)


if __name__ == "__main__":
    unittest.main()
