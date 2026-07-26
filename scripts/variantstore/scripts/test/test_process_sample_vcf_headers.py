import os
import hashlib
import unittest

import process_sample_vcf_headers as psh

PROJECT = "test-project"
DATASET = "test_dataset"


def _norm(sql):
    """Collapse all runs of whitespace to single spaces so assertions are formatting-independent."""
    return " ".join(sql.split())


class TestHeaderLoadSql(unittest.TestCase):
    """Unit tests for the anti-join INSERT SQL builders (no BigQuery client required)."""

    def test_vcf_header_lines_is_anti_join_insert(self):
        sql = _norm(psh.vcf_header_lines_insert_sql(PROJECT, DATASET))
        # Targets the right table via INSERT (not MERGE).
        self.assertIn(f"INSERT INTO `{PROJECT}.{DATASET}.vcf_header_lines` "
                      "(vcf_header_lines_hash, vcf_header_lines, is_expected_unique)", sql)
        self.assertNotIn("MERGE", sql.upper())
        # Anti-join against the target on the hash key, keeping only rows not already present.
        self.assertIn(f"LEFT JOIN `{PROJECT}.{DATASET}.vcf_header_lines` t USING (vcf_header_lines_hash)", sql)
        self.assertIn("t.vcf_header_lines_hash IS NULL", sql)
        # Skips association-only (NULL text) scratch rows and dedups the source by hash.
        self.assertIn("s.vcf_header_lines IS NOT NULL", sql)
        self.assertIn("GROUP BY s.vcf_header_lines_hash", sql)

    def test_sample_vcf_header_is_anti_join_insert(self):
        sql = _norm(psh.sample_vcf_header_insert_sql(PROJECT, DATASET))
        self.assertIn(f"INSERT INTO `{PROJECT}.{DATASET}.sample_vcf_header` "
                      "(sample_id, vcf_header_lines_hash)", sql)
        self.assertNotIn("MERGE", sql.upper())
        self.assertIn("SELECT DISTINCT s.sample_id, s.vcf_header_lines_hash", sql)
        self.assertIn(f"LEFT JOIN `{PROJECT}.{DATASET}.sample_vcf_header` t "
                      "USING (sample_id, vcf_header_lines_hash)", sql)
        self.assertIn("t.sample_id IS NULL", sql)

    def test_clean_up_scratch_is_delete(self):
        sql = _norm(psh.clean_up_scratch_sql(PROJECT, DATASET))
        self.assertTrue(sql.startswith("DELETE FROM"), sql)
        self.assertIn(f"`{PROJECT}.{DATASET}.vcf_header_lines_scratch`", sql)


@unittest.skipUnless(
    os.environ.get("GVS_HEADER_IT_PROJECT") and os.environ.get("GVS_HEADER_IT_DATASET"),
    "Set GVS_HEADER_IT_PROJECT and GVS_HEADER_IT_DATASET (a throwaway dataset) to run the BigQuery "
    "integration test.",
)
class TestHeaderLoadIntegration(unittest.TestCase):
    """Real-BigQuery test of dedup + idempotency. Skipped unless the env vars above are set.

    Creates the three header tables in the given (throwaway) dataset, loads scratch with rows that
    share a hash across samples plus per-sample unique lines, runs the load, then loads and runs it a
    SECOND time to prove re-ingest inserts nothing new.
    """

    @classmethod
    def setUpClass(cls):
        from google.cloud import bigquery
        cls.bigquery = bigquery
        cls.project = os.environ["GVS_HEADER_IT_PROJECT"]
        cls.dataset = os.environ["GVS_HEADER_IT_DATASET"]
        cls.client = bigquery.Client(project=cls.project)
        cls.fq = lambda _cls, t: f"{cls.project}.{cls.dataset}.{t}"

    def _run(self, sql):
        self.client.query(sql).result()

    def _count(self, table):
        rows = list(self.client.query(f"SELECT COUNT(*) AS c FROM `{self.fq(table)}`").result())
        return rows[0].c

    def _create_tables(self):
        self._run(f"CREATE OR REPLACE TABLE `{self.fq('vcf_header_lines_scratch')}` "
                  "(sample_id INT64, vcf_header_lines STRING, vcf_header_lines_hash STRING, is_expected_unique BOOL)")
        self._run(f"CREATE OR REPLACE TABLE `{self.fq('vcf_header_lines')}` "
                  "(vcf_header_lines_hash STRING, vcf_header_lines STRING, is_expected_unique BOOL)")
        self._run(f"CREATE OR REPLACE TABLE `{self.fq('sample_vcf_header')}` "
                  "(sample_id INT64, vcf_header_lines_hash STRING)")

    def _load_scratch(self):
        # Populate via DML INSERT (not streaming) so the later DELETE isn't blocked by a streaming buffer.
        def md5(s):
            return hashlib.md5(s.encode("utf-8")).hexdigest()
        blob, cmd1, cmd2 = "BLOB", "CMD1", "CMD2"
        rows = [
            (1, blob, md5(blob), "false"),   # shared blob, sample 1
            (2, blob, md5(blob), "false"),   # same hash, sample 2 -> dedups to one vcf_header_lines row
            (1, cmd1, md5(cmd1), "true"),
            (2, cmd2, md5(cmd2), "true"),
        ]
        values = ", ".join(f"({sid}, '{txt}', '{h}', {u})" for (sid, txt, h, u) in rows)
        self._run(f"INSERT INTO `{self.fq('vcf_header_lines_scratch')}` "
                  f"(sample_id, vcf_header_lines, vcf_header_lines_hash, is_expected_unique) VALUES {values}")

    def test_dedup_and_idempotency(self):
        self._create_tables()

        # First ingest.
        self._load_scratch()
        psh.process_sample_vcf_headers(self.project, self.dataset, client=self.client)
        self.assertEqual(self._count("vcf_header_lines"), 3, "one row per distinct hash (BLOB, CMD1, CMD2)")
        self.assertEqual(self._count("sample_vcf_header"), 4, "distinct (sample_id, hash) pairs")
        self.assertEqual(self._count("vcf_header_lines_scratch"), 0, "scratch cleaned up after load")

        # Second ingest of the same batch: anti-join INSERTs must add nothing.
        self._load_scratch()
        psh.process_sample_vcf_headers(self.project, self.dataset, client=self.client)
        self.assertEqual(self._count("vcf_header_lines"), 3, "re-ingest must not duplicate header lines")
        self.assertEqual(self._count("sample_vcf_header"), 4, "re-ingest must not duplicate associations")

    @classmethod
    def tearDownClass(cls):
        for t in ("vcf_header_lines_scratch", "vcf_header_lines", "sample_vcf_header"):
            cls.client.query(f"DROP TABLE IF EXISTS `{cls.project}.{cls.dataset}.{t}`").result()


if __name__ == "__main__":
    unittest.main()
