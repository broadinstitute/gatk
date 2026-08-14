import hashlib
import os
import unittest

import check_vcf_headers as cvh

PROJECT = "test-project"
DATASET = "test_dataset"


def _norm(sql):
    """Collapse all runs of whitespace to single spaces so assertions are formatting-independent."""
    return " ".join(sql.split())


class TestSqlBuilders(unittest.TestCase):
    """Intent-level guards on the SQL builders (no BigQuery client required).

    These assert the properties that must hold for ANY correct rewrite -- the right tables, the
    non-control/non-withdrawn cohort filter, and the predicates that pick out the reblocking and
    DRAGEN signals. They do NOT pin exact join text; behavior is covered by TestEvaluate and the
    (opt-in) integration test below.
    """

    def test_per_sample_summary_targets_and_filters(self):
        sql = _norm(cvh.per_sample_summary_sql(PROJECT, DATASET))
        self.assertIn(f"`{PROJECT}.{DATASET}.sample_info`", sql)
        self.assertIn(f"`{PROJECT}.{DATASET}.sample_vcf_header`", sql)
        self.assertIn(f"`{PROJECT}.{DATASET}.vcf_header_lines`", sql)
        # loaded, non-control cohort
        self.assertIn("is_control = FALSE", sql)
        self.assertIn("withdrawn IS NULL", sql)
        # reblocking signal lives in an is_expected_unique chunk
        self.assertIn("ReblockGVCF", sql)
        self.assertIn("is_expected_unique", sql)

    def test_dragen_breakdown_filters_on_dragen_command_line(self):
        sql = _norm(cvh.dragen_version_breakdown_sql(PROJECT, DATASET))
        # must isolate the ID=dragen command line (not the ID=HashTableBuild DRAGENCommandLine)
        self.assertIn("DRAGENCommandLine=<ID=dragen,", sql)
        self.assertIn("is_expected_unique = TRUE", sql)
        self.assertIn(cvh.SW_VERSION_REGEX, sql)
        self.assertIn("is_control = FALSE", sql)
        self.assertIn("withdrawn IS NULL", sql)

    def test_orphan_hash_sql(self):
        sql = _norm(cvh.orphan_hash_sql(PROJECT, DATASET))
        self.assertIn(f"`{PROJECT}.{DATASET}.sample_vcf_header`", sql)
        self.assertIn("LEFT JOIN", sql)
        self.assertIn("IS NULL", sql)

    def test_shared_blob_distribution_is_false_blob(self):
        sql = _norm(cvh.shared_blob_distribution_sql(PROJECT, DATASET))
        self.assertIn("is_expected_unique = FALSE", sql)


class TestTripletParsing(unittest.TestCase):
    def test_final_triplet_extracted(self):
        self.assertEqual(cvh.dragen_version_triplet("SW: 05.021.604.3.7.8"), "3.7.8")
        self.assertEqual(cvh.dragen_version_triplet("SW: 07.021.604.3.7.8"), "3.7.8")
        self.assertEqual(cvh.dragen_version_triplet("SW: 05.021.604.3.4.12"), "3.4.12")

    def test_none_and_too_short(self):
        self.assertIsNone(cvh.dragen_version_triplet(None))
        self.assertIsNone(cvh.dragen_version_triplet("SW: 3.7"))
        self.assertIsNone(cvh.dragen_version_triplet("no digits here"))


class TestExpectedDragenSpec(unittest.TestCase):
    """Parsing of the overloaded --expected_dragen_version (exact value vs LOW-HIGH range)."""

    def test_none_or_empty(self):
        self.assertEqual(cvh.parse_expected_dragen_spec(None), (None, None))
        self.assertEqual(cvh.parse_expected_dragen_spec(""), (None, None))

    def test_exact(self):
        self.assertEqual(cvh.parse_expected_dragen_spec("3.7.8"), ('exact', '3.7.8'))

    def test_range_bare_is_inclusive_both_ends(self):
        kind, payload = cvh.parse_expected_dragen_spec("3.4.12-3.7.8")
        self.assertEqual(kind, 'range')
        self.assertEqual(payload, cvh.VersionRange((3, 4, 12), True, (3, 7, 8), True))

    def test_range_interval_notation_sets_inclusivity(self):
        self.assertEqual(cvh.parse_expected_dragen_spec("[3.7.8-3.8)")[1],
                         cvh.VersionRange((3, 7, 8), True, (3, 8), False))
        self.assertEqual(cvh.parse_expected_dragen_spec("(3.7-3.8)")[1],
                         cvh.VersionRange((3, 7), False, (3, 8), False))
        self.assertEqual(cvh.parse_expected_dragen_spec("(3.4.12-3.7.8]")[1],
                         cvh.VersionRange((3, 4, 12), False, (3, 7, 8), True))

    def test_malformed_range_raises(self):
        for bad in ("3.4.12-", "-3.7.8", "3.4.12-3.5.0-3.7.8", "[3.7.8]"):
            with self.assertRaises(ValueError):
                cvh.parse_expected_dragen_spec(bad)

    def test_low_greater_than_high_raises(self):
        with self.assertRaises(ValueError):
            cvh.parse_expected_dragen_spec("3.7.8-3.4.12")

    def test_empty_interval_raises(self):
        # low == high but an endpoint is exclusive -> empty interval.
        with self.assertRaises(ValueError):
            cvh.parse_expected_dragen_spec("[3.7.8-3.7.8)")

    def test_whitespace_only_is_none(self):
        # A blank/whitespace-only value is treated as unspecified, not a never-matching exact value.
        self.assertEqual(cvh.parse_expected_dragen_spec("   "), (None, None))

    def test_exact_rejects_malformed(self):
        # Typos / stray characters in an exact value must fail loudly, not silently never match.
        for bad in ("3.7,8", "3.7.8x", "v3.7.8", "3.7 .8", "3..7", "abc"):
            with self.assertRaises(ValueError):
                cvh.parse_expected_dragen_spec(bad)

    def test_range_rejects_malformed_bound(self):
        for bad in ("3.x-3.7.8", "3.7.8-3.8y", "3.7,8-3.8", "[3.a-3.8)"):
            with self.assertRaises(ValueError):
                cvh.parse_expected_dragen_spec(bad)


class TestIdentifierValidation(unittest.TestCase):
    """project_id / dataset_name are interpolated into SQL, so they must be validated first."""

    def test_valid_identifiers_pass(self):
        # No exception for typical GVS identifiers.
        cvh._validate_bq_identifiers("broad-dsde-methods", "ng_trial_ingest_1")
        cvh._validate_bq_identifiers("gvs.internal", "dataset123")

    def test_bad_project_raises(self):
        for bad in ("proj; DROP TABLE x", "proj`", "bad project", "-leading-hyphen", ""):
            with self.assertRaises(ValueError):
                cvh._validate_bq_identifiers(bad, "good_dataset")

    def test_bad_dataset_raises(self):
        # BigQuery dataset names disallow hyphens, spaces, and punctuation.
        for bad in ("ds-with-dash", "ds;DROP", "bad dataset", "ds`", ""):
            with self.assertRaises(ValueError):
                cvh._validate_bq_identifiers("good-project", bad)


def _summary(**overrides):
    base = {
        'expected_samples': 3, 'samples_with_headers': 3, 'samples_missing_headers': 0,
        'samples_missing_unique_chunk': 0, 'samples_not_reblocked': 0,
        'example_missing_headers': [], 'example_missing_unique_chunk': [], 'example_not_reblocked': [],
    }
    base.update(overrides)
    return base


_NO_ORPHANS = {'orphan_associations': 0, 'affected_samples': 0}
_ONE_BLOB = [{'blob_hash': 'abc', 'n_samples': 3}]
_DRAGEN_378 = [{'sw_version': 'SW: 05.021.604.3.7.8', 'n_samples': 2},
               {'sw_version': 'SW: 07.021.604.3.7.8', 'n_samples': 1}]


class TestEvaluate(unittest.TestCase):
    """The pass/fail decision logic -- the heart of the tool -- with no BigQuery dependency."""

    def _names(self, checks):
        return {c.name: c for c in checks}

    def test_healthy_cohort_passes(self):
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='3.7.8')
        self.assertTrue(ok, [(c.name, c.passed) for c in checks])

    def test_consistency_only_passes_without_expected(self):
        ok, _ = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB)
        self.assertTrue(ok)

    def test_missing_headers_fails(self):
        ok, checks = cvh.evaluate(
            _summary(samples_missing_headers=1, samples_with_headers=2, example_missing_headers=['S1']),
            _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB, expected_dragen_version='3.7.8')
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['all_samples_have_headers'].passed)

    def test_not_reblocked_fails(self):
        ok, checks = cvh.evaluate(
            _summary(samples_not_reblocked=2, example_not_reblocked=['S2', 'S3']),
            _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB, expected_dragen_version='3.7.8')
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['reblocking'].passed)

    def test_not_reblocked_informational_when_allowed(self):
        # require_reblocking=False downgrades the reblocking check to informational: it still reports
        # the non-reblocked samples but no longer fails the overall run.
        ok, checks = cvh.evaluate(
            _summary(samples_not_reblocked=2, example_not_reblocked=['S2', 'S3']),
            _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB, expected_dragen_version='3.7.8',
            require_reblocking=False)
        self.assertTrue(ok)
        reblock = self._names(checks)['reblocking']
        self.assertFalse(reblock.fatal)
        self.assertFalse(reblock.passed)  # still reports the truth, just not fatal

    def test_no_samples_fails(self):
        ok, checks = cvh.evaluate(_summary(expected_samples=0, samples_with_headers=0),
                                  _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB)
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['samples_present'].passed)

    def test_mixed_dragen_triplets_fails(self):
        mixed = [{'sw_version': 'SW: 05.021.604.3.7.8', 'n_samples': 2},
                 {'sw_version': 'SW: 05.021.604.3.4.12', 'n_samples': 1}]
        ok, checks = cvh.evaluate(_summary(), mixed, _NO_ORPHANS, _ONE_BLOB)
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['dragen_version'].passed)

    def test_expected_version_mismatch_fails(self):
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='3.4.12')
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['dragen_version'].passed)

    def test_range_in_bounds_passes(self):
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='3.4.12-3.7.8')
        self.assertTrue(ok, [(c.name, c.passed) for c in checks])
        self.assertTrue(self._names(checks)['dragen_version'].passed)

    def test_range_allows_multiple_triplets_within_bounds(self):
        # A range deliberately relaxes the single-triplet requirement: 3.4.12 and 3.7.8 both pass.
        mixed = [{'sw_version': 'SW: 05.021.604.3.7.8', 'n_samples': 2},
                 {'sw_version': 'SW: 05.021.604.3.4.12', 'n_samples': 1}]
        ok, checks = cvh.evaluate(_summary(), mixed, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='3.4.0-3.7.8')
        self.assertTrue(ok, [(c.name, c.passed) for c in checks])

    def test_range_out_of_bounds_fails(self):
        out = [{'sw_version': 'SW: 05.021.604.3.4.0', 'n_samples': 1}]
        ok, checks = cvh.evaluate(_summary(), out, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='3.5.0-3.7.8')
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['dragen_version'].passed)

    def test_malformed_range_fails_fatally(self):
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='3.7.8-')
        self.assertFalse(ok)
        dragen = self._names(checks)['dragen_version']
        self.assertTrue(dragen.fatal)
        self.assertFalse(dragen.passed)

    def test_range_exclusive_high_excludes_boundary(self):
        # _DRAGEN_378 is all 3.7.8; '[3.7.7-3.7.8)' excludes 3.7.8 -> fail.
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='[3.7.7-3.7.8)')
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['dragen_version'].passed)

    def test_range_exclusive_high_includes_below_boundary(self):
        # 3.7.8 < 3.8, so '[3.7.8-3.8)' includes it -> pass.
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='[3.7.8-3.8)')
        self.assertTrue(ok, [(c.name, c.passed) for c in checks])

    def test_no_dragen_lines_is_informational_without_expected(self):
        ok, checks = cvh.evaluate(_summary(), [], _NO_ORPHANS, _ONE_BLOB)
        self.assertTrue(ok)
        dragen = self._names(checks)['dragen_version']
        self.assertFalse(dragen.fatal)  # informational: nothing to check for a non-DRAGEN cohort

    def test_no_dragen_lines_with_expected_fails(self):
        ok, checks = cvh.evaluate(_summary(), [], _NO_ORPHANS, _ONE_BLOB,
                                  expected_dragen_version='3.7.8')
        self.assertFalse(ok)
        dragen = self._names(checks)['dragen_version']
        self.assertTrue(dragen.fatal)
        self.assertFalse(dragen.passed)

    def test_orphan_hashes_fail(self):
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378,
                                  {'orphan_associations': 5, 'affected_samples': 2},
                                  _ONE_BLOB, expected_dragen_version='3.7.8')
        self.assertFalse(ok)
        self.assertFalse(self._names(checks)['referential_integrity'].passed)

    def test_shared_blob_is_informational(self):
        # multiple shared blobs must NOT fail the run (informational only)
        two_blobs = [{'blob_hash': 'a', 'n_samples': 2}, {'blob_hash': 'b', 'n_samples': 1}]
        ok, checks = cvh.evaluate(_summary(), _DRAGEN_378, _NO_ORPHANS, two_blobs,
                                  expected_dragen_version='3.7.8')
        self.assertTrue(ok)
        self.assertFalse(self._names(checks)['shared_blob_distribution'].fatal)

    def test_report_renders(self):
        ok, checks = cvh.evaluate(_summary(samples_not_reblocked=1, example_not_reblocked=['S2']),
                                  _DRAGEN_378, _NO_ORPHANS, _ONE_BLOB, expected_dragen_version='3.7.8')
        report = cvh.compose_report(PROJECT, DATASET, '3.7.8', ok, checks)
        self.assertIn("OVERALL: FAIL", report)
        self.assertIn("[FAIL] reblocking", report)
        self.assertIn("[INFO] shared_blob_distribution", report)


@unittest.skipUnless(
    os.environ.get("GVS_HEADER_IT_PROJECT") and os.environ.get("GVS_HEADER_IT_DATASET"),
    "Set GVS_HEADER_IT_PROJECT and GVS_HEADER_IT_DATASET (a throwaway dataset) to run the BigQuery "
    "integration test.",
)
class TestValidateIntegration(unittest.TestCase):
    """Real-BigQuery (or emulator) exercise of the full query stack + evaluation.

    Builds sample_info + the two header tables in a throwaway dataset, loads a small cohort, and runs
    run_checks end to end. If the backend (e.g. a limited BigQuery emulator) does not support the SQL
    features these queries rely on -- REGEXP_EXTRACT, CONTAINS_SUBSTR, ARRAY_AGG(... LIMIT), COUNTIF --
    the whole class is skipped in setUpClass so CI can never be broken by an emulator gap; a real
    BigQuery backend runs it fully.
    """

    @classmethod
    def setUpClass(cls):
        from google.cloud import bigquery
        cls.bigquery = bigquery
        cls.project = os.environ["GVS_HEADER_IT_PROJECT"]
        cls.dataset = os.environ["GVS_HEADER_IT_DATASET"]
        endpoint = os.environ.get("GVS_HEADER_IT_ENDPOINT")
        # utils.execute_with_retry reads client._default_query_job_config.labels, so the client must
        # be built with a default query job config (production's _make_client does the same).
        default_config = bigquery.QueryJobConfig(labels={}, use_legacy_sql=False)
        if endpoint:
            from google.api_core.client_options import ClientOptions
            from google.auth.credentials import AnonymousCredentials
            cls.client = bigquery.Client(
                project=cls.project,
                credentials=AnonymousCredentials(),
                client_options=ClientOptions(api_endpoint=endpoint),
                default_query_job_config=default_config,
            )
        else:
            cls.client = bigquery.Client(project=cls.project, default_query_job_config=default_config)

        # Probe: skip the class if the backend lacks any SQL feature the checks require.
        probe = ("SELECT REGEXP_EXTRACT('SW: 1.2.3', r'SW: [0-9.]+') AS v, "
                 "CONTAINS_SUBSTR('abc', 'b') AS c, "
                 "(SELECT ARRAY_AGG(x IGNORE NULLS LIMIT 2) FROM UNNEST([1, 2, 3]) AS x) AS a, "
                 "(SELECT COUNTIF(x > 1) FROM UNNEST([1, 2, 3]) AS x) AS n")
        try:
            cls.client.query(probe).result()
        except Exception as err:  # noqa: BLE001 -- any backend limitation should skip, not fail
            raise unittest.SkipTest(f"BigQuery backend lacks required SQL features; skipping: {err}")

    def _fq(self, table):
        return f"{self.project}.{self.dataset}.{table}"

    def _run(self, sql):
        self.client.query(sql).result()

    @staticmethod
    def _md5(text):
        return hashlib.md5(text.encode("utf-8")).hexdigest()

    def _create_tables(self):
        self._run(f"CREATE OR REPLACE TABLE `{self._fq('sample_info')}` "
                  "(sample_name STRING, sample_id INT64, is_loaded BOOL, is_control BOOL, withdrawn TIMESTAMP)")
        self._run(f"CREATE OR REPLACE TABLE `{self._fq('vcf_header_lines')}` "
                  "(vcf_header_lines_hash STRING, vcf_header_lines STRING, is_expected_unique BOOL)")
        self._run(f"CREATE OR REPLACE TABLE `{self._fq('sample_vcf_header')}` "
                  "(sample_id INT64, vcf_header_lines_hash STRING)")

    def _load_cohort(self, samples):
        """samples: dict sample_id -> list of (text, is_expected_unique_bool).

        Populates sample_info (non-control, non-withdrawn), the deduplicated vcf_header_lines, and the
        sample->hash map, exactly as ingest would.
        """
        si_values = ", ".join(f"('sample_{sid}', {sid}, true, false, NULL)" for sid in samples)
        self._run(f"INSERT INTO `{self._fq('sample_info')}` "
                  f"(sample_name, sample_id, is_loaded, is_control, withdrawn) VALUES {si_values}")

        header_lines = {}  # hash -> (text, is_unique)
        assocs = []        # (sample_id, hash)
        for sid, chunks in samples.items():
            for text, is_unique in chunks:
                h = self._md5(text)
                header_lines[h] = (text, is_unique)
                assocs.append((sid, h))

        hl_values = ", ".join(
            f"('{h}', '{text}', {str(is_unique).lower()})" for h, (text, is_unique) in header_lines.items())
        self._run(f"INSERT INTO `{self._fq('vcf_header_lines')}` "
                  f"(vcf_header_lines_hash, vcf_header_lines, is_expected_unique) VALUES {hl_values}")

        svh_values = ", ".join(f"({sid}, '{h}')" for sid, h in assocs)
        self._run(f"INSERT INTO `{self._fq('sample_vcf_header')}` "
                  f"(sample_id, vcf_header_lines_hash) VALUES {svh_values}")

    BLOB = '##fileformat=VCFv4.2 ##contig=<ID=chr1>'
    DRAGEN_378 = '##DRAGENCommandLine=<ID=dragen,Version="SW: 05.021.604.3.7.8">'
    DRAGEN_378_ALT = '##DRAGENCommandLine=<ID=dragen,Version="SW: 07.021.604.3.7.8">'
    DRAGEN_3412 = '##DRAGENCommandLine=<ID=dragen,Version="SW: 05.021.604.3.4.12">'
    REBLOCK = '##GATKCommandLine=<ID=ReblockGVCF,CommandLine="ReblockGVCF --output x.g.vcf">'

    def test_passes_on_good_cohort(self):
        self._create_tables()
        self._load_cohort({
            1: [(self.BLOB, False), (self.DRAGEN_378, True), (self.REBLOCK, True)],
            2: [(self.BLOB, False), (self.DRAGEN_378_ALT, True), (self.REBLOCK, True)],
        })
        ok, checks = cvh.run_checks(self.project, self.dataset,
                                    expected_dragen_version='3.7.8', client=self.client)
        self.assertTrue(ok, cvh.compose_report(self.project, self.dataset, '3.7.8', ok, checks))

    def test_fails_on_bad_cohort(self):
        self._create_tables()
        # sample 2: not reblocked AND a different DRAGEN triplet.
        self._load_cohort({
            1: [(self.BLOB, False), (self.DRAGEN_378, True), (self.REBLOCK, True)],
            2: [(self.BLOB, False), (self.DRAGEN_3412, True)],
        })
        ok, checks = cvh.run_checks(self.project, self.dataset,
                                    expected_dragen_version='3.7.8', client=self.client)
        self.assertFalse(ok)
        by_name = {c.name: c for c in checks}
        self.assertFalse(by_name['reblocking'].passed)
        self.assertFalse(by_name['dragen_version'].passed)

    @classmethod
    def tearDownClass(cls):
        for t in ("sample_info", "vcf_header_lines", "sample_vcf_header"):
            cls.client.query(f"DROP TABLE IF EXISTS `{cls.project}.{cls.dataset}.{t}`").result()


if __name__ == "__main__":
    unittest.main()
