"""
Unit tests for verify_structural_checks.py -- the independent, row-count-based checks added for
VS-1989, and the footer row count comparison added for VS-2032.

The BigQuery-reading functions are tested by patching verify_structural_checks.bigquery and asserting
on the generated SQL (mirroring the style in test_parquet_loading.py). The assessment functions are
pure and are tested directly. run_structural_checks is tested with the three query helpers patched so
no BigQuery is touched.
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
    assess_cross_family_consistency,
    assess_cardinality,
    assess_truncation_screen,
    assess_footer_row_counts,
    footer_source_uri,
    get_parquet_footer_row_counts,
    run_structural_checks,
    DEFAULT_VET_TRUNCATION_THRESHOLD,
    TRUNCATION_SCREEN_DISABLED,
)


def _row(**kwargs):
    """A stand-in BigQuery Row: attribute access returns the given values."""
    row = MagicMock()
    for k, val in kwargs.items():
        setattr(row, k, val)
    return row


OUTPUT_GCS_DIR = "gs://bucket/out"


def _run(expected_by_family, gcs_pairs_to_files=None, **kwargs):
    """run_structural_checks with the footer comparison's inputs defaulted. An empty GCS listing gives
    the footer comparison nothing to check, so tests not about it are unaffected by it."""
    return run_structural_checks("proj", "ds", expected_by_family, gcs_pairs_to_files or {},
                                 OUTPUT_GCS_DIR, **kwargs)


def _uri(table, sample_id):
    """A GCS URI shaped as GenerateParquetFiles uploads it."""
    family = table.rsplit("_", 1)[0] if table[-1].isdigit() else table
    subdir = f"{table.rsplit('_', 1)[1]}/" if table[-1].isdigit() else ""
    return f"{OUTPUT_GCS_DIR}/{family}/{subdir}{table}_{sample_id}_input.parquet"


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


class TestFooterSourceUri(unittest.TestCase):
    def test_scoped_to_the_family_directory_and_parquet_suffix(self):
        # The family directory keeps quarantine/ and sibling callsets out; the suffix keeps
        # .parquet.quarantined objects and the diagnostics JSON out.
        self.assertEqual(footer_source_uri("gs://b/out", "vet"), "gs://b/out/vet/*.parquet")

    def test_trailing_slash_is_not_doubled(self):
        self.assertEqual(footer_source_uri("gs://b/out/", "ref_ranges"), "gs://b/out/ref_ranges/*.parquet")

    def test_non_gcs_dir_rejected(self):
        with self.assertRaises(ValueError):
            footer_source_uri("/local/out", "vet")


class TestGetParquetFooterRowCounts(unittest.TestCase):
    def _client(self, mock_bq, rows_by_call):
        mock_client = MagicMock()
        mock_bq.Client.return_value = mock_client
        jobs = []
        for rows in rows_by_call:
            job = MagicMock()
            job.result.return_value = rows
            job.total_bytes_billed = 0
            job.slot_millis = 1
            jobs.append(job)
        mock_client.query.side_effect = jobs
        return mock_client

    @patch("verify_structural_checks.bigquery")
    def test_one_external_table_query_per_family(self, mock_bq):
        mock_client = self._client(mock_bq, [[], []])
        externals = []
        mock_bq.ExternalConfig.side_effect = lambda fmt: externals.append(MagicMock(format=fmt)) or externals[-1]

        get_parquet_footer_row_counts("proj", "gs://b/out", ["vet", "ref_ranges"])

        self.assertEqual(mock_client.query.call_count, 2)
        self.assertEqual([e.format for e in externals], ["PARQUET", "PARQUET"])
        self.assertEqual([e.source_uris for e in externals],
                         [["gs://b/out/vet/*.parquet"], ["gs://b/out/ref_ranges/*.parquet"]])
        # The external table is bound to the query job by name, and the query counts per file.
        table_definitions = mock_bq.QueryJobConfig.call_args_list[0].kwargs["table_definitions"]
        self.assertIs(table_definitions[vsc.FOOTER_EXTERNAL_TABLE], externals[0])
        query = mock_client.query.call_args_list[0][0][0]
        self.assertIn("_FILE_NAME", query)
        self.assertIn("COUNT(*)", query)
        self.assertIn(f"FROM {vsc.FOOTER_EXTERNAL_TABLE}", query)

    @patch("verify_structural_checks.bigquery")
    def test_returns_rows_per_file_across_families(self, mock_bq):
        self._client(mock_bq, [
            [_row(file_name="gs://b/out/vet/001/vet_001_1_x.parquet", num_rows=100)],
            [_row(file_name="gs://b/out/ref_ranges/001/ref_ranges_001_1_x.parquet", num_rows=50)],
        ])
        result = get_parquet_footer_row_counts("proj", "gs://b/out", ["vet", "ref_ranges"])
        self.assertEqual(result, {"gs://b/out/vet/001/vet_001_1_x.parquet": 100,
                                  "gs://b/out/ref_ranges/001/ref_ranges_001_1_x.parquet": 50})

    @patch("verify_structural_checks.bigquery")
    def test_no_families_skips_query(self, mock_bq):
        self.assertEqual(get_parquet_footer_row_counts("proj", "gs://b/out", []), {})
        mock_bq.Client.assert_not_called()


class TestAssessFooterRowCounts(unittest.TestCase):
    SUPER = ["vet", "ref_ranges"]
    REGULAR = ["sample_chromosome_ploidy", "vcf_header_lines_scratch"]
    FAMILIES = ["vet", "ref_ranges", "sample_chromosome_ploidy"]

    def _assess(self, file_rows, listing, partition_rows=(), regular_counts=None):
        return assess_footer_row_counts(file_rows, listing, list(partition_rows), regular_counts or {},
                                        self.FAMILIES, self.SUPER, self.REGULAR)

    def test_agreement_passes(self):
        listing = {("vet_001", 1): [_uri("vet_001", 1)],
                   ("ref_ranges_001", 1): [_uri("ref_ranges_001", 1)],
                   ("sample_chromosome_ploidy", 1): [_uri("sample_chromosome_ploidy", 1)]}
        file_rows = {_uri("vet_001", 1): 100, _uri("ref_ranges_001", 1): 50,
                     _uri("sample_chromosome_ploidy", 1): 24}
        r = self._assess(file_rows, listing, [("vet_001", 1, 100), ("ref_ranges_001", 1, 50)],
                         {"sample_chromosome_ploidy": {1: (24, 24)}})
        self.assertTrue(r["ok"])
        self.assertEqual(r["flagged_samples"], [])
        self.assertEqual({f: v["pairs_checked"] for f, v in r["per_family"].items()},
                         {"vet": 1, "ref_ranges": 1, "sample_chromosome_ploidy": 1})

    def test_more_rows_in_bigquery_is_flagged(self):
        listing = {("ref_ranges_001", 1): [_uri("ref_ranges_001", 1)]}
        r = self._assess({_uri("ref_ranges_001", 1): 50}, listing, [("ref_ranges_001", 1, 100)])
        self.assertFalse(r["ok"])
        self.assertEqual(r["flagged_samples"], [1])
        self.assertEqual(r["per_family"]["ref_ranges"]["mismatched"],
                         [{"table": "ref_ranges_001", "sample_id": 1, "footer_rows": 50, "bigquery_rows": 100}])

    def test_fewer_rows_in_bigquery_is_flagged(self):
        listing = {("vet_001", 1): [_uri("vet_001", 1)]}
        r = self._assess({_uri("vet_001", 1): 100}, listing, [("vet_001", 1, 99)])
        self.assertEqual(r["flagged_samples"], [1])

    def test_ploidy_compared_against_its_row_count(self):
        # The ploidy count comes from the regular-table read, which may carry a distinct-chromosome
        # count alongside; only the row count is compared.
        listing = {("sample_chromosome_ploidy", 1): [_uri("sample_chromosome_ploidy", 1)]}
        r = self._assess({_uri("sample_chromosome_ploidy", 1): 24}, listing,
                         regular_counts={"sample_chromosome_ploidy": {1: (48, 24)}})
        self.assertEqual(r["per_family"]["sample_chromosome_ploidy"]["mismatched"][0]["bigquery_rows"], 48)

    def test_two_files_for_one_pair_flag_even_when_the_sums_agree(self):
        # Two files each loaded once read 2N against 2N; only the file count exposes the duplicate.
        files = [_uri("vet_001", 1), _uri("vet_001", 1).replace("input", "input2")]
        listing = {("vet_001", 1): files, ("vet_001", 2): [_uri("vet_001", 2)]}
        r = self._assess({files[0]: 100, files[1]: 100, _uri("vet_001", 2): 100}, listing,
                         [("vet_001", 1, 200), ("vet_001", 2, 100)])
        self.assertFalse(r["ok"])
        self.assertEqual(r["flagged_samples"], [1])
        self.assertEqual(r["per_family"]["vet"]["multiple_files"],
                         [{"table": "vet_001", "sample_id": 1, "files": sorted(files)}])
        self.assertEqual(r["per_family"]["vet"]["mismatched"], [])

    def test_one_file_per_pair_reports_no_multiple_files(self):
        listing = {("vet_001", 1): [_uri("vet_001", 1)]}
        r = self._assess({_uri("vet_001", 1): 100}, listing, [("vet_001", 1, 100)])
        self.assertTrue(r["ok"])
        self.assertEqual(r["per_family"]["vet"]["multiple_files"], [])

    def test_sample_absent_from_bigquery_compares_against_zero(self):
        listing = {("vet_001", 1): [_uri("vet_001", 1)]}
        r = self._assess({_uri("vet_001", 1): 100}, listing)
        self.assertEqual(r["per_family"]["vet"]["mismatched"][0]["bigquery_rows"], 0)

    def test_listed_file_the_query_did_not_read_flags_its_sample(self):
        # A sample whose footer could not be counted cannot be cleared.
        listing = {("vet_001", 1): [_uri("vet_001", 1)], ("vet_001", 2): [_uri("vet_001", 2)]}
        r = self._assess({_uri("vet_001", 1): 100}, listing, [("vet_001", 1, 100), ("vet_001", 2, 100)])
        self.assertFalse(r["ok"])
        self.assertEqual(r["flagged_samples"], [2])
        self.assertEqual(r["per_family"]["vet"]["unread"],
                         [{"table": "vet_001", "sample_id": 2, "files": [_uri("vet_001", 2)]}])
        self.assertEqual(r["per_family"]["vet"]["mismatched"], [])

    def test_unlisted_file_is_reported_but_not_counted(self):
        # A file that appeared after discovery must not inflate the sample's sum.
        listing = {("vet_001", 1): [_uri("vet_001", 1)]}
        stray = _uri("vet_001", 1).replace("input", "late")
        r = self._assess({_uri("vet_001", 1): 100, stray: 100}, listing, [("vet_001", 1, 100)])
        self.assertTrue(r["ok"])
        self.assertEqual(r["unlisted_files"], [stray])

    def test_unchecked_family_is_skipped(self):
        # vcf_header_lines_scratch is not a footer-checked family; its files are not expected to have
        # been read, so they must not read as unread.
        listing = {("vcf_header_lines_scratch", 1): [f"{OUTPUT_GCS_DIR}/vcf_header_lines_scratch/x_1.parquet"]}
        r = self._assess({}, listing, regular_counts={"vcf_header_lines_scratch": {1: 10}})
        self.assertTrue(r["ok"])
        self.assertNotIn("vcf_header_lines_scratch", r["per_family"])


class TestGetPloidyRowCounts(unittest.TestCase):
    @patch("verify_structural_checks.bigquery")
    def test_query_groups_by_sample(self, mock_bq):
        mock_client = MagicMock()
        mock_bq.Client.return_value = mock_client
        mock_client.query.return_value = [
            _row(sample_id=1, n=24, n_distinct=24),
            _row(sample_id=2, n=24, n_distinct=24),
        ]

        result = get_ploidy_row_counts("proj", "ds", "sample_chromosome_ploidy")

        query = mock_client.query.call_args[0][0]
        self.assertIn("sample_chromosome_ploidy", query)
        self.assertIn("GROUP BY sample_id", query)
        self.assertIn("COUNT(*)", query)
        self.assertIn("COUNT(DISTINCT chromosome)", query)
        self.assertEqual(result, {1: (24, 24), 2: (24, 24)})

    def test_invalid_table_raises(self):
        with self.assertRaises(ValueError):
            get_ploidy_row_counts("proj", "ds", "bad table")


class TestAssessCardinality(unittest.TestCase):
    def test_uniform_counts_pass(self):
        r = assess_cardinality({i: (24, 24) for i in range(1, 6)}, set(range(1, 6)))
        self.assertTrue(r["ok"])
        self.assertEqual(r["distinct_samples"], 5)
        self.assertEqual(r["deviating_samples"], [])
        self.assertEqual(r["total_rows"], 120)

    def test_does_not_assume_24(self):
        # A callset whose per-sample count is 25 (e.g. chrM ingested) passes -- no count is privileged.
        r = assess_cardinality({i: (25, 25) for i in range(1, 4)}, set(range(1, 4)))
        self.assertTrue(r["ok"])

    def test_missing_sample_fails(self):
        r = assess_cardinality({1: (24, 24), 2: (24, 24)}, {1, 2, 3})
        self.assertFalse(r["ok"])
        self.assertEqual(r["missing_samples"], [3])

    def test_zero_count_treated_as_missing(self):
        r = assess_cardinality({1: (24, 24), 2: (0, 0)}, {1, 2})
        self.assertFalse(r["ok"])
        self.assertEqual(r["missing_samples"], [2])

    def test_extra_samples_in_table_ignored(self):
        # A sample already in the table from a prior load (99) is not expected this run: ignored.
        r = assess_cardinality({1: (24, 24), 2: (24, 24), 99: (48, 24)}, {1, 2})
        self.assertTrue(r["ok"])
        self.assertEqual(r["distinct_samples"], 2)

    def test_no_expected_is_ok(self):
        r = assess_cardinality({}, set())
        self.assertTrue(r["ok"])
        self.assertIsNone(r["min"])
        self.assertIsNone(r["max"])

    def test_duplication_flagged(self):
        # Sample 2 has 48 rows across only 24 distinct chromosomes (duplicated load).
        r = assess_cardinality({1: (24, 24), 2: (48, 24), 3: (24, 24)}, {1, 2, 3})
        self.assertFalse(r["ok"])
        self.assertEqual(r["deviating_samples"],
                         [{"sample_id": 2, "count": 48, "distinct_chromosomes": 24}])

    def test_singleton_duplication_flagged(self):
        # A single-sample cohort has no cohort consensus to appeal to, and needs none: the exact
        # comparison is within the sample.
        r = assess_cardinality({1: (48, 24)}, {1})
        self.assertFalse(r["ok"])
        self.assertEqual([d["sample_id"] for d in r["deviating_samples"]], [1])

    def test_whole_callset_duplication_flagged(self):
        # Every sample doubled -- the case a modal reference could never see, since the mode doubles too.
        r = assess_cardinality({1: (48, 24), 2: (48, 24)}, {1, 2})
        self.assertFalse(r["ok"])
        self.assertEqual([d["sample_id"] for d in r["deviating_samples"]], [1, 2])

    def test_heterogeneous_contig_counts_pass(self):
        # Legitimate contig coverage variations (23 where the gVCF omits a contig, 25 with chrM) pass.
        r = assess_cardinality({1: (23, 23), 2: (24, 24), 3: (25, 25)}, {1, 2, 3})
        self.assertTrue(r["ok"])
        self.assertEqual(r["min"], 23)
        self.assertEqual(r["max"], 25)
        self.assertEqual(r["deviating_samples"], [])

    def test_low_outlier_is_not_flagged(self):
        # The accepted gap (see assess_cardinality): a sample with far fewer rows than its peers is NOT
        # flagged, because a smaller contig set is also what a legitimate exome or karyotype looks like.
        # Whole-sample loss is caught by completeness instead.
        r = assess_cardinality({1: (24, 24), 2: (20, 20), 3: (24, 24)}, {1, 2, 3})
        self.assertTrue(r["ok"])
        self.assertEqual(r["deviating_samples"], [])

    def test_counts_without_distinct_degenerate_to_presence(self):
        # A regular table with no per-sample-unique column yields raw counts; nothing exact can be said
        # about them, so only presence is assessed.
        r = assess_cardinality({1: 24, 2: 48}, {1, 2})
        self.assertTrue(r["ok"])
        self.assertEqual(r["deviating_samples"], [])

    def test_dict_shaped_counts_accepted(self):
        r = assess_cardinality({1: {"count": 48, "distinct": 24}}, {1})
        self.assertFalse(r["ok"])
        self.assertEqual([d["sample_id"] for d in r["deviating_samples"]], [1])


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
        # A present-but-empty partition is a partial load, not an absence: it must not also be
        # double-counted as missing (regression guard for the double-count fix).
        self.assertEqual(r["per_family"]["vet"]["missing_samples"], [])
        self.assertEqual(r["per_family"]["vet"]["present"], 1)

    def test_missing_and_empty_are_disjoint(self):
        # sample 2 present (100 rows), sample 3 empty (0 rows), sample 4 absent entirely. Each must
        # land in exactly one bucket.
        part = [("vet_001", 2, 100), ("vet_001", 3, 0)]
        reg = {"sample_chromosome_ploidy": {2: 24, 3: 24, 4: 24}}
        exp = {"vet": {2, 3, 4}, "ref_ranges": set(), "sample_chromosome_ploidy": {2, 3, 4}}
        r = self._run(part, reg, exp)
        vet = r["per_family"]["vet"]
        self.assertFalse(vet["ok"])
        self.assertEqual(vet["missing_samples"], [4])
        self.assertEqual(vet["empty_partition_samples"], [3])
        self.assertEqual(vet["present"], 1)

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


class TestAssessCrossFamilyConsistency(unittest.TestCase):
    """The cross-family gaps assess_family_completeness judges vacuously: within the co-produced data
    group (vet / ref_ranges / ploidy, emitted together per sample), a sample present in some members but
    absent from another, and an entirely absent member. The check is presence-activated -- dormant unless
    a group member has files this run -- so a headers-only ingest, which produces none of these, passes."""

    CO_PRODUCED = ["vet", "ref_ranges", "sample_chromosome_ploidy"]

    def test_identical_family_sets_pass(self):
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = assess_cross_family_consistency(exp, self.CO_PRODUCED)
        self.assertTrue(r["ok"])
        self.assertEqual(r["union_size"], 2)

    def test_sample_missing_from_one_family_flagged(self):
        # sample 2 was produced for vet and ploidy but its ref_ranges file never landed -- the per-sample
        # cross-family gap. completeness would pass ref_ranges vacuously (it only expects {1}); this
        # catches it.
        exp = {"vet": {1, 2}, "ref_ranges": {1}, "sample_chromosome_ploidy": {1, 2}}
        r = assess_cross_family_consistency(exp, self.CO_PRODUCED)
        self.assertFalse(r["ok"])
        self.assertEqual(r["union_size"], 2)
        self.assertEqual(r["per_family"]["ref_ranges"]["missing_samples"], [2])
        self.assertTrue(r["per_family"]["vet"]["ok"])
        self.assertTrue(r["per_family"]["sample_chromosome_ploidy"]["ok"])

    def test_entirely_absent_member_flagged(self):
        # No vet file was produced for any sample, so vet has no key in expected_by_family -- but a sibling
        # member (ref_ranges / ploidy) is present, so the check activates and vet is compared as the empty
        # set and flagged missing every sample the other members carry. (Regression for the whole-family
        # hole: a union over only present families would drop vet and pass, authorizing deletion of Parquet
        # that lacks all variant data.)
        exp = {"ref_ranges": {1, 2, 3}, "sample_chromosome_ploidy": {1, 2, 3}}
        r = assess_cross_family_consistency(exp, self.CO_PRODUCED)
        self.assertFalse(r["ok"])
        self.assertEqual(r["union_size"], 3)
        self.assertEqual(r["per_family"]["vet"]["missing_samples"], [1, 2, 3])
        self.assertTrue(r["per_family"]["ref_ranges"]["ok"])
        self.assertTrue(r["per_family"]["sample_chromosome_ploidy"]["ok"])

    def test_headers_only_run_with_scratch_samples_passes(self):
        # A supported headers-only ingest lists only vcf_header_lines_scratch in GCS -- with real samples --
        # while the data prefixes are configured but produce nothing. vcf_header_lines_scratch is not a
        # co-produced member, so no member is present, the check is dormant, and the run passes. (Regression
        # for the R10 whole-family fix rejecting the headers-only path: taking the union over every
        # configured prefix reported every header sample missing from all three data families.)
        exp = {"vcf_header_lines_scratch": {1, 2, 3}}
        r = assess_cross_family_consistency(exp, self.CO_PRODUCED)
        self.assertTrue(r["ok"])
        self.assertEqual(r["union_size"], 0)
        self.assertEqual(r["per_family"], {})

    def test_empty_run_passes(self):
        # No family produced anything at all: no member present, dormant, passes.
        r = assess_cross_family_consistency({}, self.CO_PRODUCED)
        self.assertTrue(r["ok"])
        self.assertEqual(r["union_size"], 0)


class TestScreensNameFlaggedSamples(unittest.TestCase):
    """The truncation screen reports the samples it objects to in flagged_samples, and
    run_structural_checks unions them with the footer comparison's. That list is what the workflow quarantines, so it must be exactly the flagged samples
    -- too few and flagged Parquet is deleted, too many and clean Parquet is needlessly retained."""

    PREFIXES = ["vet", "ref_ranges"]

    def test_truncation_outliers_are_named(self):
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 1)]
        r = assess_truncation_screen(part, "vet", set(range(1, 10)), self.PREFIXES, 1.6)
        self.assertEqual(r["flagged_samples"], [9])

    def test_clean_screen_names_nothing(self):
        part = [("vet_001", i, 100) for i in range(1, 6)]
        self.assertEqual(
            assess_truncation_screen(part, "vet", set(range(1, 6)), self.PREFIXES, 1.6)["flagged_samples"], [])

    def test_singleton_names_the_one_sample_it_could_not_clear(self):
        # A singleton load has no peers, so the screen cannot clear its one sample and names it -- there
        # are no outliers to report, but the sample is still what gets quarantined.
        part = [("vet_001", 7, 100)]
        r = assess_truncation_screen(part, "vet", {7}, self.PREFIXES, 1.6)
        self.assertTrue(r["singleton_flagged"])
        self.assertEqual(r["outliers"], [])
        self.assertEqual(r["flagged_samples"], [7])

    def test_unscreened_family_names_nothing(self):
        part = [("vet_001", i, 100) for i in range(1, 6)]
        r = assess_truncation_screen(part, "ref_ranges", set(range(1, 6)), self.PREFIXES, 1.6)
        self.assertEqual(r["flagged_samples"], [])

    def test_multiple_outliers_are_sorted_by_sample_id(self):
        # outliers themselves are ordered by ratio; flagged_samples is a set to quarantine, so it is
        # ordered by sample_id for stable diffing of the work list.
        part = [("vet_001", i, 100) for i in range(1, 9) if i != 3]
        part += [("vet_001", 9, 10), ("vet_001", 3, 20)]
        r = assess_truncation_screen(part, "vet", set(range(1, 10)), self.PREFIXES, 1.6)
        self.assertEqual([o["sample_id"] for o in r["outliers"]], [9, 3])
        self.assertEqual(r["flagged_samples"], [3, 9])


class TestAssessTruncationScreen(unittest.TestCase):
    """Flags a grossly under-rowed vet partition."""

    def test_flags_low_outlier(self):
        # sample 9 has 1 row where its peers have 100 -- a truncated partition. floor = 100/1.6 = 62.5.
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 1)]
        r = assess_truncation_screen(part, "vet", set(range(1, 10)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["median"], 100)
        self.assertEqual(r["baseline"], 100)
        self.assertEqual([o["sample_id"] for o in r["outliers"]], [9])
        self.assertAlmostEqual(r["outliers"][0]["ratio"], 0.01)

    def test_uniform_no_outliers(self):
        part = [("vet_001", i, 100) for i in range(1, 6)]
        r = assess_truncation_screen(part, "vet", set(range(1, 6)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["outliers"], [])

    def test_mild_dip_within_threshold_not_flagged(self):
        # A sample at 70% of the median is above the floor (median/1.6 = 62.5%) and must not flag.
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 70)]
        r = assess_truncation_screen(part, "vet", set(range(1, 10)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["outliers"], [])

    def test_zero_row_partition_not_screened(self):
        # A 0-row partition is the completeness check's empty case, not truncation's: it is excluded
        # here (only present, non-empty partitions are judged) so it is not double-reported.
        part = [("vet_001", i, 100) for i in range(1, 6)] + [("vet_001", 6, 0)]
        r = assess_truncation_screen(part, "vet", set(range(1, 7)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["samples_screened"], 5)
        self.assertEqual(r["outliers"], [])

    def test_other_family_not_screened(self):
        part = [("vet_001", i, 100) for i in range(1, 6)]
        r = assess_truncation_screen(part, "ref_ranges", set(range(1, 6)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["samples_screened"], 0)
        self.assertEqual(r["outliers"], [])

    def test_only_expected_samples_screened(self):
        part = [("vet_001", i, 100) for i in range(1, 6)] + [("vet_001", 99, 1)]
        r = assess_truncation_screen(part, "vet", set(range(1, 6)), ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["samples_screened"], 5)
        self.assertEqual(r["outliers"], [])

    def test_singleton_is_flagged(self):
        part = [("vet_001", 1, 100)]
        r = assess_truncation_screen(part, "vet", {1}, ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["samples_screened"], 1)
        self.assertEqual(r["outliers"], [])
        self.assertTrue(r["singleton_flagged"])

    def test_two_sample_load_with_truncation_flagged(self):
        # Regression for two-sample median averaging: with row counts 100 and 200, median is 150
        # which would give ratio 0.667 (> 0.625). Using the upper baseline catches sample 1 at 0.5x.
        part = [("vet_001", 1, 100), ("vet_001", 2, 200)]
        r = assess_truncation_screen(part, "vet", {1, 2}, ["vet", "ref_ranges"], 1.6)
        self.assertEqual(r["samples_screened"], 2)
        self.assertEqual(r["median"], 150)
        self.assertEqual(r["baseline"], 200)
        self.assertEqual([o["sample_id"] for o in r["outliers"]], [1])
        self.assertAlmostEqual(r["outliers"][0]["ratio"], 0.5)


class TestRunStructuralChecks(unittest.TestCase):
    """run_structural_checks with the three BigQuery helpers patched (no BQ touched)."""

    def _patch(self, partition_rows, ploidy_counts, file_rows=None):
        self._patch_regular(partition_rows, {"sample_chromosome_ploidy": ploidy_counts}, file_rows)

    def _patch_regular(self, partition_rows, counts_by_table, file_rows=None):
        """Like _patch but returns per-table regular counts, keyed on the table-name argument."""
        p1 = patch("verify_structural_checks.get_partition_row_counts", return_value=partition_rows)
        p2 = patch(
            "verify_structural_checks.get_ploidy_row_counts",
            side_effect=lambda project_id, dataset_name, table: counts_by_table.get(table, {}),
        )
        self.footer_query = MagicMock(return_value=file_rows or {})
        p3 = patch("verify_structural_checks.get_parquet_footer_row_counts", self.footer_query)
        for p in (p1, p2, p3):
            p.start()
            self.addCleanup(p.stop)

    def test_all_good(self):
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: 24, 2: 24})
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])
        self.assertTrue(r["cross_family_ok"])
        self.assertFalse(r["row_count_mismatch_flagged"])
        self.assertTrue(r["details"]["footer_row_counts"]["ok"])
        self.assertIn("backfill_caveat", r["details"]["cardinality"]["sample_chromosome_ploidy"])

    def _listing(self, samples, ploidy=True):
        """A GCS listing of one vet, one ref_ranges and (optionally) one ploidy file per sample."""
        tables = ["vet_001", "ref_ranges_001"] + (["sample_chromosome_ploidy"] if ploidy else [])
        return {(t, sid): [_uri(t, sid)] for t in tables for sid in samples}

    def test_flagged_samples_is_the_union_of_footer_and_truncation(self):
        # Sample 3's ref_ranges was loaded twice (footer 50, BigQuery 100); sample 4's vet partition is
        # truncated (1 against peers of 100, and its footer agrees). Each check sees one of them; the
        # quarantine list must be both.
        samples = range(1, 5)
        part = ([("vet_001", i, 100) for i in (1, 2, 3)] + [("vet_001", 4, 1)]
                + [("ref_ranges_001", i, 50) for i in (1, 2, 4)] + [("ref_ranges_001", 3, 100)])
        listing = self._listing(samples)
        file_rows = {_uri("vet_001", i): (1 if i == 4 else 100) for i in samples}
        file_rows.update({_uri("ref_ranges_001", i): 50 for i in samples})
        file_rows.update({_uri("sample_chromosome_ploidy", i): 24 for i in samples})
        self._patch(part, {i: (24, 24) for i in samples}, file_rows)
        exp = {"vet": set(samples), "ref_ranges": set(samples),
               "sample_chromosome_ploidy": set(samples)}
        r = _run(exp, listing)
        self.assertTrue(r["row_count_mismatch_flagged"])
        self.assertTrue(r["truncation_flagged"])
        self.assertEqual(r["flagged_samples"], [3, 4])
        # The mismatch subset is what GvsImportGenomes aborts on, so the truncated sample must not
        # appear in it -- the truncation threshold is uncalibrated and must not fail a completed ingest.
        self.assertEqual(r["row_count_mismatch_samples"], [3])
        self.assertEqual(r["truncation_flagged_samples"], [4])
        # Neither touches the exact checks -- and so all_loaded.
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])
        self.assertTrue(r["cross_family_ok"])

    def test_doubled_ref_ranges_is_flagged(self):
        # The gap VS-2032 exists to close: a ref_ranges partition loaded twice, with every other family
        # clean. Nothing but the footer comparison can see it.
        samples = range(1, 4)
        part = ([("vet_001", i, 100) for i in samples]
                + [("ref_ranges_001", 1, 50), ("ref_ranges_001", 2, 7000), ("ref_ranges_001", 3, 14000)])
        file_rows = {_uri("vet_001", i): 100 for i in samples}
        file_rows.update({_uri("ref_ranges_001", 1): 50, _uri("ref_ranges_001", 2): 7000,
                          _uri("ref_ranges_001", 3): 7000})
        file_rows.update({_uri("sample_chromosome_ploidy", i): 24 for i in samples})
        self._patch(part, {i: (24, 24) for i in samples}, file_rows)
        exp = {"vet": set(samples), "ref_ranges": set(samples),
               "sample_chromosome_ploidy": set(samples)}
        r = _run(exp, self._listing(samples))
        self.assertEqual(r["row_count_mismatch_samples"], [3])
        mismatched = r["details"]["footer_row_counts"]["per_family"]["ref_ranges"]["mismatched"]
        self.assertEqual(mismatched, [{"table": "ref_ranges_001", "sample_id": 3,
                                       "footer_rows": 7000, "bigquery_rows": 14000}])

    def test_footer_query_scoped_to_configured_families_with_files(self):
        # A wildcard URI matching nothing fails the BigQuery query, so a family with no files this run
        # must not be queried; nor may one this run does not load.
        part = [("vet_001", 1, 100)]
        self._patch(part, {})
        _run({"vet": {1}}, {("vet_001", 1): [_uri("vet_001", 1)]},
             superpartitioned_table_prefixes=["vet", "ref_ranges"])
        self.footer_query.assert_called_once_with("proj", OUTPUT_GCS_DIR, ["vet"])

    def test_headers_only_run_queries_no_footers(self):
        self._patch_regular([], {"vcf_header_lines_scratch": {1: 10}})
        _run({"vcf_header_lines_scratch": {1}},
             regular_table_prefixes=["sample_chromosome_ploidy", "vcf_header_lines_scratch"])
        self.footer_query.assert_called_once_with("proj", OUTPUT_GCS_DIR, [])

    def test_clean_run_flags_no_samples(self):
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: 24, 2: 24})
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertEqual(r["flagged_samples"], [])
        self.assertEqual(r["row_count_mismatch_samples"], [])

    def test_waived_screens_still_report_flagged_samples(self):
        # allow_flagged_vet_loads is a policy decision made downstream in compute_quarantine, so the
        # truncation screen keeps reporting what it saw and the log can say what is being waived.
        part = [("vet_001", i, 100) for i in range(1, 3)] + [("vet_001", 3, 1)]
        part += [("ref_ranges_001", i, 50) for i in range(1, 4)]
        self._patch(part, {i: 24 for i in range(1, 4)})
        exp = {"vet": {1, 2, 3}, "ref_ranges": {1, 2, 3}, "sample_chromosome_ploidy": {1, 2, 3}}
        r = _run(exp, allow_flagged_vet_loads=True)
        self.assertEqual(r["flagged_samples"], [3])
        self.assertEqual(r["truncation_flagged_samples"], [3])
        self.assertTrue(r["allow_flagged_vet_loads"])

    def test_cross_family_gap_fails(self):
        # sample 2's ref_ranges file was never produced, so GCS lists it only for vet/ploidy. Every
        # listed pair is present in BigQuery, so completeness passes (ref_ranges only expected {1});
        # the cross-family check catches the gap.
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", 1, 50)]
        self._patch(part, {1: 24, 2: 24})
        exp = {"vet": {1, 2}, "ref_ranges": {1}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertTrue(r["completeness_ok"])
        self.assertFalse(r["cross_family_ok"])
        self.assertEqual(
            r["details"]["cross_family_consistency"]["per_family"]["ref_ranges"]["missing_samples"], [2]
        )

    def test_entirely_absent_family_fails(self):
        # No vet files at all; ref_ranges + ploidy present for every sample. Completeness passes (nothing
        # is "expected" for vet, so it is never iterated), but its co-produced siblings are present, so the
        # cross-family check activates, compares vet as the empty set, and fails all_loaded -- blocking
        # deletion of Parquet that lacks all variant data.
        part = [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: 24, 2: 24})
        exp = {"ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertTrue(r["completeness_ok"])
        self.assertFalse(r["cross_family_ok"])
        self.assertEqual(
            r["details"]["cross_family_consistency"]["per_family"]["vet"]["missing_samples"], [1, 2]
        )

    def test_headers_only_run_passes(self):
        # A supported headers-only ingest: the WDL configures the data prefixes (vet / ref_ranges /
        # sample_chromosome_ploidy) alongside vcf_header_lines_scratch, but only header files land, so GCS
        # lists samples for vcf_header_lines_scratch alone. Completeness passes (only that family is
        # expected and every sample has header rows), cardinality skips it (variable per-sample counts),
        # and the cross-family check stays dormant because no co-produced member is present -- so the run
        # is not falsely reported incomplete. (Regression for the R10 whole-family fix, which took the
        # union over every configured prefix and reported every header sample missing from all three data
        # families.)
        self._patch_regular(
            [],
            {
                "vcf_header_lines_scratch": {1: 10, 2: 12, 3: 8},
                "sample_chromosome_ploidy": {},
            },
        )
        exp = {"vcf_header_lines_scratch": {1, 2, 3}}
        r = _run(
            exp,
            superpartitioned_table_prefixes=["vet", "ref_ranges"],
            regular_table_prefixes=["sample_chromosome_ploidy", "vcf_header_lines_scratch"],
        )
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])
        self.assertTrue(r["cross_family_ok"])
        self.assertEqual(r["details"]["cross_family_consistency"]["per_family"], {})

    def test_combined_header_and_data_missing_header_fails_cross_family(self):
        # When headers and data are loaded together (e.g. GvsQuickstartIntegration with load_vcf_headers=true),
        # vcf_header_lines_scratch is included in the cross-family check. If sample 3's header Parquet never
        # landed, cross_family_ok fails so safe_to_delete_parquet is blocked.
        part = ([("vet_001", i, 100) for i in (1, 2, 3)]
                + [("ref_ranges_001", i, 50) for i in (1, 2, 3)])
        counts_by_table = {
            "sample_chromosome_ploidy": {1: 24, 2: 24, 3: 24},
            "vcf_header_lines_scratch": {1: 30, 2: 45},
        }
        self._patch_regular(part, counts_by_table)
        exp = {"vet": {1, 2, 3}, "ref_ranges": {1, 2, 3},
               "sample_chromosome_ploidy": {1, 2, 3}, "vcf_header_lines_scratch": {1, 2}}

        r = _run(
            exp,
            superpartitioned_table_prefixes=["vet", "ref_ranges"],
            regular_table_prefixes=["sample_chromosome_ploidy", "vcf_header_lines_scratch"],
        )
        self.assertFalse(r["cross_family_ok"])
        self.assertEqual(
            r["details"]["cross_family_consistency"]["per_family"]["vcf_header_lines_scratch"]["missing_samples"],
            [3],
        )

    def test_heterogeneous_ploidy_passes(self):
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: (24, 24), 2: (23, 23)})  # sample 2's gVCF legitimately covers 23 contigs
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])

    def test_low_ploidy_outlier_passes_cardinality(self):
        # The accepted gap: a sample with fewer ploidy rows than its peers and no repeated chromosome
        # passes, because a smaller contig set is indistinguishable from a legitimate one by row count.
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: (24, 24), 2: (20, 20)})  # sample 2 has 20 contigs, none repeated
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])
        self.assertEqual(r["details"]["cardinality"]["sample_chromosome_ploidy"]["deviating_samples"], [])

    def test_duplicated_ploidy_fails_cardinality_by_chromosome_collision(self):
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: (24, 24), 2: (48, 24)})  # sample 2 duplicated ploidy load (48 rows, 24 chromosomes)
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertTrue(r["completeness_ok"])
        self.assertFalse(r["cardinality_ok"])
        self.assertEqual(
            r["details"]["cardinality"]["sample_chromosome_ploidy"]["deviating_samples"][0]["sample_id"],
            2,
        )

    def test_nonzero_truncated_vet_partition_flagged_not_gated_by_default(self):
        # Regression for the Copilot finding: a vet partition present with a *nonzero* but grossly
        # truncated row count (1 where peers have 100) passes completeness (it is not empty), yet the
        # truncation screen surfaces it. The exact checks stay green (so all_loaded stays true); the
        # flag only quarantines the sample, unless --allow-flagged-vet-loads is set.
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 1)]
        part += [("ref_ranges_001", i, 50) for i in range(1, 10)]
        self._patch(part, {i: 24 for i in range(1, 10)})
        exp = {"vet": set(range(1, 10)), "ref_ranges": set(range(1, 10)),
               "sample_chromosome_ploidy": set(range(1, 10))}

        r = _run(exp, allow_flagged_vet_loads=False)
        self.assertTrue(r["truncation_flagged"])
        # The hard gates do not see it -- present and non-empty, uniform ploidy.
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])
        truncation = r["details"]["truncation_screen"]["vet"]
        self.assertEqual([o["sample_id"] for o in truncation["outliers"]], [9])

    def test_singleton_is_flagged_by_run_structural_checks(self):
        # A single-sample load has no cohort consensus, so the truncation screen cannot clear it.
        part = [("vet_001", 1, 100), ("ref_ranges_001", 1, 50)]
        self._patch(part, {1: 24})
        exp = {"vet": {1}, "ref_ranges": {1}, "sample_chromosome_ploidy": {1}}

        r = _run(exp, allow_flagged_vet_loads=False)
        self.assertTrue(r["truncation_flagged"])
        self.assertEqual(r["flagged_samples"], [1])
        # Factual gates remain ok (as long as presence, family completeness, and ploidy cardinality are good)
        self.assertTrue(r["completeness_ok"])
        self.assertTrue(r["cardinality_ok"])

    def test_truncation_screen_covers_vet_only(self):
        # ref_ranges is deliberately unscreened (wide per-sample distribution), so its detail is a
        # zero-sample screen even when a ref_ranges partition is tiny.
        part = ([("vet_001", i, 100) for i in (1, 2)]
                + [("ref_ranges_001", 1, 50), ("ref_ranges_001", 2, 1)])
        self._patch(part, {1: 24, 2: 24})
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertFalse(r["truncation_flagged"])
        self.assertIn("vet", r["details"]["truncation_screen"])
        self.assertNotIn("ref_ranges", r["details"]["truncation_screen"])

    def test_uniformly_duplicated_ploidy_fails_cardinality(self):
        # Every sample doubled -- the case no modal reference could catch, since the mode doubles with
        # the callset. The within-sample comparison catches it and gates deletion.
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        self._patch(part, {1: (48, 24), 2: (48, 24)})
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2}, "sample_chromosome_ploidy": {1, 2}}
        r = _run(exp)
        self.assertFalse(r["cardinality_ok"])
        card = r["details"]["cardinality"]["sample_chromosome_ploidy"]
        self.assertEqual([d["sample_id"] for d in card["deviating_samples"]], [1, 2])

    def test_variable_cardinality_regular_table_not_gated(self):
        # With VCF-header loading, the WDL passes vcf_header_lines_scratch alongside ploidy as a
        # regular table. Its per-sample row count legitimately varies (one row per header line), so it
        # must be completeness-checked but NOT held to the uniform-cardinality mode -- otherwise a valid
        # headers ingest fails. Regression guard for the header-enabled prefix list.
        part = ([("vet_001", i, 100) for i in (1, 2, 3)]
                + [("ref_ranges_001", i, 50) for i in (1, 2, 3)])
        counts_by_table = {
            "sample_chromosome_ploidy": {1: 24, 2: 24, 3: 24},   # uniform
            "vcf_header_lines_scratch": {1: 30, 2: 45, 3: 12},   # legitimately non-uniform
        }
        self._patch_regular(part, counts_by_table)
        exp = {"vet": {1, 2, 3}, "ref_ranges": {1, 2, 3},
               "sample_chromosome_ploidy": {1, 2, 3}, "vcf_header_lines_scratch": {1, 2, 3}}

        r = _run(
            exp,
            regular_table_prefixes=["sample_chromosome_ploidy", "vcf_header_lines_scratch"],
        )

        # The varying header counts must not fail cardinality, because the header table is never
        # cardinality-checked -- only ploidy is.
        self.assertTrue(r["cardinality_ok"])
        self.assertIn("sample_chromosome_ploidy", r["details"]["cardinality"])
        self.assertNotIn("vcf_header_lines_scratch", r["details"]["cardinality"])
        # Completeness still covers the header table (all three samples present).
        self.assertTrue(r["completeness_ok"])
        self.assertIn("vcf_header_lines_scratch", r["details"]["family_completeness"]["per_family"])

    def test_missing_header_sample_still_fails_completeness(self):
        # Dropping the header table from the cardinality check must not blind us to a lost header load:
        # a sample present in ploidy but absent from headers still fails completeness.
        part = [("vet_001", i, 100) for i in (1, 2)] + [("ref_ranges_001", i, 50) for i in (1, 2)]
        counts_by_table = {
            "sample_chromosome_ploidy": {1: 24, 2: 24},
            "vcf_header_lines_scratch": {1: 30},   # sample 2 missing
        }
        self._patch_regular(part, counts_by_table)
        exp = {"vet": {1, 2}, "ref_ranges": {1, 2},
               "sample_chromosome_ploidy": {1, 2}, "vcf_header_lines_scratch": {1, 2}}

        r = _run(
            exp,
            regular_table_prefixes=["sample_chromosome_ploidy", "vcf_header_lines_scratch"],
        )

        self.assertFalse(r["completeness_ok"])
        header = r["details"]["family_completeness"]["per_family"]["vcf_header_lines_scratch"]
        self.assertEqual(header["missing_samples"], [2])
        self.assertTrue(r["cardinality_ok"])  # ploidy still uniform


class TestTruncationThreshold(unittest.TestCase):
    """VS-1989: the truncation floor is uncalibrated, so it has to be retunable and switchable off."""

    PREFIXES = ["vet", "ref_ranges"]

    def _patch(self, partition_rows, ploidy_counts):
        for p in (patch("verify_structural_checks.get_partition_row_counts", return_value=partition_rows),
                  patch("verify_structural_checks.get_ploidy_row_counts", return_value=ploidy_counts),
                  patch("verify_structural_checks.get_parquet_footer_row_counts", return_value={})):
            p.start()
            self.addCleanup(p.stop)

    def _truncated_callset(self):
        # Sample 9 at 1 row against peers of 100 is 0.01x: flagged at any sane floor, missed by none.
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 1)]
        part += [("ref_ranges_001", i, 50) for i in range(1, 10)]
        self._patch(part, {i: 24 for i in range(1, 10)})
        return {"vet": set(range(1, 10)), "ref_ranges": set(range(1, 10)),
                "sample_chromosome_ploidy": set(range(1, 10))}

    def test_disabled_screen_does_not_flag_a_truncated_sample(self):
        exp = self._truncated_callset()
        r = _run(exp, vet_truncation_threshold=TRUNCATION_SCREEN_DISABLED)
        self.assertFalse(r["truncation_flagged"])
        self.assertEqual(r["flagged_samples"], [])
        # Reported as switched off rather than silently absent, so a green result is not mistaken for
        # a screen that looked and found nothing.
        self.assertTrue(r["details"]["truncation_screen"]["vet"]["disabled"])
        self.assertEqual(r["details"]["truncation_screen"]["vet"]["samples_screened"], 0)

    def test_truncation_threshold_moves_the_floor_on_its_own(self):
        # Sample 9 at 70 against a median of 100 is 0.7x: above the default 0.625 floor, below the
        # 1/1.25 = 0.8 floor a looser truncation threshold sets.
        part = [("vet_001", i, 100) for i in range(1, 9)] + [("vet_001", 9, 70)]
        part += [("ref_ranges_001", i, 50) for i in range(1, 10)]
        self._patch(part, {i: 24 for i in range(1, 10)})
        exp = {"vet": set(range(1, 10)), "ref_ranges": set(range(1, 10)),
               "sample_chromosome_ploidy": set(range(1, 10))}

        self.assertFalse(_run(exp)["truncation_flagged"])

        r = _run(exp, vet_truncation_threshold=1.25)
        self.assertTrue(r["truncation_flagged"])
        self.assertEqual(r["flagged_samples"], [9])

    def test_screen_helper_reports_disabled_without_reading_partitions(self):
        part = [("vet_001", 1, 100), ("vet_001", 2, 1)]
        r = assess_truncation_screen(part, "vet", {1, 2}, self.PREFIXES, TRUNCATION_SCREEN_DISABLED)
        self.assertTrue(r["disabled"])
        self.assertEqual(r["flagged_samples"], [])
        self.assertIsNone(r["baseline"])


class TestTruncationThresholdValidation(unittest.TestCase):
    """run_structural_checks rejects a nonsensical truncation ratio before any BigQuery read, but
    accepts the disable sentinel.

    An invalid threshold must raise before get_partition_row_counts is reached, so no patching is
    needed; if the guard were removed the call would instead try to touch BigQuery and error with a
    different type, failing these assertRaisesRegex checks."""

    def _call(self, threshold):
        _run({"vet": {1}}, vet_truncation_threshold=threshold)

    def test_one_rejected(self):
        with self.assertRaisesRegex(ValueError, "vet_truncation_threshold"):
            self._call(1)

    def test_below_one_rejected(self):
        with self.assertRaisesRegex(ValueError, "vet_truncation_threshold"):
            self._call(0.5)

    def test_negative_rejected(self):
        with self.assertRaisesRegex(ValueError, "vet_truncation_threshold"):
            self._call(-1)

    def test_nan_rejected(self):
        with self.assertRaisesRegex(ValueError, "vet_truncation_threshold"):
            self._call(float("nan"))

    def test_infinite_rejected(self):
        with self.assertRaisesRegex(ValueError, "vet_truncation_threshold"):
            self._call(float("inf"))

    def test_zero_accepted_as_the_disable_sentinel(self):
        # 0 is meaningful here: it switches the screen off, so it must get past validation.
        with patch("verify_structural_checks.get_partition_row_counts", return_value=[]), \
             patch("verify_structural_checks.get_ploidy_row_counts", return_value={}), \
             patch("verify_structural_checks.get_parquet_footer_row_counts", return_value={}):
            r = _run({"vet": {1}}, vet_truncation_threshold=TRUNCATION_SCREEN_DISABLED)
        self.assertFalse(r["truncation_flagged"])



if __name__ == "__main__":
    unittest.main()
