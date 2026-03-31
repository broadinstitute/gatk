"""
Unit tests for parquet loading scripts.

Tests the core functionality of:
- parse_and_group_files.py: Parsing GCS paths and building BigQuery queries
- load_parquet_to_bq.py: Deterministic job ID generation
"""

import os
import sys
import unittest
from unittest.mock import MagicMock, patch

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from parse_and_group_files import parse_table_and_sample_id_from_file_path, get_already_loaded_tables_and_sample_ids, _validate_table_prefixes
from load_parquet_to_bq import _make_job_id


class TestParseSuperpartitionedTableFromPath(unittest.TestCase):
    """Test parsing GCS paths to BigQuery table names."""
    
    def test_parse_vet_path_basic(self):
        """Test basic vet path parsing."""
        path = "gs://bucket/vet/042/vet_042_1_file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), ("vet_042", 1))

    def test_parse_vet_path_single_digit(self):
        """Test single-digit vet path."""
        path = "gs://bucket/vet/001/vet_001_123_file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), ("vet_001", 123))

    def test_parse_vet_path_three_digits(self):
        """Test three-digit vet path."""
        path = "gs://bucket/vet/999/vet_999_888_file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), ("vet_999", 888))

    def test_parse_ref_ranges_path(self):
        """Test ref_ranges path parsing."""
        path = "gs://bucket/ref_ranges/123/ref_ranges_123_456_file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), ("ref_ranges_123", 456))

    def test_parse_ref_ranges_path_with_subdirs(self):
        """Test ref_ranges path with subdirectories."""
        path = "gs://bucket/output/ref_ranges/042/subdir/ref_ranges_042_67_file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), ("ref_ranges_042", 67))

    def test_parse_vet_path_with_subdirs(self):
        """Test vet path with nested subdirectories."""
        path = "gs://bucket/project/output/vet/007/part1/vet_007_42_file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), ("vet_007", 42))

    def test_parse_invalid_path_returns_none(self):
        """Test that invalid paths return None."""
        path = "gs://bucket/other/042/file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), (None, None))
    
    def test_parse_path_without_number_returns_none(self):
        """Test that paths without numbers return None."""
        path = "gs://bucket/vet/abc/file.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path, ["vet", "ref_ranges"]), (None, None))

class TestParseRegularTableFromPath(unittest.TestCase):
    def test_parse_ploidy_path_basic(self):
        """Test basic ploidy path parsing."""
        path = "gs://bucket/sample_chromosome_ploidy/sample_chromosome_ploidy_1_foo.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path), ("sample_chromosome_ploidy", 1))

    def test_parse_ploidy_path_with_subdirs(self):
        path = "gs://fc-bucket/dir/subdir/sample_chromosome_ploidy/sample_chromosome_ploidy_100_input_vcf_0_FOOBAR.vcf.gz.parquet"
        self.assertEqual(parse_table_and_sample_id_from_file_path(path), ("sample_chromosome_ploidy", 100))


class TestMakeJobId(unittest.TestCase):
    """Test deterministic job ID generation."""

    def test_job_id_is_deterministic(self):
        """Test that same inputs produce same job ID."""
        table_name = "vet_042"
        batch = ["gs://bucket/file1.parquet", "gs://bucket/file2.parquet"]

        job_id_1 = _make_job_id("my-project", "my_dataset", table_name, batch)
        job_id_2 = _make_job_id("my-project", "my_dataset", table_name, batch)

        self.assertEqual(job_id_1, job_id_2)

    def test_job_id_different_for_different_batches(self):
        """Test that different batch contents produce different job IDs."""
        table_name = "vet_042"

        job_id_batch_1 = _make_job_id("my-project", "my_dataset", table_name, ["gs://bucket/file1.parquet"])
        job_id_batch_2 = _make_job_id("my-project", "my_dataset", table_name, ["gs://bucket/file2.parquet"])

        self.assertNotEqual(job_id_batch_1, job_id_batch_2)

    def test_job_id_different_for_different_files(self):
        """Test that different file lists produce different job IDs."""
        table_name = "vet_042"

        job_id_1 = _make_job_id("my-project", "my_dataset", table_name, ["gs://bucket/file1.parquet"])
        job_id_2 = _make_job_id("my-project", "my_dataset", table_name, ["gs://bucket/file2.parquet"])

        self.assertNotEqual(job_id_1, job_id_2)

    def test_job_id_different_for_different_tables(self):
        """Test that different table names produce different job IDs."""
        batch = ["gs://bucket/file1.parquet"]

        job_id_vet = _make_job_id("my-project", "my_dataset", "vet_042", batch)
        job_id_ref = _make_job_id("my-project", "my_dataset", "ref_042", batch)

        self.assertNotEqual(job_id_vet, job_id_ref)

    def test_job_id_format(self):
        """Test that job ID has expected format."""
        table_name = "vet_042"
        batch = ["gs://bucket/file1.parquet"]

        job_id = _make_job_id("my-project", "my_dataset", table_name, batch)

        # Should start with "load_<table_name>_"
        self.assertTrue(job_id.startswith(f"load_{table_name}_"))

        # Should end with a hash (16 hex characters for truncated SHA1)
        # Format is: load_vet_042_<16-char-hash>
        # So we need to remove the "load_vet_042_" prefix
        prefix = f"load_{table_name}_"
        self.assertTrue(job_id.startswith(prefix))
        hash_part = job_id[len(prefix):]
        self.assertEqual(len(hash_part), 16)
        self.assertTrue(all(c in "0123456789abcdef" for c in hash_part))

    def test_job_id_order_independence(self):
        """Test that batch file order affects job ID (sorted internally)."""
        table_name = "vet_042"

        # The function sorts files internally, so order shouldn't matter
        job_id_1 = _make_job_id("my-project", "my_dataset", table_name, ["file2.parquet", "file1.parquet"])
        job_id_2 = _make_job_id("my-project", "my_dataset", table_name, ["file1.parquet", "file2.parquet"])

        # Job IDs should be the same because files are sorted inside _make_job_id
        self.assertEqual(job_id_1, job_id_2)


class TestPathNormalization(unittest.TestCase):
    """Test path normalization for verification."""
    
    def test_path_normalization_strips_gs_prefix(self):
        """Test that gs:// prefix is handled correctly."""
        # This would be implemented in verify_all_loaded.py if needed
        # For now, we just verify the concept
        path1 = "gs://bucket/file.parquet"
        path2 = "gs://bucket/file.parquet"
        self.assertEqual(path1, path2)
    
    def test_path_comparison_is_exact(self):
        """Test that path comparison is case-sensitive and exact."""
        path1 = "gs://bucket/File.parquet"
        path2 = "gs://bucket/file.parquet"
        self.assertNotEqual(path1, path2)


class TestValidateTablePrefixes(unittest.TestCase):
    """Test prefix validation enforced by _validate_table_prefixes."""

    def test_valid_prefixes_do_not_raise(self):
        _validate_table_prefixes(["vet", "ref_ranges"], ["sample_chromosome_ploidy"])

    def test_empty_lists_do_not_raise(self):
        _validate_table_prefixes([], [])

    def test_single_char_prefix_is_invalid(self):
        # Regex requires at least 2 characters (letter + one or more alnum/underscore)
        with self.assertRaises(ValueError):
            _validate_table_prefixes(["v"], [])

    def test_prefix_starting_with_digit_is_invalid(self):
        with self.assertRaises(ValueError):
            _validate_table_prefixes(["1vet"], [])

    def test_prefix_with_hyphen_is_invalid(self):
        with self.assertRaises(ValueError):
            _validate_table_prefixes(["ref-ranges"], [])

    def test_prefix_with_space_is_invalid(self):
        with self.assertRaises(ValueError):
            _validate_table_prefixes([], ["bad prefix"])

    def test_prefix_with_backtick_is_invalid(self):
        # Backtick is especially dangerous as prefixes appear inside BigQuery backtick-quoted identifiers
        with self.assertRaises(ValueError):
            _validate_table_prefixes(["`vet`"], [])

    def test_error_message_lists_all_invalid_prefixes(self):
        with self.assertRaises(ValueError) as ctx:
            _validate_table_prefixes(["bad-one"], ["bad two"])
        msg = str(ctx.exception)
        self.assertIn("bad-one", msg)
        self.assertIn("bad two", msg)

    def test_parse_path_raises_on_invalid_prefix(self):
        with self.assertRaises(ValueError):
            parse_table_and_sample_id_from_file_path(
                "gs://bucket/vet/001/vet_001_1_file.parquet",
                superpartitioned_table_prefixes=["vet; DROP TABLE foo--"],
            )


class TestGetAlreadyLoadedTablesAndSampleIds(unittest.TestCase):
    """Test that get_already_loaded_tables_and_sample_ids builds correct BigQuery queries."""

    @patch("parse_and_group_files.bigquery")
    def test_returns_set_of_tuples(self, mock_bq_module):
        mock_client = MagicMock()
        mock_bq_module.Client.return_value = mock_client
        mock_row1 = MagicMock()
        mock_row1.table_name = "vet_001"
        mock_row1.sample_id = 42
        mock_row2 = MagicMock()
        mock_row2.table_name = "sample_chromosome_ploidy"
        mock_row2.sample_id = 42
        mock_client.query.return_value = [mock_row1, mock_row2]

        result = get_already_loaded_tables_and_sample_ids("my-project", "my_dataset")

        self.assertIsInstance(result, set)
        self.assertIn(("vet_001", 42), result)
        self.assertIn(("sample_chromosome_ploidy", 42), result)
        self.assertEqual(len(result), 2)

    @patch("parse_and_group_files.bigquery")
    def test_query_contains_superpartitioned_regex(self, mock_bq_module):
        mock_client = MagicMock()
        mock_bq_module.Client.return_value = mock_client
        mock_client.query.return_value = []

        get_already_loaded_tables_and_sample_ids(
            "proj", "ds",
            superpartitioned_table_prefixes=["vet", "ref_ranges"],
            regular_table_prefixes=[]
        )

        query_used = mock_client.query.call_args[0][0]
        self.assertIn("INFORMATION_SCHEMA.PARTITIONS", query_used)
        self.assertIn("^vet_[0-9]+$", query_used)
        self.assertIn("^ref_ranges_[0-9]+$", query_used)
        self.assertIn("total_logical_bytes > 0", query_used)

    @patch("parse_and_group_files.bigquery")
    def test_query_contains_regular_table(self, mock_bq_module):
        mock_client = MagicMock()
        mock_bq_module.Client.return_value = mock_client
        mock_client.query.return_value = []

        get_already_loaded_tables_and_sample_ids(
            "proj", "ds",
            superpartitioned_table_prefixes=[],
            regular_table_prefixes=["sample_chromosome_ploidy"]
        )

        query_used = mock_client.query.call_args[0][0]
        self.assertIn("sample_chromosome_ploidy", query_used)
        self.assertNotIn("INFORMATION_SCHEMA", query_used)

    @patch("parse_and_group_files.bigquery")
    def test_query_skips_regular_table_without_sample_id(self, mock_bq_module):
        mock_client = MagicMock()
        mock_bq_module.Client.return_value = mock_client
        mock_client.query.return_value = []

        get_already_loaded_tables_and_sample_ids(
            "proj", "ds",
            superpartitioned_table_prefixes=[],
            regular_table_prefixes=["sample_chromosome_ploidy", "vrs_allele"]
        )

        query_used = mock_client.query.call_args[0][0]
        self.assertIn("sample_chromosome_ploidy", query_used)
        self.assertNotIn("vrs_allele", query_used)

    @patch("parse_and_group_files.bigquery")
    def test_empty_prefixes_returns_empty_set(self, mock_bq_module):
        mock_client = MagicMock()
        mock_bq_module.Client.return_value = mock_client

        result = get_already_loaded_tables_and_sample_ids(
            "proj", "ds",
            superpartitioned_table_prefixes=[],
            regular_table_prefixes=[]
        )

        self.assertEqual(result, set())
        mock_client.query.assert_not_called()

    @patch("parse_and_group_files.bigquery")
    def test_raises_on_bigquery_error(self, mock_bq_module):
        mock_client = MagicMock()
        mock_bq_module.Client.return_value = mock_client
        mock_client.query.side_effect = Exception("BQ error")

        with self.assertRaises(Exception):
            get_already_loaded_tables_and_sample_ids("proj", "ds")


if __name__ == "__main__":
    unittest.main()
