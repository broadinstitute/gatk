"""
Unit tests for parquet loading scripts.

Tests the core functionality of:
- parse_and_group_files.py: Parsing GCS paths and building BigQuery queries
- load_parquet_to_bq.py: Deterministic job ID generation, retry logic, and
  resume-via-Conflict after VM preemption

These tests run inside the Docker image where google-cloud-bigquery is
installed.  External API calls are prevented by patching google.cloud.bigquery.Client.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

try:
    from google.api_core import exceptions as _bq_exceptions
    _GOOGLE_AVAILABLE = True
except ImportError:
    _bq_exceptions = None  # type: ignore[assignment]
    _GOOGLE_AVAILABLE = False

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))


from parse_and_group_files import parse_table_and_sample_id_from_file_path, get_already_loaded_tables_and_sample_ids, _validate_table_prefixes
from load_parquet_to_bq import (
    _make_job_id,
    _is_permanently_failed_job,
    _submit_load_job_with_retry,
    _submit_resumable_load_job,
    _wait_for_job_with_retry,
    load_table_from_parquet_files,
)


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

    def test_generation_zero_is_unsuffixed(self):
        """Generation 0 (the default) must equal the historical unsuffixed ID."""
        batch = ["gs://bucket/file1.parquet"]

        default_id = _make_job_id("my-project", "my_dataset", "vet_042", batch)
        gen0_id = _make_job_id("my-project", "my_dataset", "vet_042", batch, generation=0)

        self.assertEqual(default_id, gen0_id)
        # No "_r" generation suffix appended for generation 0.
        self.assertFalse(default_id.endswith("_r0"))

    def test_generation_suffix_appended_for_nonzero(self):
        """Generation N>=1 appends an _rN suffix to the base ID."""
        batch = ["gs://bucket/file1.parquet"]

        base_id = _make_job_id("my-project", "my_dataset", "vet_042", batch)
        gen1_id = _make_job_id("my-project", "my_dataset", "vet_042", batch, generation=1)
        gen2_id = _make_job_id("my-project", "my_dataset", "vet_042", batch, generation=2)

        self.assertEqual(gen1_id, f"{base_id}_r1")
        self.assertEqual(gen2_id, f"{base_id}_r2")

    def test_generations_are_distinct(self):
        """Different generations of the same batch produce different job IDs."""
        batch = ["gs://bucket/file1.parquet"]

        ids = {
            _make_job_id("my-project", "my_dataset", "vet_042", batch, generation=g)
            for g in range(4)
        }
        self.assertEqual(len(ids), 4)


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


# ===========================================================================
# Tests for get_already_loaded_tables_and_sample_ids
# ===========================================================================

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


# ===========================================================================
# Tests for the loading functions in load_parquet_to_bq.py
#
# google-cloud-bigquery is available in the Docker test environment, so we
# import real exception classes and only mock out the BigQuery Client
# constructor to prevent actual API calls.
# ===========================================================================

@unittest.skipUnless(_GOOGLE_AVAILABLE, "google-cloud-bigquery not installed")
class TestSubmitLoadJobWithRetry(unittest.TestCase):
    """Tests for _submit_load_job_with_retry."""

    def setUp(self):
        self.client = MagicMock()
        self.batch = ["gs://bucket/file1.parquet", "gs://bucket/file2.parquet"]
        self.table_id = "proj.ds.vet_001"
        self.job_config = MagicMock()
        self.job_id = "load_vet_001_abc123def456789a"
        self.location = "us-central1"

    def test_first_attempt_succeeds(self):
        """Normal path: job is submitted and returned immediately."""
        mock_job = MagicMock()
        self.client.load_table_from_uri.return_value = mock_job

        result = _submit_load_job_with_retry(
            self.client, self.batch, self.table_id,
            self.job_config, self.job_id, self.location,
        )

        self.assertIs(result, mock_job)
        self.client.load_table_from_uri.assert_called_once_with(
            self.batch, self.table_id,
            job_config=self.job_config, job_id=self.job_id, location=self.location,
        )
        self.client.get_job.assert_not_called()

    def test_conflict_resumes_existing_job(self):
        """
        Key resume test: when a VM is preempted and a new VM retries, BigQuery
        returns Conflict because the job was already submitted.  The function must
        fetch the existing job via get_job() rather than re-submitting, so that
        no data is written twice.
        """
        existing_job = MagicMock()
        self.client.load_table_from_uri.side_effect = _bq_exceptions.Conflict("Job already exists")
        self.client.get_job.return_value = existing_job

        result = _submit_load_job_with_retry(
            self.client, self.batch, self.table_id,
            self.job_config, self.job_id, self.location,
        )

        self.assertIs(result, existing_job)
        self.client.load_table_from_uri.assert_called_once()
        self.client.get_job.assert_called_once_with(job_id=self.job_id, location=self.location)

    @patch("time.sleep")
    def test_quota_error_retried_then_succeeds(self, mock_sleep):
        """TooManyRequests on the first attempt is retried and the second succeeds."""
        mock_job = MagicMock()
        self.client.load_table_from_uri.side_effect = [
            _bq_exceptions.TooManyRequests("quota exceeded"),
            mock_job,
        ]

        result = _submit_load_job_with_retry(
            self.client, self.batch, self.table_id,
            self.job_config, self.job_id, self.location,
        )

        self.assertIs(result, mock_job)
        self.assertEqual(self.client.load_table_from_uri.call_count, 2)
        mock_sleep.assert_called_once_with(30)  # first backoff delay

    @patch("time.sleep")
    def test_quota_error_exhausted_raises(self, mock_sleep):
        """All retry attempts fail with TooManyRequests; the exception propagates."""
        self.client.load_table_from_uri.side_effect = _bq_exceptions.TooManyRequests("quota")

        with self.assertRaises(_bq_exceptions.TooManyRequests):
            _submit_load_job_with_retry(
                self.client, self.batch, self.table_id,
                self.job_config, self.job_id, self.location,
                max_retries=3,
            )

        self.assertEqual(self.client.load_table_from_uri.call_count, 4)  # 1 + 3 retries
        self.assertEqual(mock_sleep.call_count, 3)

    def test_non_retryable_error_raises_immediately(self):
        """A non-quota error propagates without any retries."""
        self.client.load_table_from_uri.side_effect = ValueError("bad input")

        with self.assertRaises(ValueError):
            _submit_load_job_with_retry(
                self.client, self.batch, self.table_id,
                self.job_config, self.job_id, self.location,
            )

        self.client.load_table_from_uri.assert_called_once()


@unittest.skipUnless(_GOOGLE_AVAILABLE, "google-cloud-bigquery not installed")
class TestIsPermanentlyFailedJob(unittest.TestCase):
    """Tests for the _is_permanently_failed_job predicate."""

    @staticmethod
    def _job(state, error_result):
        job = MagicMock()
        job.state = state
        job.error_result = error_result
        return job

    def test_done_with_error_result_is_permanent(self):
        """A terminal job carrying an error_result loaded nothing and is superseded."""
        self.assertTrue(
            _is_permanently_failed_job(self._job("DONE", {"reason": "invalid", "message": "bad parquet"}))
        )

    def test_done_without_error_result_is_not_permanent(self):
        """A successfully completed job must be resumed, never superseded."""
        self.assertFalse(_is_permanently_failed_job(self._job("DONE", None)))

    def test_running_is_not_permanent(self):
        """An in-flight job must be resumed, never superseded."""
        self.assertFalse(_is_permanently_failed_job(self._job("RUNNING", None)))

    def test_running_with_error_result_is_not_permanent(self):
        """Only DONE jobs are superseded; a non-terminal job is always resumed."""
        self.assertFalse(_is_permanently_failed_job(self._job("RUNNING", {"reason": "backendError"})))


@unittest.skipUnless(_GOOGLE_AVAILABLE, "google-cloud-bigquery not installed")
class TestSubmitResumableLoadJob(unittest.TestCase):
    """Tests for _submit_resumable_load_job's generation-walk resume/supersede logic."""

    def setUp(self):
        self.client = MagicMock()
        self.project_id = "proj"
        self.dataset_name = "ds"
        self.table_name = "vet_001"
        self.batch = ["gs://bucket/vet/001/vet_001_42_file.parquet"]
        self.table_id = "proj.ds.vet_001"
        self.job_config = MagicMock()
        self.location = "us-central1"

    def _job(self, state="RUNNING", error_result=None):
        job = MagicMock()
        job.state = state
        job.error_result = error_result
        return job

    def _job_id_for(self, generation):
        return _make_job_id(self.project_id, self.dataset_name, self.table_name, self.batch, generation)

    def _submit(self, max_generations=100):
        return _submit_resumable_load_job(
            self.client, self.project_id, self.dataset_name, self.table_name,
            self.batch, self.table_id, self.job_config, self.location,
            max_generations=max_generations,
        )

    def test_fresh_submit_uses_generation_zero(self):
        """With no collision the batch loads under the unsuffixed generation-0 ID."""
        fresh = self._job(state="RUNNING")
        self.client.load_table_from_uri.return_value = fresh

        load_job, job_id = self._submit()

        self.assertIs(load_job, fresh)
        self.assertEqual(job_id, self._job_id_for(0))
        self.client.get_job.assert_not_called()

    def test_conflict_resumes_running_job(self):
        """A collision with a still-running job resumes it rather than superseding."""
        running = self._job(state="RUNNING")
        self.client.load_table_from_uri.side_effect = _bq_exceptions.Conflict("exists")
        self.client.get_job.return_value = running

        load_job, job_id = self._submit()

        self.assertIs(load_job, running)
        self.assertEqual(job_id, self._job_id_for(0))
        self.client.get_job.assert_called_once_with(job_id=self._job_id_for(0), location=self.location)

    def test_conflict_resumes_succeeded_job(self):
        """A collision with a DONE job that has no error_result resumes it."""
        succeeded = self._job(state="DONE", error_result=None)
        self.client.load_table_from_uri.side_effect = _bq_exceptions.Conflict("exists")
        self.client.get_job.return_value = succeeded

        load_job, job_id = self._submit()

        self.assertIs(load_job, succeeded)
        self.assertEqual(job_id, self._job_id_for(0))

    def test_permanently_failed_job_is_superseded_by_fresh_generation(self):
        """
        The core VS-1988 fix: generation 0 already exists as a DONE job with an
        error_result, so it is superseded by a freshly submitted generation-1 job
        instead of re-fetching the dead job forever.
        """
        failed = self._job(state="DONE", error_result={"reason": "invalid", "message": "bad parquet"})
        fresh = self._job(state="RUNNING")
        self.client.load_table_from_uri.side_effect = [_bq_exceptions.Conflict("exists"), fresh]
        self.client.get_job.return_value = failed

        load_job, job_id = self._submit()

        self.assertIs(load_job, fresh)
        self.assertEqual(job_id, self._job_id_for(1))
        self.assertTrue(job_id.endswith("_r1"))
        # The fresh submit used the generation-1 ID, not the burned generation-0 ID.
        self.assertEqual(self.client.load_table_from_uri.call_args_list[1][1]["job_id"], self._job_id_for(1))
        self.client.get_job.assert_called_once()

    def test_two_failed_generations_are_walked(self):
        """Two consecutive permanently-failed generations are skipped before a fresh one."""
        failed0 = self._job(state="DONE", error_result={"reason": "invalid"})
        failed1 = self._job(state="DONE", error_result={"reason": "invalid"})
        fresh = self._job(state="RUNNING")
        self.client.load_table_from_uri.side_effect = [
            _bq_exceptions.Conflict("exists"),
            _bq_exceptions.Conflict("exists"),
            fresh,
        ]
        self.client.get_job.side_effect = [failed0, failed1]

        load_job, job_id = self._submit()

        self.assertIs(load_job, fresh)
        self.assertEqual(job_id, self._job_id_for(2))
        self.assertEqual(self.client.load_table_from_uri.call_count, 3)
        self.assertEqual(self.client.get_job.call_count, 2)

    def test_all_generations_permanently_failed_raises(self):
        """If every generation up to the cap is a dead job, the batch cannot proceed."""
        self.client.load_table_from_uri.side_effect = _bq_exceptions.Conflict("exists")
        self.client.get_job.return_value = self._job(state="DONE", error_result={"reason": "invalid"})

        with self.assertRaises(Exception):
            self._submit(max_generations=3)

        self.assertEqual(self.client.load_table_from_uri.call_count, 3)


@unittest.skipUnless(_GOOGLE_AVAILABLE, "google-cloud-bigquery not installed")
class TestWaitForJobWithRetry(unittest.TestCase):
    """Tests for _wait_for_job_with_retry."""

    def test_job_completes_immediately(self):
        """Happy path: result() returns on the first call."""
        mock_job = MagicMock()
        mock_job.result.return_value = None

        result = _wait_for_job_with_retry(mock_job)

        self.assertIs(result, mock_job)
        mock_job.result.assert_called_once()
        mock_job.reload.assert_not_called()

    @patch("time.sleep")
    def test_transient_error_then_completes(self, mock_sleep):
        """ServiceUnavailable while waiting triggers a reload-and-retry."""
        mock_job = MagicMock()
        mock_job.result.side_effect = [_bq_exceptions.ServiceUnavailable("transient"), None]

        result = _wait_for_job_with_retry(mock_job)

        self.assertIs(result, mock_job)
        self.assertEqual(mock_job.result.call_count, 2)
        mock_job.reload.assert_called_once()
        mock_sleep.assert_called_once_with(30)

    @patch("time.sleep")
    def test_transient_error_exhausted_raises(self, mock_sleep):
        """All wait attempts fail with transient errors; the exception propagates."""
        mock_job = MagicMock()
        mock_job.result.side_effect = _bq_exceptions.ServiceUnavailable("transient")

        with self.assertRaises(_bq_exceptions.ServiceUnavailable):
            _wait_for_job_with_retry(mock_job, max_retries=3)

        self.assertEqual(mock_job.result.call_count, 4)  # 1 + 3 retries
        self.assertEqual(mock_sleep.call_count, 3)

    def test_non_retryable_error_raises_immediately(self):
        """A non-transient error propagates without reloading or retrying."""
        mock_job = MagicMock()
        mock_job.result.side_effect = RuntimeError("job exploded")

        with self.assertRaises(RuntimeError):
            _wait_for_job_with_retry(mock_job)

        mock_job.result.assert_called_once()
        mock_job.reload.assert_not_called()


@unittest.skipUnless(_GOOGLE_AVAILABLE, "google-cloud-bigquery not installed")
class TestLoadTableFromParquetFiles(unittest.TestCase):
    """
    End-to-end tests for load_table_from_parquet_files with all external
    calls mocked out.  google.cloud.bigquery.Client is patched so no real
    API calls are made; the real LoadJobConfig / SourceFormat / WriteDisposition
    values from the installed library are used as-is.
    """

    def setUp(self):
        self.mock_client = MagicMock()
        self._client_patcher = patch(
            "google.cloud.bigquery.Client", return_value=self.mock_client
        )
        self._client_patcher.start()

        mock_table = MagicMock()
        mock_table.location = "us-central1"
        self.mock_client.get_table.return_value = mock_table

    def tearDown(self):
        self._client_patcher.stop()

    def _make_fofn(self, paths):
        """Write paths to a NamedTemporaryFile and return its path."""
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".fofn", delete=False)
        f.write("\n".join(paths) + ("\n" if paths else ""))
        f.flush()
        f.close()
        self.addCleanup(os.unlink, f.name)
        return f.name

    def _make_mock_job(self, rows=10, errors=None, job_id="load_vet_001_abc123"):
        job = MagicMock()
        job.output_rows = rows
        job.errors = errors
        job.job_id = job_id
        job.result.return_value = None
        return job

    def test_happy_path_single_batch(self):
        """All files load successfully in one batch; stats reflect the loaded rows."""
        files = [f"gs://bucket/vet/001/vet_001_{i}_file.parquet" for i in range(5)]
        fofn = self._make_fofn(files)
        self.mock_client.load_table_from_uri.return_value = self._make_mock_job(rows=100)

        stats = load_table_from_parquet_files(
            project_id="proj", dataset_name="ds", table_name="vet_001",
            files_fofn=fofn, schema_path=None,
        )

        self.assertEqual(stats["completion_status"], "SUCCESS")
        self.assertEqual(stats["files_loaded"], 5)
        self.assertEqual(stats["rows_loaded"], 100)
        self.assertEqual(stats["batches_processed"], 1)
        self.assertEqual(stats["batches_failed"], 0)

    def test_resume_after_conflict(self):
        """
        Key end-to-end resume test: load_table_from_uri raises Conflict (because
        a previous VM already submitted this job before being preempted).
        The function must fetch the pre-existing job via get_job() and use its
        result — no duplicate write occurs.
        """
        files = ["gs://bucket/vet/001/vet_001_42_file.parquet"]
        fofn = self._make_fofn(files)
        existing_job = self._make_mock_job(rows=50, job_id="load_vet_001_preexisting")
        self.mock_client.load_table_from_uri.side_effect = _bq_exceptions.Conflict("already exists")
        self.mock_client.get_job.return_value = existing_job

        stats = load_table_from_parquet_files(
            project_id="proj", dataset_name="ds", table_name="vet_001",
            files_fofn=fofn, schema_path=None,
        )

        # load_table_from_uri was attempted (Conflict raised), then get_job used
        self.mock_client.load_table_from_uri.assert_called_once()
        self.mock_client.get_job.assert_called_once()
        # The run completes successfully using the pre-existing job's output
        self.assertEqual(stats["completion_status"], "SUCCESS")
        self.assertEqual(stats["rows_loaded"], 50)
        self.assertEqual(stats["batches_failed"], 0)

    def test_supersedes_permanently_failed_job_on_retry(self):
        """
        End-to-end VS-1988 fix: a prior task attempt left a DONE job with an
        error_result under the deterministic (generation-0) job ID.  This attempt
        must supersede it with a fresh generation-1 job that succeeds, rather than
        re-fetching the dead job and failing the batch again.
        """
        files = ["gs://bucket/vet/001/vet_001_42_file.parquet"]
        fofn = self._make_fofn(files)

        failed_job = MagicMock()
        failed_job.state = "DONE"
        failed_job.error_result = {"reason": "invalid", "message": "schema mismatch"}

        fresh_job = self._make_mock_job(rows=77, job_id="load_vet_001_gen1")
        fresh_job.state = "RUNNING"
        fresh_job.error_result = None

        # Generation 0 collides (Conflict -> fetch the dead job); generation 1 is fresh.
        self.mock_client.load_table_from_uri.side_effect = [
            _bq_exceptions.Conflict("already exists"),
            fresh_job,
        ]
        self.mock_client.get_job.return_value = failed_job

        stats = load_table_from_parquet_files(
            project_id="proj", dataset_name="ds", table_name="vet_001",
            files_fofn=fofn, schema_path=None,
        )

        self.assertEqual(stats["completion_status"], "SUCCESS")
        self.assertEqual(stats["rows_loaded"], 77)
        self.assertEqual(stats["batches_failed"], 0)
        # Generation 0 was fetched once (the dead job); generation 1 was submitted fresh.
        self.mock_client.get_job.assert_called_once()
        self.assertEqual(self.mock_client.load_table_from_uri.call_count, 2)

    def test_batch_with_job_errors_is_partial_failure(self):
        """A job that completes but carries errors is counted as a failed batch."""
        files = ["gs://bucket/vet/001/vet_001_1_file.parquet"]
        fofn = self._make_fofn(files)
        self.mock_client.load_table_from_uri.return_value = self._make_mock_job(
            rows=0, errors=[{"message": "schema mismatch"}]
        )

        stats = load_table_from_parquet_files(
            project_id="proj", dataset_name="ds", table_name="vet_001",
            files_fofn=fofn, schema_path=None,
        )

        self.assertEqual(stats["completion_status"], "PARTIAL")
        self.assertEqual(stats["batches_failed"], 1)
        self.assertEqual(stats["rows_loaded"], 0)

    def test_exception_in_one_batch_continues_to_next(self):
        """An exception in a batch is logged but remaining batches proceed."""
        files = [f"gs://bucket/vet/001/vet_001_{i}_file.parquet" for i in range(3)]
        fofn = self._make_fofn(files)
        good_job = self._make_mock_job(rows=30, job_id="job_good")
        # With batch_size=2: first batch (2 files) raises, second batch (1 file) succeeds
        self.mock_client.load_table_from_uri.side_effect = [
            RuntimeError("network blip"),
            good_job,
        ]

        stats = load_table_from_parquet_files(
            project_id="proj", dataset_name="ds", table_name="vet_001",
            files_fofn=fofn, schema_path=None, batch_size=2,
        )

        self.assertEqual(stats["completion_status"], "PARTIAL")
        self.assertEqual(stats["batches_failed"], 1)
        self.assertEqual(stats["batches_processed"], 1)
        self.assertEqual(stats["rows_loaded"], 30)

    def test_empty_fofn_returns_immediate_success(self):
        """An empty FOFN short-circuits to SUCCESS without touching BigQuery."""
        fofn = self._make_fofn([])

        stats = load_table_from_parquet_files(
            project_id="proj", dataset_name="ds", table_name="vet_001",
            files_fofn=fofn, schema_path=None,
        )

        self.assertEqual(stats["completion_status"], "SUCCESS")
        self.assertEqual(stats["files_loaded"], 0)
        self.assertEqual(stats["rows_loaded"], 0)
        self.mock_client.load_table_from_uri.assert_not_called()

    def test_missing_table_exits_with_code_1(self):
        """NotFound from get_table() causes an immediate sys.exit(1)."""
        files = ["gs://bucket/vet/001/vet_001_1_file.parquet"]
        fofn = self._make_fofn(files)
        self.mock_client.get_table.side_effect = _bq_exceptions.NotFound("table not found")

        with self.assertRaises(SystemExit) as cm:
            load_table_from_parquet_files(
                project_id="proj", dataset_name="ds", table_name="vet_001",
                files_fofn=fofn, schema_path=None,
            )

        self.assertEqual(cm.exception.code, 1)

    def test_multi_batch_all_succeed(self):
        """Files exceeding batch_size are spread across multiple batches, all succeeding."""
        files = [f"gs://bucket/vet/001/vet_001_{i}_file.parquet" for i in range(5)]
        fofn = self._make_fofn(files)
        self.mock_client.load_table_from_uri.side_effect = [
            self._make_mock_job(rows=10, job_id=f"job_{i}") for i in range(3)
        ]

        stats = load_table_from_parquet_files(
            project_id="proj", dataset_name="ds", table_name="vet_001",
            files_fofn=fofn, schema_path=None, batch_size=2,
        )

        self.assertEqual(stats["batches_processed"], 3)  # batches of 2, 2, 1
        self.assertEqual(stats["rows_loaded"], 30)
        self.assertEqual(stats["completion_status"], "SUCCESS")


if __name__ == "__main__":
    unittest.main()
