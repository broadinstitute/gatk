import contextlib
import importlib.util
import io
import os
from pathlib import Path
import unittest
from urllib.error import URLError


VARIANTSTORE = Path(os.environ.get(
    "GVS_VARIANTSTORE_DIR", Path(__file__).resolve().parents[2]))
SCRIPT = VARIANTSTORE / "utils/diagnose_failed_workflow.py"
SPEC = importlib.util.spec_from_file_location("diagnose_failed_workflow", SCRIPT)
diagnose = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(diagnose)


def call(status, shard=-1, attempt=1, **kwargs):
    return dict(executionStatus=status, shardIndex=shard, attempt=attempt, **kwargs)


class DiagnoseFailedWorkflowTest(unittest.TestCase):
    def test_failed_details_precede_retryable_failures_and_activity(self):
        for status in ("Failed", "Running"):
            with self.subTest(status=status):
                metadata = {"status": status, "calls": {
                    "W.Retry": [call("RetryableFailure")],
                    "W.Active": [call("Running")],
                    "W.Child": [call("Failed", subWorkflowMetadata={"calls": {
                        "Child.Preempted": [call("Preempted")],
                        "Child.Fatal": [call("Failed", shard=42,
                                            failures=[{"message": "fatal error detail"}])],
                    }})],
                }}
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    diagnose.report(metadata)
                text = output.getvalue()
                self.assertIn("Failure 1: W.Child", text)
                self.assertLess(text.index("fatal error detail"), text.index("Failure 2:"))
                self.assertLess(text.index("Child.Fatal"), text.index("W.Retry"))
                self.assertLess(text.index("Child.Fatal"), text.index("Child.Preempted"))
                if status == "Running":
                    self.assertLess(text.index("fatal error detail"),
                                    text.index("Current task/subworkflow activity"))

    def test_nested_shards_and_retries(self):
        metadata = {"calls": {"Outer.Child": [
            call("Failed", shard=2, subWorkflowMetadata={"calls": {"Child.Task": [
                call("RetryableFailure", shard=0, attempt=1),
                call("Done", shard=0, attempt=2),
                call("Failed", shard=1, attempt=1),
                call("Failed", shard=1, attempt=2, returnCode=137),
            ]}})
        ]}}
        failures, recovered = diagnose.collect_failures(metadata)
        self.assertEqual(1, len(failures))
        self.assertEqual(2, recovered)
        path, final, _ = failures[0]
        self.assertIn("Outer.Child [shard=2", path[0])
        self.assertIn("Child.Task [shard=1, attempt=2]", path[1])
        self.assertEqual(137, final["returnCode"])

    def test_expression_failure_and_unexpanded_child(self):
        for extra in (
            {"subWorkflowMetadata": {"failures": [{"message": "bad expression"}]}},
            {"subWorkflowId": "child-id"},
        ):
            with self.subTest(extra=extra):
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    result = diagnose.report({
                        "status": "Failed", "calls": {"Outer.Child": [call("Failed", **extra)]}
                    })
                self.assertIn("Failure 1:", output.getvalue())
                if "subWorkflowId" in extra:
                    self.assertEqual(2, result)
                    self.assertIn("incomplete", output.getvalue())
                else:
                    self.assertEqual(0, result)
                    self.assertIn("bad expression", output.getvalue())

    def test_log_failure_is_not_silent(self):
        class Client:
            def log_tail(self, uri, limit):
                raise URLError("permission denied")

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as errors:
            result = diagnose.report({
                "status": "Failed", "calls": {
                    "W.Task": [call("Failed", stderr="gs://bucket/stderr")]
                },
            }, Client())
        self.assertEqual(2, result)
        self.assertIn("permission denied", errors.getvalue())

    def test_message_chain_and_log_tail(self):
        class Client:
            def log_tail(self, uri, limit):
                return "older\nactual error"

        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = diagnose.report({
                "status": "Failed", "calls": {"W.Task": [call(
                    "Failed", stderr="gs://bucket/stderr",
                    failures=[{"message": "outer", "causedBy": [{"message": "inner"}]}],
                )]},
            }, Client(), tail_lines=1)
        self.assertEqual(0, result)
        self.assertIn("caused by: inner", output.getvalue())
        self.assertIn("actual error", output.getvalue())
        self.assertNotIn("older", output.getvalue())

    def test_running_without_errors_lists_nested_latest_activity(self):
        metadata = {"status": "Running", "calls": {"Outer.Child": [
            call("Running", shard=3, subWorkflowMetadata={"calls": {"Child.Task": [
                call("Done", shard=0), call("Running", shard=1, jobId="backend-job"),
                call("QueuedInCromwell", shard=2), call("NotStarted", shard=4),
            ]}}),
        ]}}
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = diagnose.report(metadata)
        text = output.getvalue()
        self.assertEqual(0, result)
        self.assertIn("Running without recorded errors", text)
        self.assertIn("Outer.Child [shard=3", text)
        self.assertIn("Child.Task [shard=1", text)
        self.assertIn("backend-job", text)
        self.assertIn("QueuedInCromwell", text)
        self.assertIn("NotStarted", text)
        self.assertNotIn("Child.Task [shard=0", text)

    def test_running_diagnoses_history_even_after_success_or_retry(self):
        for latest_status in ("Running", "Done"):
            with self.subTest(latest_status=latest_status):
                class Client:
                    def log_tail(self, uri, limit):
                        return "original failure detail"

                metadata = {"status": "Running", "calls": {"W.Task": [
                    call("RetryableFailure", shard=7, stderr="gs://bucket/log",
                         failures=[{"message": "worker died"}]),
                    call(latest_status, shard=7, attempt=2),
                    call("Failed", shard=8, failures=[{"message": "bad input"}]),
                    call("Preempted", shard=9),
                ]}}
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    result = diagnose.report(metadata, Client())
                text = output.getvalue()
                self.assertEqual(0, result)
                self.assertIn("Failed/preempted attempts (all recorded attempts): 3", text)
                self.assertIn("Failed: 1", text)
                self.assertIn("Retryable/preempted: 2", text)
                self.assertIn("worker died", text)
                self.assertIn("original failure detail", text)
                self.assertIn("bad input", text)
                self.assertNotIn("Running without recorded errors", text)
                active, _ = diagnose.collect_activity(metadata)
                self.assertEqual(1 if latest_status == "Running" else 0, len(active))

    def test_running_unexpanded_child_never_reports_error_free(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = diagnose.report({"status": "Running", "calls": {
                "W.Child": [call("Running", subWorkflowId="child-id")],
            }})
        self.assertEqual(2, result)
        self.assertIn("unexpanded subworkflow", output.getvalue())
        self.assertNotIn("Running without recorded errors", output.getvalue())

    def test_running_workflow_error_without_task_failure(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            diagnose.report({"status": "Running", "failures": [{"message": "evaluation failed"}]})
        self.assertIn("evaluation failed", output.getvalue())
        self.assertNotIn("Running without recorded errors", output.getvalue())

    def test_other_workflow_status_rejected(self):
        with self.assertRaisesRegex(ValueError, "Expected a Failed or Running workflow"):
            diagnose.report({"status": "Succeeded"})

    def test_execution_directory(self):
        uuid = "12345678-1234-1234-1234-123456789abc"
        self.assertEqual(
            ("bucket", uuid, uuid),
            diagnose.execution_ids(f"gs://bucket/submissions/{uuid}/W/{uuid}/call-child/"),
        )
        with self.assertRaises(ValueError):
            diagnose.execution_ids("gs://bucket/arbitrary-output-folder")


if __name__ == "__main__":
    unittest.main()
