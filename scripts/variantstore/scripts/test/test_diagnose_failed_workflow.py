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

    def test_nonfailed_workflow_rejected(self):
        with self.assertRaisesRegex(ValueError, "Expected a failed workflow"):
            diagnose.report({"status": "Running"})

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
