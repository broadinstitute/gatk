#!/usr/bin/env python3
"""Diagnose a failed Terra workflow using expanded Cromwell metadata.

Uses Python's standard library and the active gcloud identity. Does not modify
the workflow or download inputs/outputs. See diagnose_failed_workflow.md.
"""

import argparse
import json
import re
import subprocess
import sys
from collections import defaultdict
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlsplit
from urllib.request import Request, urlopen


API = "https://api.firecloud.org/api"
UUID = r"[0-9a-fA-F]{8}(?:-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}"


def execution_ids(uri):
    """Extract the top-level submission/workflow, even from a nested call path."""
    parsed = urlsplit(uri)
    match = re.match(
        rf"^/submissions/({UUID})/[^/]+/({UUID})(?:/|$)", parsed.path
    )
    if parsed.scheme != "gs" or not parsed.netloc or not match:
        raise ValueError(
            "Expected gs://BUCKET/submissions/SUBMISSION_UUID/"
            "WORKFLOW_NAME/WORKFLOW_UUID[/...]"
        )
    return parsed.netloc, match[1], match[2]


class Terra:
    def __init__(self):
        self.token = subprocess.run(
            ["gcloud", "auth", "print-access-token"],
            check=True, capture_output=True, text=True,
        ).stdout.strip()
        if not self.token:
            raise ValueError("gcloud returned an empty access token")

    def request(self, url, headers=None):
        return urlopen(
            Request(url, headers={
                "Authorization": "Bearer " + self.token, **(headers or {})
            }),
            timeout=120,
        )

    def get(self, path):
        with self.request(API + path) as response:
            return json.load(response)

    def workspace(self, bucket):
        matches = [
            entry["workspace"] for entry in self.get("/workspaces")
            if entry["workspace"].get("bucketName") == bucket
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Found {len(matches)} accessible workspaces for bucket {bucket}; "
                "supply --workspace NAMESPACE/NAME"
            )
        return matches[0]["namespace"], matches[0]["name"]

    def metadata(self, workspace, submission, workflow):
        namespace, name = workspace
        return self.get(
            f"/workspaces/{quote(namespace, safe='')}/{quote(name, safe='')}"
            f"/submissions/{submission}/workflows/{workflow}"
            "?expandSubWorkflows=true"
        )

    def log_tail(self, uri, limit):
        parsed = urlsplit(uri)
        if parsed.scheme != "gs" or not parsed.netloc or not parsed.path.strip("/"):
            raise ValueError(f"Not a GCS log URI: {uri}")
        url = (
            "https://storage.googleapis.com/storage/v1/b/"
            f"{quote(parsed.netloc, safe='')}/o/"
            f"{quote(parsed.path[1:], safe='')}?alt=media"
        )
        try:
            with self.request(url, {"Range": f"bytes=-{limit}"}) as response:
                # Never silently report the beginning of an oversized log as its tail.
                data = response.read(limit + 1)
                if len(data) > limit:
                    raise ValueError(f"Server ignored bounded log range for {uri}")
        except HTTPError as error:
            if error.code == 416:
                return "(empty log)"
            raise
        return data.decode("utf-8", errors="replace").strip() or "(empty log)"


def messages(failures, indent=""):
    for failure in failures:
        yield indent + failure.get("message", "(no failure message)")
        yield from messages(failure.get("causedBy", []), indent + "  caused by: ")


def collect_failures(metadata, parents=()):
    """Return final failed calls; retain parent scatter context and retry outcomes."""
    failures = []
    recovered = 0
    for name, attempts in metadata.get("calls", {}).items():
        shards = defaultdict(list)
        for attempt in attempts:
            shards[attempt.get("shardIndex", -1)].append(attempt)
        for shard, history in shards.items():
            final = max(history, key=lambda call: call.get("attempt", 1))
            label = f"{name} [shard={shard if shard != -1 else 'unscattered'}, attempt={final.get('attempt', 1)}]"
            path = parents + (label,)
            recovered += sum(
                call is not final
                and call.get("executionStatus") in ("Failed", "RetryableFailure", "Preempted")
                for call in history
            )
            child = final.get("subWorkflowMetadata")
            if child is not None:
                nested, retries = collect_failures(child, path)
                failures.extend(nested)
                recovered += retries
                if final.get("executionStatus") == "Failed" and not nested:
                    # Expression/input evaluation can fail without any failed task.
                    failures.append((path, final, child.get("failures", [])))
            elif final.get("executionStatus") in ("Failed", "RetryableFailure"):
                failures.append((path, final, []))
    return failures, recovered


def report(metadata, client=None, tail_bytes=65536, tail_lines=40):
    if metadata.get("status") != "Failed":
        raise ValueError(
            f"Expected a failed workflow, got status {metadata.get('status')!r}. "
            "This tool deliberately diagnoses failed runs only."
        )
    print(f"Workflow: {metadata.get('workflowName', '(unknown)')} ({metadata.get('id')})")
    print("Status: Failed")
    print(f"Execution directory: {metadata.get('workflowRoot', '(not recorded)')}")
    print("\nWorkflow failure chain:")
    for message in messages(metadata.get("failures", [])):
        print("  " + message)
    failed, recovered = collect_failures(metadata)
    print(f"\nFinal failed calls: {len(failed)}")
    print(f"Earlier failed/preempted attempts superseded by retries: {recovered}")
    incomplete = False
    for index, (path, call, child_failures) in enumerate(failed, 1):
        print(f"\nFailure {index}: " + " -> ".join(path))
        print(f"  Status: {call.get('executionStatus')}")
        print(f"  Return code: {call.get('returnCode', '(not recorded)')}")
        print(f"  Call directory: {call.get('callRoot', '(not recorded)')}")
        for message in messages(call.get("failures", []) + child_failures):
            print("  " + message)
        if call.get("subWorkflowId"):
            incomplete = True
            print("  WARNING: unexpanded subworkflow; child task diagnosis is incomplete.")
        job = call.get("jobId")
        if job:
            print(f"  Backend job: {job}")
            batch = re.fullmatch(r"projects/([^/]+)/locations/[^/]+/jobs/([^/]+)", job)
            if batch:
                print(
                    "  Batch log command: gcloud logging read "
                    f"""'logName="projects/{batch[1]}/logs/batch_task_logs" """
                    f"""AND resource.labels.job_id="{batch[2]}"' """
                    f"--project {batch[1]} --freshness=365d --format='value(textPayload)'"
                )
        for stream in ("stderr", "stdout"):
            uri = call.get(stream)
            if not uri:
                continue
            print(f"  {stream}: {uri}")
            if client is not None:
                try:
                    tail = client.log_tail(uri, tail_bytes)
                except (HTTPError, URLError, ValueError, TimeoutError) as error:
                    incomplete = True
                    print(f"  WARNING: cannot read {stream}: {error}", file=sys.stderr)
                else:
                    print(f"  Last up to {tail_lines} lines (bounded to {tail_bytes} bytes):")
                    for line in tail.splitlines()[-tail_lines:]:
                        print("    " + line)
    if not failed:
        print(
            "\nNo failed task was recorded. Inspect the workflow failure chain above: "
            "workflow evaluation can fail before a task starts."
        )
    print(
        "\nThese are recorded failures, not an inferred unique root cause. "
        "Empty/missing task logs may require the backend Batch logs; "
        "log tails may omit earlier errors."
    )
    return 2 if incomplete else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--execution-dir", help="Terra gs:// execution directory")
    source.add_argument("--metadata", help="Local expanded Cromwell metadata JSON")
    parser.add_argument("--workspace", help="Terra NAMESPACE/NAME; otherwise resolved from bucket")
    parser.add_argument("--no-log-content", action="store_true", help="Only show log URIs")
    parser.add_argument("--tail-lines", type=int, default=40)
    parser.add_argument("--tail-bytes", type=int, default=65536)
    args = parser.parse_args()
    if args.tail_lines <= 0 or args.tail_bytes <= 0:
        parser.error("--tail-lines and --tail-bytes must be positive")
    try:
        if args.metadata:
            with open(args.metadata) as handle:
                metadata = json.load(handle)
            client = None if args.no_log_content else Terra()
        else:
            bucket, submission, workflow = execution_ids(args.execution_dir)
            client = Terra()
            if args.workspace:
                parts = args.workspace.split("/", 1)
                if len(parts) != 2 or not all(parts):
                    raise ValueError("--workspace must be NAMESPACE/NAME")
                workspace = tuple(parts)
            else:
                workspace = client.workspace(bucket)
            metadata = client.metadata(workspace, submission, workflow)
        return report(
            metadata, None if args.no_log_content else client,
            args.tail_bytes, args.tail_lines,
        )
    except subprocess.CalledProcessError as error:
        print(f"ERROR: gcloud authentication failed: {error.stderr.strip()}", file=sys.stderr)
    except (OSError, ValueError, URLError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
