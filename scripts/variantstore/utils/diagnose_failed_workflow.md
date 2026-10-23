# Diagnose a failed or running Terra workflow

Run `diagnose_failed_workflow.py` with a Terra execution directory:

```bash
python3 scripts/variantstore/utils/diagnose_failed_workflow.py \
  --execution-dir 'gs://WORKSPACE_BUCKET/submissions/SUBMISSION_UUID/WORKFLOW_NAME/WORKFLOW_UUID'
```

Requires Python 3.9+, `gcloud auth login`, Terra access, and permission to read the execution bucket. No additional Python packages or WDL file are required. This is a local utility, not part of the Variants Docker image.

The script resolves the workspace from the bucket and requests expanded Cromwell metadata through Terra. To skip workspace discovery, supply `--workspace 'NAMESPACE/NAME'`. A nested call directory is accepted, but diagnosis covers its enclosing top-level workflow.

Workflows whose recorded status is `Failed` or `Running` are accepted, with the same command-line options. For each reported failed/preempted attempt, the report includes the task name, zero-based shard index (or `unscattered`), attempt number, enclosing subworkflow/shard context, recorded error chain, return code, backend job ID, and stdout/stderr locations. `Failed` attempts are counted separately from `RetryableFailure` and `Preempted`: a retryable failure does not establish that retries were exhausted. Multiple failures are reported without guessing which was the unique original cause. Workflow-level failures are shown even if no task started.

For failed workflows, diagnosis covers the latest attempt of each task/shard; earlier failures superseded by retries are counted separately. For running workflows, all recorded failed/preempted attempts are diagnosed, including failures followed by an active or successful retry. Thus historical failures do not disappear merely because a retry started, but are not necessarily current blockers.

After the workflow header and summary, detailed reports for attempts marked `Failed` appear first, ahead of retryable/preempted attempts and the current-activity listing. Ordering within each failure category is preserved.

Running workflows also list the latest nonterminal task attempts, including nested subworkflow/shard context and recorded backend job IDs. Waiting/queued states are shown as such: this is a metadata snapshot, not confirmation that a VM is executing the command. Only when there are no recorded failures and no unexpanded subworkflows does the report say the workflow is running without recorded errors. This is not a guarantee that logs contain no errors. If no active attempts are visible, the report says so rather than guessing which task is running.

By default, the report reads the last 64 KiB of each failed call's stdout and stderr and prints at most 40 lines per stream. Adjust with `--tail-bytes` and `--tail-lines`; these are tails, not a complete search for errors. Use `--no-log-content` to show only metadata and log URIs. Empty logs and infrastructure failures may require the printed GCP Batch log command; log retention and logging permissions still apply. Dataproc/Hail driver logs may require separate inspection through the task's recorded logs.

For already downloaded metadata, use `--metadata expanded-metadata.json --no-log-content` to work entirely offline. Metadata must have been fetched with `expandSubWorkflows=true`; unexpanded failed subworkflows are explicitly flagged as incomplete.

Exit code 0 means the report completed, not that the workflow succeeded. Exit code 2 means authentication, metadata retrieval, validation, log access, or subworkflow expansion was incomplete. No cloud resources are modified. Reports can contain sensitive error/log data: do not commit them without review and redaction.
