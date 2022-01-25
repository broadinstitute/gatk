WF_ID=$1
TOKEN=$(gcloud auth print-access-token)

PARAMS="includeKey=start&includeKey=end&includeKey=executionStatus&includeKey=preemptible&includeKey=shardIndex&includeKey=wdl-task-name&includeKey=jes&includeKey=jobId&includeKey=executionEvents&includeKey=outputs"

curl -sL -X 'GET' \
  "https://api.firecloud.org/api/workflows/v1/${WF_ID}/metadata?${PARAMS}" \
  -H 'accept: application/json' \
  -H "Authorization: Bearer ${TOKEN}"
