
We have found two tools that are useful in visualizing BigQuery performance:
1. [BQ Visualizer](#bq-visualizer)
2. [Shuffle usage with a Terra Notebook](#shuffle-usage-with-a-terra-notebook)

## BQ Visualizer

### Identify the Query
The first job is to find the query in question and get the job id (e.g. `spec-ops-aou:5b673c2f-8073-4c22-b5bb-b567eceda62f`).  The BigQuery console is very helpful here, using the `Query History` tab to find your query.

This part is painful for AoU because the BigQuery jobs are run as the service account user and our pmi-ops.org accounts don't have access to see those jobs.  This means we can't use the GUI to look for these job ids.

You need to use the BQ command line tool to list jobs

```
bq --project_id <project_id> ls -j -n <number-of-jobs-to-show>
```

We have added several labels to our queries such as:

	gvs_query_name:populate-final-export-pet
	gvs_tool_name:gvs_prepare_callset
	id:aa9b6910

But the bq cli doesn't support filtering jobs by label (feature request since 2017, not yet implemented)

However, we can do our own filtering.  First, we capture the full JSON job details for a large number of jobs.  Enough to be confident that the query we want is in there.  10,000 is a fairly large number for these purposes.  We can capture that to a file, and then process it with JQ to find the one we're interested in

```
# capture 10,000 jobs to all.json
bq --project_id aou-genomics-curation-prod ls -j -n 100000 --format prettyjson > all.json

# extract several fields from the JSON where the job has the label of the gvs_prepare_callset for gvs_tool_name
cat all.json | jq -r '.[] | select(.configuration.labels.gvs_tool_name == "gvs_prepare_callset") | [ .id, .configuration.query.destinationTable.tableId, .configuration.labels.id, .configuration.labels.gvs_query_name, .errorResult.message ] | @csv '
```

### Download Job Metadata

Given a job id, downloading the job information is fairly easy:

```
JOB_ID=spec-ops-aou:5b673c2f-8073-4c22-b5bb-b567eceda62f
bq show --format prettyjson -j $JOB_ID > ${JOB_ID}.json
```

### Visualize w/ BQ Visualizer

Visit https://bqvisualiser.appspot.com/, and upload your JSON!  This is run by google professional services, but do not log in with your credentials to get queries since the app would also have access to all the data in your BigQuery projects.  Upload the JSON instead.

## Shuffle usage with a Terra Notebook

Sometimes it's necessary look at the amount of shuffle used for an entire GvsPrepareCallset run.  One can extract the metadata by filtering by the run id which is output in the logs:

```
cat all.json | jq '[.[] | select(.configuration.labels.id == "6b66ffb0")]' > 22k_success_6b66ffb0.json
```

There is a notebook in the [Spec Ops GVS Analysis](https://app.terra.bio/#workspaces/broad-dsp-spec-ops-fc/Spec%20Ops%20GVS%20Analysis) workspace which can be used to plot this data 
