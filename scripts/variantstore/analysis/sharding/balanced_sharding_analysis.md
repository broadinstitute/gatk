## Data Preparation

In order to perform the analysis we generate several input files:

1. `${PREFIX}.task.info.txt` - contains task level information from a run of the `GvsExtractCallset.wdl` pipeline
2. `${PREFIX}.intervals.summary.txt` - Sharded intervals summary for this pipeline run
3. `${PREFIX}.op.details.txt` - Google Pipelines API operation info for the above tasks, for start/stop times

The below assumes `WFID` contains the workflow id of the Extract run in Terra, and `PREFIX` is a unique file prefix

There were a few important runs we analyzed:

| **Prefix**| **Description** | **Workflow ID** | 
| --- | --- | ----------- |
| 40k.original | 40k Ref Ranges Scaling Run | d6c6bb84-fdbb-40bb-8935-cf77d899bb73 |
| 10k.original | 10k Original Run | 259f1102-1c5e-4eda-95ac-f446c65d2d30 |
| 10k.balanced | 10k Balanced Run | 6f7f2e7f-a0bd-4119-ad41-8f8f918654a3 |


### Generating Input Files

First gather the Cromwell Metadata, and generate the task info
```
./get_workflow_metadata.sh ${WFID} > ${PREFIX}.md
./extract_task_info.sh ${PREFIX}.md > ${PREFIX}.task.info.txt
```

To generate the intervals information
```
mkdir ${PREFIX}_interval_files
./extract_intervals_info.sh ${PREFIX}.md > ${PREFIX}.interval_files.list

# copy all the interval files to a temp directory
cat ${PREFIX}.interval_files.list | gsutil -m cp -I ${PREFIX}_interval_files

# summarize the intervals
./summarize_intervals.sh ${PREFIX}_interval_files > ${PREFIX}.intervals.summary.txt

# clean up
rm  ${PREFIX}_interval_files/*.interval_list
rmdir  ${PREFIX}_interval_files
```

To generate the Operations Id information

**NOTE**: this is only necessary because if Cromwell restarts during the workflow run the task start/stop times in the metadata are incorrect.
**NOTE**: you must have your application default credentials authenticated as a user with access to the PAPI operations API

```
cat ${PREFIX}.task.info.txt | cut -f7 | python get_op_info.py > ${PREFIX}.op.details.txt
```

Typically, these resulting files will be pushed to a Google Bucket, referred to as `BUCKET_PATH` below

```
BUCKET_PATH=gs://fc-80261891-d43b-4aeb-b868-3b7489d0ebf9/gvs/analysis2/

gsutil cp *op.details.txt ${BUCKET_PATH}
gsutil cp *task.info.txt ${BUCKET_PATH}
gsutil cp *intervals.summary.txt ${BUCKET_PATH}
```

## Visualizations

See the "Balanced Sharding.ipynb" Jupyter notebook for plots


## Modeling




