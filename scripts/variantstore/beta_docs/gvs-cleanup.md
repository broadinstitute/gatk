# Genomic Variant Store Cleanup

Running the GVS leaves behind intermediate files, BigQuery tables of variant data, and the final variant outputs. These data accrue approximately $X/genome/month in storage costs. Here, we introduce steps to clean up those outputs to reduce unnecessary storage cost if you no longer need them.

## Use Case 1: One-Off Callset
### Assumptions
- the samples in the callset are going to be joint called once (no sub-cohorts, no re-using the model)
- the outputs (VCFs, indexes, interval lists, manifest, sample name list) have all been copied to an independent location (not in the workspace bucket)
- no need to maintain provenance
- the workspace where the pipeline was run is only being used for this GVS callset

### Actions
#### Delete BigQuery Dataset
**CAUTION -- This will result in permanent deletion of data.**

In the Google console, within the Explorer panel, select the project where you created the dataset.  Expand the Actions menu by clicking on the three vertical dots next to the dataset name and click Delete.  A dialog will pop up confirming that you want to delete the dataset.  Double-check that you have selected the right dataset to delete before typing "delete" and clicking the button to confirm.
#### Delete Workspace
**CAUTION -- This will result in permanent deletion of files.**

Navigate to the workspace you used to create the callset by going to https://app.terra.bio/#workspaces and selecting the workspace.  While on the Workspace Dashboard page, click on the three vertical dots at the top right of the window and select the "Delete" option.  Double-check that this is what you want to do before typing "Delete Workspace" and clicking the button to confirm.
