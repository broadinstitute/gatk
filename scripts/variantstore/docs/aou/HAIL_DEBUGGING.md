# Debugging Issues with Hail/Clusters in WDL

These guidelines assume that you are attempting to do something related to GVS using a WDL to run some Hail in a cluster on Terra and it has failed in a way that leads you to believe something went wrong with either the cluster or Hail (or both?).

## Get Cluster Log Files

If the cluster has already been deleted, most diagnostic information has been lost.  Run the WDL in question with the `leave_hail_cluster_running_at_end` set to `true`,  and then, once the failure has occurred, make note of the cluster ID (should be something like `vds-cluster-52e55963-ff46`), region (probably `us-central1`) and workspace Google project ID and run:

```
gcloud dataproc clusters diagnose <cluster ID> --region=<region> --project=<project ID>
```

This might take a minute or two and will produce an archive of a ton of cluster logs.  Once you've copied this someplace else, feel free to delete the custer using your `@firecloud.org` account in the Google Cloud Console at https://console.cloud.google.com/dataproc/clusters for the workspace Google project ID.
