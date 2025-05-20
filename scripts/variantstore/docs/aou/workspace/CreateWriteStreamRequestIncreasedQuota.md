## Request Additional Quota for Bulk Ingest

### Overview

The GVS bulk ingest process is prone to hitting BigQuery `CreateWriteStream` quota limits when loading large numbers of
VCFs, which certainly applies for AoU callsets. If an ingest job hits this quota it will immediately fail and
exit with an error. To avoid this, we want to calculate the rate at which a single ingest job issues `CreateWriteStream`
requests and limit the scattering of ingest tasks to not exceed that rate.

### Ingest rate

It cannot be assumed that whatever ingest rate was observed in the previous callset will be the same for the next
callset. The performance of ingest jobs can vary from callset to callset due to various factors (code changes,
changes in how jobs are invoked, etc).
For example, in the Echo callset an ingest process took about 10-15 minutes per sample (headers + references +
variants). Foxtrot ingest breaks up ingestion into a header phase at about 3:30 per sample, followed by a reference and
variant phase at under 9 minutes per sample. Before kicking off a full callset it's a good idea to do at least a
small ingest run to measure ingest performance.

Each ingest job takes a group of samples and loads them into BigQuery. In order to calculate the rate at which a single
ingest job issues `CreateWriteStream` requests, we need to know the rate at which the job loads samples into BigQuery.
This information can be found in the logs of the ingest job. For example, looking for one particular line in the stderr
file for an ingest job:

```
$ grep 'gatk --java-options -Xmx2g CreateVariantIngestFiles' stderr | awk '{print $2}'
```

This will show the times at which the ingest process kicked off the actual GATK command to load data into BigQuery.
This corresponds to the same point in the ingest process for each sample, so by differencing consecutive timestamps
(readily done in a spreadsheet), we can see how long it took to load each sample.

### WriteStream Utilization

The code in the GVS ingest process has evolved over time, using different numbers of write streams at different points
in its history. In Delta and earlier two write streams were used per sample, while Echo added a third write stream for
header data. We still intend to load header data for Foxtrot, but now headers are loaded in a separate phase before the
variant and reference data. So for Foxtrot we will estimate write stream consumption based on the more intensive
two-streams-per-sample variant and reference data loading phase of the Foxtrot ingest process.

### Requesting a Quota Increase

When the time comes to kick off ingest, request a quota increase from Google. Unfortunately we don't seem to be able to
ask for this quota increase in advance; during Echo Google only allowed us to run with the quota we specified for a very
limited time. You can point to [this ticket](https://console.cloud.google.com/support/cases/detail/v2/47548796?project=broad-dsde-methods)
as an example. We should ask for our quota of CreateWriteStream requests in us-central1 to be adjusted
to 45k per 4 hours. Mathematically this is about the same as the 167 requests per minute value, but apparently we had
issues with the "per minute" granularity due to the "bursty" way our ingest processes run.

Note that the various CreateWriteStream quotas in the Google Cloud Console don't appear to accurately report actual
usage. In all the visualizations this page offers, usage always appears to be zero. Similarly, this page does not appear
to accurately reflect the quota that is actually in effect for the project after a quota increase and will continue to
show the default quota of 167 requests per minute.

### Historical Notes

For Echo:
We calculated that the `load_data_batch` should be set to 1245 without a quota increase.

From previously run ingest tasks, we were able to determine how many samples (genomes) are loaded by a single thread in
an hour by examining the logs for one thread. That gave us 10-15 minutes as an estimate.
From there we determined the theoretical CreateWriteStream usage of a single thread per hour, and then calculated how
many threads we could support based upon that.

In Echo there are 414,830 non-control samples that need to be loaded during the bulk ingest process.

We used 3 streams per sample (variant data, reference data and header info).
Based on Google documentation, the default is 10,000 streams per hour, and we needed three streams per sample.

We are going to scatter 414830 samples by 333 wide. So 1245 is our batch size
414830 / 333 = 1245.

For our scale test:
Google increased the quota to 4k per 4 hours.
Roughly 2.5K / hour or 60K samples / day. Scattered 500 wide while using two WriteAPI streams -- one for variant data
and one for reference data
