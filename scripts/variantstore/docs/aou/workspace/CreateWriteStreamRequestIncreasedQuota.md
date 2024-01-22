## Request Additional Quota for Bulk Ingest
BigQuery Storage API CreateWriteStream quota for small regions per minute per region for region : us-central1
(the region for our project is us-central1 so that is the quota we needed to request an update on)

Note that once the quota has been adjusted, it will not necessarily be visible as anything other than the default on this page.

| Quota                                                                                       | Default                                                                                                | Notes                                                                                                                                                          | 
|---------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `CreateWriteStream` requests  | 10,000 streams every hour, per project per region | You can call CreateWriteStream up to 10,000 times per hour per project per region. Consider using the default stream if you don't need exactly-once semantics. |
For more information, see [Google's quota documentation](https://cloud.google.com/bigquery/quotas)

## File a Google Support Ticket
You can point to [this ticket](https://console.cloud.google.com/support/cases/detail/v2/47548796?project=broad-dsde-methods) as an example. It is summarized below.

Case title: In the project terra-e40eac26 we are hitting theBigQuery Write API quota "RESOURCE_EXHAUSTED"..." (see below in details)
Category: Technical > BigQuery > API - Storage API and Streaming
Project ID: projects/broad-dsde-methods
Observed error message:
com.google.api.gax.rpc.ResourceExhaustedException: io.grpc.StatusRuntimeException: RESOURCE_EXHAUSTED: Exceeds 'CreateWriteStream requests' quota, user_id: project0000002c6ed4d198_us (status: INSUFFICIENT_TOKENS), you can issue a raise quota request through Google Cloud Console. 
Be sure to include this full error message in the request description. Before requesting an increase, please ensure your existing streams have enough utilization in terms of write traffic. For guidelines see https://cloud.google.com/bigquery/docs/write-api-best-practices#limit_the_rate_of_stream_creation. Entity: projects/aou-genomics-curation-prod/datasets/delcho_v2/tables/ref_ranges_001

Case description:
In the project terra-e40eac26 we were using the BigQuery Write API to load data into BigQuery and had many failures with the message "RESOURCE_EXHAUSTED: Exceeds 'CreateWriteStream requests' quota, user_id: project0000002c6ed4d198_us (status: INSUFFICIENT_TOKENS), you can issue a raise quota request through Google Cloud Console." 
When I used the Google Cloud Console to investigate, the Quota "CreateWriteStream requests quota for us region per minute," showed 0% usage for the past seven days. We are contacting you to get this quota increased.
Note that this support ticket is just like 46139869, for which I believe the ultimate outcome was a quota increase: "Our product team have completed the quota change to 30K per 4 hours based on our discussion".
Also note that our usage for this project is very bursty, we anticipate having a high volume of CreateWriteStream requests for several days and then it should die down to a much more limited amount.


## Calculate Quota To be Requested
The default quota for this value is 50. The number that you will want for the quota is likely closer to 500 and will require some calculations to land on precisely.
This is the value that will be input for `load_data_batch` and is necessary to input for large amounts of data.

The calculations initially got us to 500 write streams for two distinct types of stream. 
Currently, however, there are three possible types of write streams that we use. One for variant information (into the vet_xxx tables), a second for reference information (into the ref_ranges_xxx tables) and lastly a sample information stream that will put the headers into place (if that feature is turned on).

Google has told us that it would be technically feasible to apply 45k per 4 hours, but we've settled on 30k per 4 hours.

For Echo:
We calculated that the `load_data_batch` should be set to 1245 without a quota increase.

In Echo there are 414,830 non-control samples that need to be loaded during the bulk ingest process.
We are going to be using 3 streams per sample (variant data, reference data and header info -- sometimes we only load the first 2).
Based on Google docs, the default is 10,000 streams per hour, and we need three streams per sample.


We are going to scatter 414830 samples by 333 wide. So 1245 is our batch size
414830 / 333 = 1245.

Our quota increase request will be 30k instead of 10k. If we get it, we can tripe our `load_data_batch` size.




For scale test:
Google increased the quota to 4k per 4 hours.
Roughly 2.5K / hour or 60K samples / day. Scattered 500 wide while using two WriteAPI streams -- one for variant data and one for reference data
