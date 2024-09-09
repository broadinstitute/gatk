# Storage vs Regeneration Costs

## Prepare tables

These numbers assume the `GvsPrepareRangesCallset` workflow is invoked with the `only_output_vet_tables` input set
to `true`. If this is not the case, meaning the prepare version of the ref ranges table was also generated, all costs
below should be multiplied by about 4:

* Running
  GvsPrepareRangesCallset: [$429.18](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0)
```
-- Look in the cost observability table for the bytes scanned for the appropriate run of `GvsPrepareRanges`.
SELECT
  ROUND(event_bytes * (5 / POW(1024, 4)), 2) AS cost, -- $5 / TiB on demand https://cloud.google.com/bigquery/pricing#on_demand_pricing
  call_start_timestamp
FROM
  `aou-genomics-curation-prod.aou_wgs_fullref_v2.cost_observability`
WHERE
  step = 'GvsPrepareRanges'
ORDER BY call_start_timestamp DESC
```
* Storing prepare data: $878.39 / month
    * Assuming compressed pricing, multiply the number of physical bytes by $0.026 / GiB.

## Avro files

The Avro files generated from the Delta callset onward are very large, several times the size of the final Hail VDS.
For the ~250K sample Delta callset the Avro files consumed nearly 80 TiB of GCS storage while the delivered VDS was
"only" about 26 TiB.

Approximate figures for the ~250K sample Delta callset:

* Avro storage cost: $1568 / month (might be lower if we can get a colder bucket to copy them into)
    * `76.61 TiB  gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/submissions/c86a6e8f-71a1-4c38-9e6a-f5229520641e/GvsExtractAvroFilesForHail/efb3dbe8-13e9-4542-8b69-02237ec77ca5/call-OutputPath/avro`
* [Avro generation cost](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0):
  $3000, 12 hours runtime.

## Hail VariantDataset (VDS)

The Hail VDS generated for the Delta callset consumes about 26 TiB of space in GCS at a cost of approximately $500 /
month. Recreating the VDS from Avro files would take around 10 hours at about $100 / hour in cluster time for a total of
about $1000. Note that re-creating the VDS requires Avro files; if we have not retained the Avro files per the step
above, we would need to regenerate those as well which would add significantly to the cost.

Approximate figures for the ~250K samples Delta callset:

* VDS storage cost: ~$500 / month. Note AoU should have exact copies of the VDSes we have delivered for Delta, though
  it's not certain that these copies will remain accessible to the Variants team in the long term. The delivered VDSes are put here `gs://prod-drc-broad/` and we have noted that we need them to remain there for hot-fixes. The Variants team has
  generated five versions of the Delta VDS so far, one of which (the original) still exist:
    * First version of the callset, includes many samples that were later
      removed `gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/vds/2022-10-19/dead_alleles_removed_vs_667_249047_samples/gvs_export.vds`
* VDS regeneration cost: $1000 (~10 hours @ ~$100 / hour cluster cost) + $3000 to regenerate Avro files if necessary.

