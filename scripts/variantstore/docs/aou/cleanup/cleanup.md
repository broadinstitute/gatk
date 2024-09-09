# AoU callset cleanup

## Overview

The current Variants policy for AoU callsets is effectively to retain all versions of all artifacts forever. As the  storage costs for these artifacts can be significant (particularly from Delta onward), the Variants team would like to make the cost of retaining artifacts more clear so conscious choices can be made about what to keep and what to delete.

As a general rule, any artifacts that have clearly become obsolete (e.g. VDSes with known issues that have been  superseded by corrected versions, obsolete sets of prepare tables, etc.) should be deleted ASAP. If it's not clear to the Variants team whether an artifact should be cleaned up or not, we should calculate the monthly cost to preserve the artifact (e.g. the sum of all relevant GCS or BigQuery storage costs) as well as the cost to regenerate the artifact.

Reach out to leadership with these numbers for his verdict on whether to keep or delete.

## Specific AoU GVS Artifacts

During the course of creating AoU callsets several large and expensive artifacts are created:

* Pilot workspace / dataset
    * for the Delta callset the Variants team created an AoU 10K workspace and dataset to pilot the Hail/VDS creation
    * these processes will mature to the point where some or all of the contents of this workspace and dataset should be deleted
* Production BigQuery dataset
    * for each previous callset, there was (at least) one new dataset created
    * the dream is to keep the same dataset for multiple callsets and just add new samples, regenerate the filter and create new deliverables, but that has yet to happen because of new features requested for each callset (e.g. update to Dragen version, addition of ploidy data, different requirements to use Hail...etc.)
* Prepare tables
    * needed for VCF or PGEN extract
    * only variant tables are used by `GvsCallsetStatistics.wdl` for callset statistics deliverable
* Sites-only VCFs
    * the VAT is created from a sites-only VCF and its creation is the most resource-intensive part of the VAT pipeline
    * clean up any failed runs once the VAT has been delivered and accepted
* Avro files (Delta onward)
    * huge, several times larger than the corresponding Hail VDS
    * as long as the VDS has been delivered and accepted, they can be deleted
* Hail VariantDataset (VDS)
    * we used to create it in the Terra workspace and then copy it
    * WDL was updated to create in the "deliverables" GCS bucket so there is only one copy of each one 
    * clean up any failed runs once the VDS has been delivered and accepted
* PGEN/VCF Intermediate Files
    * PGEN: multiple versions of the PGEN files are created by GvsExtractCallsetPgenMerged.wdl because it delivers files split by chromosome
    * VCF: only one version of the VCF files and indices are created, but check for failed runs

## Internal sign-off protocol

The Variants team currently has the following VDS internal sign-off protocol:

1. Generate a VDS for the candidate callset into the "delivery" bucket.
1. Open up the VDS in a beefy notebook and confirm the "shape" looks right.
1. Run `GvsPrepareRangesCallset.wdl` to generate a prepare table of VET data
1. Run `GvsCallsetStatistics.wdl` to generate callset statistics for the candidate callset using the prepare VET created in the preceding step
1. Copy the output of `GvsCallsetStatistics.wdl` into the "delivery" bucket.
1. Email the paths to the VDS and callset statistics to Lee/Wail for QA / approval

## Storage versus regeneration costs

### Prepare tables

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

### Avro files

The Avro files generated from the Delta callset onward are very large, several times the size of the final Hail VDS.
For the ~250K sample Delta callset the Avro files consumed nearly 80 TiB of GCS storage while the delivered VDS was
"only" about 26 TiB.

Approximate figures for the ~250K sample Delta callset:

* Avro storage cost: $1568 / month (might be lower if we can get a colder bucket to copy them into)
    * `76.61 TiB  gs://fc-secure-fb908548-fe3c-41d6-adaf-7ac20d541375/submissions/c86a6e8f-71a1-4c38-9e6a-f5229520641e/GvsExtractAvroFilesForHail/efb3dbe8-13e9-4542-8b69-02237ec77ca5/call-OutputPath/avro`
* [Avro generation cost](https://docs.google.com/spreadsheets/d/1fcmEVWvjsx4XFLT9ZUsruUznnlB94xKgDIIyCGu6ryQ/edit#gid=0):
  $3000, 12 hours runtime.

### Hail VariantDataset (VDS)

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

