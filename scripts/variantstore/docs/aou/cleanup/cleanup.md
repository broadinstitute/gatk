# AoU callset cleanup

## Overview

The current Variants policy for AoU callsets is effectively to retain all versions of all artifacts forever. As the
storage costs for these artifacts can be significant (particularly from Delta onward), the Variants team would like to
make the cost of retaining artifacts more clear so conscious choices can be made about what to keep and what to delete.

As a general rule, any artifacts that have clearly become obsolete (e.g. VDSes with known issues that have been
superseded by corrected versions, obsolete sets of prepare tables, etc.) should be deleted ASAP. If it's not clear to
the Variants team whether an artifact should be cleaned up or not, we should calculate the monthly cost to preserve the
artifact (e.g. the sum of all relevant GCS or BigQuery storage costs) as well as the cost to regenerate the artifact.
Reach out to Lee with these numbers for his verdict on whether to keep or delete.

## Specific AoU GVS Artifacts

During the course of creating AoU callsets several large and expensive artifacts are created:

* Pilot workspace / dataset
    * For the AoU Delta callset the Variants team created an AoU 10K workspace and dataset to pilot the Hail-related
      processes we were using for the first time. At some point these processes will mature to the point where some or
      all of the contents of this workspace and dataset should be deleted. This is likely not an issue for discussion
      with Lee as it is internal to the Variants team's development process, but we should be mindful to clean up the
      parts of this that we are done using promptly.
* Production BigQuery dataset
    * The plan from Delta forward is to use the same production BigQuery dataset for all future callsets. This was also
      the plan historically as well, but in practice that didn't work out for various reasons. For Delta in particular
      the Variants team was forced to create a new dataset due to the use of drop state NONE for Hail compatibility. If
      we are forced to create another BigQuery dataset in the future, discuss with Lee to determine what to do with the
      previous dataset(s). In the AoU Delta timeframe, the Variants team additionally reached out to a contact at VUMC
      to determine if various datasets from the Alpha / Beta eras were still in use or could possibly be deleted.
* Prepare tables
    * These tables are used for callset statistics only starting with Delta but were previously used for extract as well
      in pre-Delta callsets. Per 2023-01-25 meeting with Lee it appears this VET prepare table and associated callset
      statistics tables can be deleted once the callset statistics have been generated and handed off.
* Terra workspace
    * It seems that VAT workflows may generate large amounts of data under the submissions "folder". e.g. ~10 TiB of
      data under this folder in the AoU 10K workspace (!). At the time of this writing the VAT process is not fully
      defined so this item may especially benefit from updates.
* Avro files (Delta onward)
    * These are huge, several times larger than the corresponding Hail VDS. It's not clear that there's any point to
      keeping these files around unless there was a bug in the Hail GVS import code that would require a patch and
      re-import. Per 2023-01-25 meeting with Lee we have now deleted the Delta versions of these Avro files. Per the
      preceding comments, going forward these files can be deleted once the Variants team feels reasonably confident
      that they won't be needed for the current callset any longer.
* Hail VariantDataset (VDS) (Delta onward)
    * The Variants team creates a copy of the VDS and then delivers a copy to the AoU preprod datasets bucket. That copy
      of the VDS seems to stay in the delivery location for at least a few days, but it's not clear if that copy gets
      cleaned up after AoU later copies the VDS to a production bucket. The Variants team should not rely on this copy
      of the VDS being available long-term. Per 2023-01-25 meeting with Lee, we have retained the oldest (with AI/AN +
      controls) and most recent versions (without AI/AN or controls, corrected phasing and GT) of the Delta VDS. This
      can serve as our team's guidance for how to handle multiple VDS versions going forward, though of course we can
      always ask Lee for explicit guidance.

## Internal sign-off protocol

The Variants team currently has the following VDS internal sign-off protocol:

* Generate a VDS for the candidate callset
* Run validation on this VDS
* Run `GvsPrepareRangesCallset` to generate a prepare table of VET data
* Generate callset statistics for the candidate callset using the prepare VET created in the preceding step
* Forward VDS and callset statistics to Lee for QA / approval

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

