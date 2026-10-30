# BigQuery Cost Report Script

Small utility to estimate BigQuery pipeline costs from each dataset's `cost_observability` table.

## What it does

For each input `project.dataset`, the script queries:
- `CreateAltAlleles` (Query Scanned)
- `GvsCreateFilterSet` (Query Scanned)
- `GvsCreateFilterSet` (Storage API Scanned)
- `GvsPrepareRanges` (all bytes, billed as Query Scanned)
- `GvsExtractCallset` (Storage API Scanned)

Then it prints per-stage GiB/TiB and estimated USD cost, plus table total.

If any stage query fails, the script aborts and exits non-zero so partial totals are not reported as complete.

## Requirements

- Python 3.9+
- BigQuery access to each `<project.dataset>.cost_observability`
- Application Default Credentials configured (`gcloud auth application-default login`)

Install:

```bash
pip install -r requirements.txt
```

## Usage

Single dataset:

```bash
python bq_cost_report.py aou-genomics-curation-prod.parquet_20k_scale_test_2
```

Multiple datasets:

```bash
python bq_cost_report.py project_a.dataset_x project_b.dataset_y
```

Optional billing/job project:

```bash
python bq_cost_report.py --project my-billing-project project_a.dataset_x
```

Input validation:
- Each dataset must be in `project.dataset` format.
- Invalid identifiers are rejected before query execution.

## Example output

```text
========================================================================
  Table: aou-genomics-curation-prod.parquet_20k_scale_test_2.cost_observability
========================================================================

  Populate Alt Allele
    Step  : CreateAltAlleles
    Type  : Query Scanned  ($6.25/TiB)
    Size  :    14,624.13 GiB  ->    14.2814 TiB
    Cost  : $   89.2590

  Create Filter Set - Query Scanned
    Step  : GvsCreateFilterSet
    Type  : Query Scanned  ($6.25/TiB)
    Size  :    14,624.13 GiB  ->    14.2814 TiB
    Cost  : $   89.2590

  Create Filter Set - Storage API
    Step  : GvsCreateFilterSet
    Type  : Storage API Scanned  ($1.10/TiB)
    Size  :       132.91 GiB  ->     0.1298 TiB
    Cost  : $    0.1428

  Prepare Ranges
    Step  : GvsPrepareRanges
    Type  : Query Scanned  ($6.25/TiB)
    Size  :    33,155.18 GiB  ->    32.3781 TiB
    Cost  : $  202.3633

  Extract - Storage API
    Step  : GvsExtractCallset
    Type  : Storage API Scanned  ($1.10/TiB)
    Size  :    38,327.06 GiB  ->    37.4288 TiB
    Cost  : $   41.1717

  --------------------------------------------------
  Table total:  $422.1958
========================================================================
```

Note: values depend on data in `cost_observability` and current pricing constants in `bq_cost_report.py`.
