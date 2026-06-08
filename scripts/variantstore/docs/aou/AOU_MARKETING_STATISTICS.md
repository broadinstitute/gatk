# Introduction

For the AoU v9 srWGS callsets (aka Foxtrot), the Variants team was asked to provide some specific statistics that were
intended for use in the marketing materials for the release. In the event that the Variants team is required to generate
these statistics for future releases, this document serves as a reference for how to generate them.


# Statistics

Much of the material here is drawn from [this](https://broadworkbench.atlassian.net/browse/VS-1565) JIRA ticket.

## Total number of genetic variants from the VDS

This is conceptually similar to the following in BigQuery:

```bigquery
    SELECT DISTINCT filters, COUNT(*)
    FROM `<project>.<dataset>.filter_set_sites`
    WHERE filter_set_name = '<filter_set_name>'
    GROUP BY 1
    ORDER BY 1
```

The reason this statistic specifies "from the VDS" is to only count variant sites that are not carried solely by control
or withdrawn samples, which the above query would not filter out.

To calculate this statistic, use the `scripts/variantstore/scripts/filter_set_site_stats_from_vds.py` script. If run
with `--skip-non-ref-filter` this should be efficient enough to run in Terra Hail cluster with two workers even for
AoU-scale VDSes, no massive autoscaling cluster required. `--skip-non-ref-filter` should be safe to use with any VDS
that has been produced with the standard GVS VDS-making scripts that clean up monomorphic reference sites.

## Number of variants within coding regions

This statistic requires an appropriate Exome split MatrixTable; ask in `#dsp-variants` where to find this. Then within
a Jupyter notebook:

```python
import hail as hl
v9_exomeMT = hl.read_matrix_table("<path to split exome MT>")
v9_exomeMT.count_rows()
```
## Total number of variants new to the All of Us dataset (compared to the previous dataset)

This is a fairly simple bit of Hail code but it processes a *lot* of data, so it has been wrapped up in a WDL that will
summon an autoscaling Hail cluster. See `scripts/variantstore/wdl/GvsCountVariantsNewToFoxtrot.wdl` for details, the
only required inputs are the paths to the previous and current VDSes.

## Total number of variants only found in the All of Us dataset over time (compared to gnomAD)

By far the most challenging part of generating this statistic is collecting the VID data for gnomAD into BigQuery.
There is an existing copy of gnomAD 4.1 VIDs in BigQuery at `aou-genomics-curation-prod.foxtrot.gnomad_vids`. Assuming
this is still the current version of gnomAD, run the following query:

```bigquery
SELECT count(distinct vid) FROM `<project>.<dataset>.<vat_table>`
where vid NOT in (
  select vid from `aou-genomics-curation-prod.foxtrot.gnomad_vids`
)
and gvs_all_sc >= 3;
```

If this is *not* the current version of gnomAD, then the following *conceptual* steps are required. Note that the scale
of the data here is quite large (gnomAD 4.1 VCFs are ~800 GiB) and could benefit from a scatter in a WDL:

1. Download the gnomAD sites-only VCFs (all chromosomes) from: `gs://gcp-public-data--gnomad/release/<version>/vcf/joint/`.
1. For each file (chromosome) do something like
    ```
    gunzip $file | grep -v '^#' | cut -f 1,2,4,5 | tr $'\t' '-' | sed 's/chr//g'
    ```
    and put that output into a new file, one per chromosome. This is then a list of the 'vids' (corresponding to the vids
    in our VAT) from gnomAD.
1. Concatenate all the chromosome-specific files together into one large file and load that into BigQuery as a table
    of vids. e.g. `gnomad_vids_<new_version>`.
1. Execute the above query, replacing `aou-genomics-curation-prod.foxtrot.gnomad_vids` with the new table name.
