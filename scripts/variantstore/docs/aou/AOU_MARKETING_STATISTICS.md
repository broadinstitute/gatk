# Introduction

For the AoU v9 srWGS callset (aka Foxtrot), the Variants team was asked to provide some specific statistics for use in
the marketing materials for the release. In the event that the Variants team is required to generate any of these
statistics in the future, this document serves as a reference for how to do that.

# Statistics

Much of the material here is drawn from [this JIRA ticket](https://broadworkbench.atlassian.net/browse/VS-1565).

## Total number of genetic variants - from the VDS

This statistic isn't actually counting "genetic variants", but rather counting variant *sites* grouped by their filters.
This is conceptually similar to the following in BigQuery:

```bigquery
    SELECT DISTINCT filters, COUNT(*)
    FROM `<project>.<dataset>.filter_set_sites`
    WHERE filter_set_name = '<filter_set_name>'
    GROUP BY 1
    ORDER BY 1
```

The reason this statistic specifies "from the VDS" is to only count variant sites that are not carried solely by control
or withdrawn samples, which the simple SQL query above would not filter out.

To calculate this statistic, use the `scripts/variantstore/scripts/filter_set_site_stats_from_vds.py` script. If run
with `--skip-non-ref-filter` this should be efficient enough to run in Terra Hail cluster with two workers even for
AoU-scale VDSes, no massive autoscaling cluster required. `--skip-non-ref-filter` should be safe to use with any VDS
that has been produced with the standard GVS VDS-making scripts that clean up monomorphic reference sites.

## Number of variants within coding regions

This statistic requires an appropriate split Exome MatrixTable; ask in `#dsp-variants` where to find this. Then within
a Jupyter notebook:

```python
import hail as hl
v9_exomeMT = hl.read_matrix_table("<path to split Exome MT>")
v9_exomeMT.count_rows()
```

## Total number of variants new to the All of Us dataset (compared to the previous dataset)

This is another statistic that purports to count variants but is really counting variant *sites* (loci) that are not
found in the previous dataset.
This is a fairly simple bit of Hail code but it processes a *lot* of data, so it has been wrapped up in a WDL that will
summon an autoscaling Hail cluster. See `scripts/variantstore/wdl/GvsCountVdsNovelLoci.wdl` for details. The
only required inputs are the paths to the previous and current VDSes.

## Total number of variants only found in the All of Us dataset over time (compared to gnomAD)

By far the most challenging part of generating this statistic is collecting the VID data for gnomAD into BigQuery.
If the version of gnomAD you want to compare against is still 4.1, then you are in luck: there is an existing copy of
gnomAD 4.1 VIDs in BigQuery at `aou-genomics-curation-prod.foxtrot.gnomad_vids`. In this case, you can simply run
the following:

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
1. For each file (chromosome) do something like the following to convert gnomAD VCF data into the same VID format as
   used in the VAT:
    ```
    zcat $file | grep -v '^#' | cut -f 1,2,4,5 | tr $'\t' '-' | sed 's/chr//g' > ${file%.vcf.gz}.tsv
    ```
    Note that the current version of this logic assumes that multi-allelics are split in the input VCF, which is true
    for gnomAD 4.1.
1. Stage all of these TSV files in a GCS bucket.
1. Create a new BigQuery table to hold the VID-formatted gnomAD data:
   ```bigquery
    CREATE TABLE `<your_project_id>.<your_dataset_name>.gnomad_vids_<gnomad version>` (
         vid STRING
    );
   ```
1. Load all the VID-formatted gnomAD TSVs into the new BigQuery table:
   ```shell
   bq load \
       --source_format=CSV \
       --field_delimiter='\t' \
       <your_project_id>:<your_dataset_name>.gnomad_vids_<gnomad version> \
       "gs://bucket/path/to/vid/formatted/gnomad/data/*.tsv" \
       vid:STRING
   ```
1. Run the query above joining the VAT table to this new gnomAD VID table.

# Notes

## Total number of genetic variants - from the VDS

This is counting loci, not variants. That appears to be the intent given the “number of sites” terminology in the
"Genomics counts" spreadsheet, but the description “total number of genetic variants” is misleading. Note this does not
identify or filter loci that are present only due to alleles with failing allele filters.

## Number of variants within coding regions

This counts variants by counting the rows in the split exome MT. This MT may include all variants within the exome
intervals, including those with failing allele or site filters.

## Total number of variants new to the All of Us dataset (compared to the previous dataset)

This is counting novel loci, not novel variants. Given a lack of context it’s not clear if this was intentional or not.
This code could be modified to count novel variants if that’s actually what’s wanted here, as well as filtering on
failing allele or site filters.

## Total number of variants only found in the All of Us dataset over time (compared to gnomAD)

This is comparing the count of hard-filtered variants seen in at least 3 samples (gvs_all_sc >= 3) found uniquely in the
VAT with respect to gnomAD. This is in contrast with the previous statistics that made no attempt to exclude variants
or loci on the basis of failing allele or site filters, or on the number of samples carrying an allele.
