# Tracing Back Orphaned VIDs in the Echo VAT

## Introduction

Investigation done in a Terra notebook terminal.
The procedure described in this section can be applied to any of the unmapped VIDs, although the example is worked with
one particular VID.

## General procedure

This first example will utilize VID `2-15219938-C-CTATA`, a small insertion. Search for variants with the same "shape"
(insert/deletion size) nearby, downstream in the sites-only VCF. In this particular case we search for an insertion
of four bases within 10 bases of the VID position:

```shell
# Parse a VID into components.
parse_vid() {
    local vid="$1"
    # The variables below are assigned and visible to the caller.
    IFS='-' read -r chr pos ref alt <<< "$vid"
    echo "chr=$chr pos=$pos ref=$ref alt=$alt"
}

vat_vid="2-15219938-C-CTATA"
parse_vid $vat_vid

search_range=10
insert_len=$((${#alt}-${#ref}))
start_pos=$((pos + 1))
end_pos=$((pos + search_range - 1))

# The `ILEN` filter is used to select out variants of a specific size, insertions being positive and deletions being negative.
bcftools query --include "(ILEN=${insert_len})" --format "${chr}\-%POS\-%REF\-%ALT" \
  --regions chr${chr}:${start_pos}-${end_pos} sites-only.vcf
```

This returns:

```
2-15219939-T-TATAT
```

Let's confirm that the left-aligned form of this variant corresponds to our VID of interest. Select out the variant
and left align:

```shell
bcftools view --include "(ILEN=${ilen})" --regions chr${chr}:${start_pos}-${end_pos} sites-only.vcf |
    bcftools norm -f Homo_sapiens_assembly38.fasta  | bgzip > left_aligned.vcf.gz
```

Index then query:
```shell
bcftools index left_aligned.vcf.gz
bcftools query --regions chr${chr}:${pos} --format "${chr}\-%POS\-%REF\-%ALT" left_aligned.vcf.gz
```

This returns:
```
2-15219938-C-CTATA
```

Which is the VID we were looking for. Now use the handy `sample_data_table` in BigQuery (scraped from the Echo callset's
`sample` Terra data table) to find paths to input reblocked gVCFs and the unreblocked gVCFs from which they were made:

```shell
# From our findings above
gvs_vid="2-15219939-T-TATAT"
parse_vid $gvs_vid

bq --apilog=false query --project_id=aou-genomics-curation-prod --format=prettyjson --use_legacy_sql=false "

  SELECT dt.reblocked_gvcf, dt.gvcf_path
  FROM echo.alt_allele aa
  JOIN echo.sample_info si ON aa.sample_id = si.sample_id
  JOIN echo.sample_data_table dt ON dt.research_id = si.sample_name
  WHERE
    location = ${chr} * 1000000000000 + ${pos}
    AND ref = '${ref}'
    AND allele='${alt}'

" | jq '.[0]' | tee gvcfs.json
```

This prints a "pretty" JSON object while also saving it to a file called `gvcfs.json`. The output will look something like this:

```json
{
  "reblocked_gvcf": "gs://path/to/my/reblocked_gvcf.gvcf.vcf.gz",
  "gvcf_path": "gs://different/path/to/my/unreblocked_gvcf.gvcf.vcf.gz"
}
```

To copy these files and their indexes to the :

```shell

for gvcf in $(jq -r '.[] | values' gvcfs.json)
do
    gcloud storage cp "${gvcf}*" .
done
```

Assign these variables for convenience:

```shell
reblocked_gvcf=$(basename $(jq -r '.reblocked_gvcf' gvcfs.json))
unreblocked_gvcf=$(basename $(jq -r '.gvcf_path' gvcfs.json))
```
