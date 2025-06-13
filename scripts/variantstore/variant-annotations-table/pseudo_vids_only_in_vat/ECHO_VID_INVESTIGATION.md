# Tracing Back Orphaned VIDs in the Echo VAT

## Introduction

Early in this investigation it became clear that many orphaned "pseudo vids" that did not correlate with entries in the
GVS `alt_allele` or `filter_set_info` tables actually existed elsewhere in these tables in synonymous but non-left
aligned forms. The goal of this analysis was to confirm if all the "pseudo vids" had non-left aligned representations
and to establish the correspondences between them.

In fact correspondences were found in GVS for all "pseudo vids", with the data describing their relation to existing GVS
entries uploaded to the table `aou-genomics-curation-prod.echo.pseudo_vid_mapping`. In all cases it appears that the
non-left aligned variants came directly from the reblocked gVCFs used as input to GVS, which in turn took this data from
the original unreblocked gVCFs. The `README.md` file that is a sibling to this document describes the procedure used to
populate this `pseudo_vid_mapping` table.

Further analysis allowed for the creation of the `aou-genomics-curation-prod.echo.pseudo_vid_sample_id` table which
links orphaned "pseudo vids" to referring `sample_id`s in the `echo.alt_allele` table. This in turn enabled the
`SearchGVCFsForUnmappedVIDs.wdl` workflow to be created, which captured relevant lines in associated gVCFs to create
the `aou-genomics-curation-prod.echo.pseudo_vid_gvcf_content` table.

All the investigation steps described below were performed in a Terra notebook terminal within the AoU security perimeter.

## General procedure - `2-15219938-C-CTATA`

This first detailed example uses VID `2-15219938-C-CTATA`, a small insertion. Search for variants with the same "shape"
(insert/deletion size) nearby, downstream in the sites-only VCF. In this particular case we search for an insertion of
four bases within 20 bases of the VID position:

```shell
# Parse a VID into components.
parse_vid() {
    local vid="$1"
    # The variables below are assigned and visible to the caller.
    IFS='-' read -r chr pos ref alt <<< "$vid"

    insert_len=$((${#alt}-${#ref}))
    start_pos=$((pos + 1))
    end_pos=$((pos + search_range))

    # The "pseudo vid"s are either short (< 10 base) inserts, or longer (up to ~400 base) deletions.
    # Choose the search range accordingly.
    if [[ ${insert_len} -gt 0 ]]
    then
        search_range=20
    else
        search_range=200
    fi

    echo "chr=$chr pos=$pos ref=$ref alt=$alt insert_len=$insert_len start_pos=$start_pos end_pos=$end_pos search_range=$search_range"
}

search_sites_only() {
    # The `ILEN` filter is used to select out variants of a specific size, insertions being positive and deletions being negative.
    bcftools query --include "(ILEN=${insert_len})" --format "${chr}\-%POS\-%REF\-%ALT" \
        --regions chr${chr}:${start_pos}-${end_pos} sites-only.vcf
}

vat_vid="2-15219938-C-CTATA"
parse_vid $vat_vid
search_sites_only
```

This returns:
```
2-15219939-T-TATAT
```

Note that this particular example returns a single value, but several of these "pseudo vids" in Echo will return
multiple results. Each "pseudo vid" would need to be processed separately per the process described below.

Let's confirm that the left-aligned form of this variant corresponds to our VID of interest. Select out the variant
and left align:

```shell
search_sites_only_left_aligned()  {
    # Make a left aligned version from this query
    bcftools view --include "(ILEN=${insert_len})" --regions chr${chr}:${start_pos}-${end_pos} sites-only.vcf |
        bcftools norm -f Homo_sapiens_assembly38.fasta  2>/dev/null | bgzip > left_aligned.vcf.gz

    # Index then query specifically at the VID position
    bcftools index left_aligned.vcf.gz
    bcftools query --regions chr${chr}:${pos} --format "${chr}\-%POS\-%REF\-%ALT" left_aligned.vcf.gz
}

search_sites_only_left_aligned
```

This returns:
```
2-15219938-C-CTATA
```

Which is the VID we were looking for. Now download the input reblocked gVCF and the unreblocked gVCF from which it was made
and see how the data appears there:

```shell
# From our findings above
gvs_vid="2-15219939-T-TATAT"
parse_vid $gvs_vid

download_gvcfs() {
    if [[ "${chr}" == "X" ]]
    then
        alt_allele_chr=23
    elif [[ "${chr}" == "Y" ]]
    then
        alt_allele_chr=24
    else
        alt_allele_chr=${chr}
    fi

    # Use the handy `sample_data_table` in BigQuery (scraped from the Echo callset's `sample` Terra data table) to find
    # paths to input reblocked gVCFs and the unreblocked gVCFs from which they were made:
    bq --apilog=false query --project_id=aou-genomics-curation-prod --format=prettyjson --use_legacy_sql=false "

        SELECT dt.reblocked_gvcf, dt.gvcf_path
        FROM echo.alt_allele aa
        JOIN echo.sample_info si ON aa.sample_id = si.sample_id
        JOIN echo.sample_data_table dt ON dt.research_id = si.sample_name
        WHERE
            location = ${alt_allele_chr} * 1000000000000 + ${pos}
            AND ref = '${ref}'
            AND allele='${alt}'

    " | jq '.[0]' | tee gvcfs.json

    # A GCS-enabled htslib/bcftools can query GCS-hosted gVCFs directly (see SearchGVCFsForUnmappedVIDs.wdl
    # as an example), but the bcftools currently in the default Terra environment does not support GCS.
    for gvcf in $(jq -r '.[] | values' gvcfs.json)
    do
        gcloud storage cp "${gvcf}*" .
    done
    # Make sure to give the indexes a more recent modification time than the gVCFs so bcftools doesn't get upset.
    touch *.tbi

    reblocked_gvcf=$(basename $(jq -r '.reblocked_gvcf' gvcfs.json))
    unreblocked_gvcf=$(basename $(jq -r '.gvcf_path' gvcfs.json))
}

search_reblocked_gvcf() {
    bcftools query --regions chr${chr}:${pos} --format "%CHROM\t%POS\t%REF\t%ALT" ${reblocked_gvcf}
}

search_unreblocked_gvcf() {
    bcftools query --regions chr${chr}:${pos} --format "%CHROM\t%POS\t%REF\t%ALT" ${unreblocked_gvcf}
}

```

Now look in these files using queries similar to the ones we ran before against the sites-only VCF. First look in the
reblocked gVCF that is the actual input to GVS:

```shell
download_gvcfs
search_reblocked_gvcf
```

This returns:
```
chr2    15219939        T       TATAT,<NON_REF>
```

So we can see that the reblocked gVCF has the same position, ref, and alt as in the `alt_allele` table in GVS. Next the
unreblocked gVCF:

```shell
search_unreblocked_gvcf
```

This returns:

```
chr2    15219939        TGGCCGGGCAGAGGGCTCCTCACTTCCCAGTAGGGGCGGCCGGGCAGAGGCGCCCCTCACCTCCCGGACGGGGCGGCTGGCCAGGCGGGGGGCTGATCCCCCCACCTCCCTCCCGGACGGGGCGGCTGGCCGGGCGGGGGGCTGACCCCCCCCACCTCCCTCCTGGACGGGGCGGCTGGCCGGGCGGGGGGCTGACCCCCCCACCTCCCTCCCGGACGGGGCGGCTGGCCGGGC GGGGGGCTGACCCCCCCACCTCCCTCCCGGACGGGGCGGCTGGCCGGGCAGAGGGGCTCCTCACTTCCCAGTAGGGGCGGCCGGGCAGAGGCGCCCCTCACCTCCCGGACGGGGCGGCTATAT     T,TATATGGCCGGGCAGAGGGCTCCTCACTTCCCAGTAGGGGCGGCCGGGCAGAGGCGCCCCTCACCTCCCGGACGGGGCGGCTGGCCAGGCGGGGGGCTGATCCCCCCACCTCCCTCCCGGACGGGGCGGCTGGCCGGGCGGGGGGCTGACCCCCCCCACCTCCCTCCTGGACGGGGCGGCTGGCCGGGCGGGGGGCTGACCCCCCCACCTCCCTCCCGGACGGGGCGGCTGGCCGGGCGGGGGGCTGACCCCCCCACCTCCCTCCCGGACGGGGCGGCTGGCCGGGCAGAGGGGCTCCTCACTTCCCAGTAGGGGCGGCCGGGCAGAGGCGCCCCTCACCTCCCGGACGGGGCGGCTATAT,<NON_REF>
```

So apparently this was a hetvar site in the original unreblocked gVCF, but in the reblocking process the alleles seem to have been assigned different positions.

## Large deletion 21-10769701-TCCTGAAA...-T

This example happens to be the largest deletion among the orphaned VIDs.

```shell
long_deletion="21-10769701-TCCTGAAA...-T"
parse_vid $long_deletion

search_sites_only
```

Produces output:
```
21-10769704-TGAAA...-T
```

So this finds a deletion of the same size as what we're looking for 3 bases downstream of the left-aligned position.

Applying the same left alignment steps and querying as above:

Searching left aligned:

```shell
search_sites_only_left_aligned
```

produces:

```
21-10769701-TCCTGAAA...-T
```

recovering the original VID we were looking for.

In this case, both files show exactly the same data recorded in GVS, i.e. unlike the previous example this site was not
a hetvar in the unreblocked gVCF.
