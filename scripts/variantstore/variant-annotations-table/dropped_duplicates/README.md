# Introduction

The purpose of ticket VS-1683 was to investigate what looked like some questionable handling of duplicate data in
the VAT-making process, [these lines](https://github.com/broadinstitute/gatk/blob/882585aadc026f9263a43d91a970275b8954116b/scripts/variantstore/variant-annotations-table/GvsCreateVATfromVDS.wdl#L611-L630)
overall and specifically [these lines](https://github.com/broadinstitute/gatk/blob/882585aadc026f9263a43d91a970275b8954116b/scripts/variantstore/variant-annotations-table/GvsCreateVATfromVDS.wdl#L628-L630).
From the previous `pseudo_vids_only_in_vat` work, we know that there were recurring issues with variants in input gVCFs
that were not left-aligned. That non-left-aligned data was also what necessitated special duplicate handling, but given
some of the characteristics of the data that resulted from this investigation, we may want to change the way we're
handling these duplicates.

Because the underlying issues here and with "pseudo VIDs" are so similar, the work in this directory is able to reuse
some of the previous work in the sibling directory `pseudo_vids_only_in_vat`.

## If this code drops duplicates, how did we end up with "pseudo VIDs" with multiple synonyms?

Actually our "pseudo VIDs" do not represent duplicates as far as this deduplication code is concerned, even when
multiple synonyms of a canonical left-aligned representation are found. This is because the deduplication code
[first filters out alleles with `AC=0`](https://github.com/broadinstitute/gatk/blob/882585aadc026f9263a43d91a970275b8954116b/scripts/variantstore/variant-annotations-table/GvsCreateVATfromVDS.wdl#L618)
before [running normalization](https://github.com/broadinstitute/gatk/blob/882585aadc026f9263a43d91a970275b8954116b/scripts/variantstore/variant-annotations-table/GvsCreateVATfromVDS.wdl#L623)
that can introduce duplicates. For all "pseudo VIDs", there is exactly one synonym with `AC>0`; all other synonyms have
`AC=0`.

## How do the observed allele counts of the "pseudo VIDs" compare to the allele counts of these dropped duplicates?

These dropped duplicates appear to have allele counts roughly 1000 times higher than the "pseudo VIDs". Looking at Echo
pseudo VIDs with the highest allele counts:

```sql
SELECT
  SUM(CAST (REGEXP_EXTRACT(info_field, "AC=([0-9]+);") AS int64)) AS affected_allele_count,
  left_aligned_location,
  left_aligned_ref,
  left_aligned_alt,
FROM
  `aou-genomics-curation-prod.echo.pseudo_vid_mapping`
GROUP BY
  left_aligned_location,
  left_aligned_ref,
  left_aligned_alt
ORDER BY
  affected_allele_count DESC
LIMIT
  5
```

```
affected_allele_count	left_aligned_location	left_aligned_ref	left_aligned_alt
787	16000089732863	GCCCTGACCCTGCTGTACCCCGCTGTACCCTGCGCCCTCGCCCTCTGCTGTGTTCACCCCGACCCTGCTGTACCCTGCGCCCTCATCCTCTGCTGTGCCCTCCCCCTCTGCTGTGTTCGCCCTGACCCTGCTGTACCCTGCGCCCTCGCCCTCTGCTGTGTTTGCCCTGACCCTGCTGTACCCTGCGCCCTCGCCCTCTGCTGTGTTCA	G
204	24000011676582	CAATGGACTGGAATGGACAGGAAAGGAATGGACTGGATTGGAAAGGATTTGAAGGGAATGTGCTAGAGTGGAACGGACTCGAATGGAATGGAAACGAATGGAATGGAATGAAATGGAATGGAATGGAATGGAAAGGAATAGCCTGGAATGAAATGGGATGGAAAGCAATTGAATGGAATGGAATGGAGTCGAATGGAATAAAATCGAATGGAATGGCATCGAATGGAATGGAGTGGAATGGAATGGACTCGAATGG	C
203	19000035375413	G	GA
201	1000001028899	T	TG
201	4000049142988	ATTCCATTCCGTTCCTTTTATCTTGGGTTGATTCCATTCCATTCAATTCCTTTCCGTTCCGTTCCGTTCCATTCCATTCCATTCCATTCCATTCCATTGCATTCCACTGGGGTTGATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATACCCTTCGGGTTGATTCCTTTCCAATCCATTCCATTCTATACCATTCCACTTCATTCCGTTCCATTCCTTTCGGGTTGATTCTGTTCCATTCCATGCCCTTTTATTCCATTCCATTCCATTCCATACCATTCCACCAAAGTTGATTGCATGTTATTCCATTCCATTCCATTCCATTCCATTCCT	A
```

Pointing this query to the `vs_1683_vat_dropped_duplicates` table instead returns:

```
affected_allele_count	left_aligned_location	left_aligned_ref	left_aligned_alt
737565	11000134839958	G	GC
700817	4000076917565	A	AG
333129	3000192641140	C	CT
293037	13000112864345	G	GC
243444	11000018467958	TATATATATACGTATATATATACGTATATATATGC	T
```

Wrapping the above query with another `SUM`:

```sql
SELECT
  SUM(affected_allele_count)
FROM (
  SELECT
    SUM(CAST (REGEXP_EXTRACT(info_field, "AC=([0-9]+);") AS int64)) AS affected_allele_count,
    left_aligned_location,
    left_aligned_ref,
    left_aligned_alt,
  FROM
    `aou-genomics-curation-prod.echo.pseudo_vid_mapping`
  GROUP BY
    left_aligned_location,
    left_aligned_ref,
    left_aligned_alt
  ORDER BY
    affected_allele_count DESC )
```

returns totals of 13081 for pseudo VID AC and 3862660 for dropped duplicate AC.

## Procedure to find this data

I started by copying interesting things like `track_dropped.tsv` and `stderr` from the successful Echo VAT creation:

```shell
gsutil -m rsync -r -x '^(?!.*(/track_dropped.tsv$|/stderr$)).*' gs://<Echo GCS Bucket>/submissions/<Submission ID>/GvsCreateVATfromVDS/<Workflow ID>/call-RemoveDuplicatesFromSitesOnlyVCF/ .
```

Fancy `find` and `sort` to locate and order non-empty `track_dropped.tsv` files:

```shell
find shard-* -name track_dropped.tsv -size +0 | xargs cat | sort -k 1.4n -k 2n  > input.data
```

Use the Gemini-written script `fish_for_synonyms.py` to find all synonyms (including canonical left alignments!) of
the variants that were dropped:

```shell
 cat input.data | python fish_for_synonyms.py > commands.sh
 sh commands.sh > out.vcf
```

Normalize:

```shell
bcftools norm out.vcf -f Homo_sapiens_assembly38.fasta > out_norm.vcf
```

Sneakily make a version of the normalized VCF that only includes that sites that were dropped:

```shell
cat <(bcftools head sites-only.vcf.gz) <(grep -wFf input.data out_norm.vcf) > out_norm_dropped.vcf
```

A bit of explanation for the above: the first `<(...)` is just going to write a VCF header. The second `<(...)` is
using the `input.data` file as a set of strings to search in the `out_norm.vcf`. This works because `input.data` is
actually [the first 5 columns of the VCF in order](https://github.com/broadinstitute/gatk/blob/882585aadc026f9263a43d91a970275b8954116b/scripts/variantstore/variant-annotations-table/GvsCreateVATfromVDS.wdl#L630),
and contains only those entries that were present in the `track_dropped.tsv`s.

Then run the `compare_vcfs.py` script (horrible name) from the pseudo VID investigation to make a BigQuery TSV load file:

```shell
python compare_vcfs.py out.vcf out_norm_dropped.vcf > vat_dropped_duplicates.tsv
```

And load following the pseudo VID instructions (different target table, of course):

```shell
bq load --project_id aou-genomics-curation-prod --source_format=CSV --skip_leading_rows=1 --field_delimiter="\t" \
                    echo.vs_1683_vat_dropped_duplicates vat_dropped_duplicates.tsv \
                            vid:STRING,chr:STRING,input_location:INTEGER,input_position:INTEGER,input_ref:STRING,input_alt:STRING,left_aligned_location:INTEGER,left_aligned_position:INTEGER,left_aligned_ref:STRING,left_aligned_alt:STRING,info_field:STRING
```