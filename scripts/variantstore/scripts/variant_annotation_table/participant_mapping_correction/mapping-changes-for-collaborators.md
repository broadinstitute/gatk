# Changes to the Foxtrot VID-to-participant mapping table

This note describes what changed between the mapping table you received previously and the corrected one, and what to expect when you compare them. It is written to be read on its own; nothing else in this directory is needed.

## The short version

The mapping table now lists, for each VID, exactly the participants whose genotypes contribute to that VID's `gvs_all_ac` in the VAT. Previously it listed every participant with any record of that allele in the callset's source tables, including records the callset itself does not count.

**Entries were only ever removed, never added.** No VID gained a participant. Any count you derive from the table can only have gone down.

## What was removed

Two classes of genotype were being mapped even though the callset does not count them:

**GQ 0 no-calls.** A genotype with `GQ` 0 is a no-call in the VDS — the participant has no genotype at that site and contributes nothing to `AC`. The underlying record still exists in the source table and was still being mapped.

**Genotypes failing `FT`.** A genotype that fails the filter model stops contributing to `AC`. These were also still being mapped.

The second one is the one worth internalizing, because it does not behave the way people expect.

## `FT` is a property of the genotype, not of the allele

A genotype is judged by its **worst** non-reference allele. So a `1/2` het-var genotype is filtered when *either* of its two alleles fails — including when the allele you are looking at passes cleanly and its partner allele, at the same site, is the one that failed.

The consequence: a participant can disappear from a VID whose own allele has no quality problem at all.

This is easy to miss because het-vars are invisible in a split representation. Splitting a multi-allelic site rewrites a `1/2` as `0/1` at each allele's row, so a split view shows those participants as ordinary heterozygotes. At `13-32368001-C-CTT`, a split count gives 711 carriers where the unsplit truth is 283 `0/1` plus 428 `1/2`; all 126 participants removed there by `FT` are het-vars. In a split view they look like ordinary het carriers vanishing for no reason.

If your analysis counts carriers from a split representation and ignores `FT`, it will disagree with the corrected mapping table, and the disagreement will be concentrated at multi-allelic sites.

## How much to expect

Measured against the previously delivered table:

| Quantity                               | Value             |
|----------------------------------------|-------------------|
| Entries in the delivered table         | 2,559,382,818,198 |
| Entries removed                        | 2,341,606,369     |
| Fraction of entries removed            | 0.0915%           |
| VIDs whose entry changed               | 5,437,323         |
| Fraction of VIDs changed               | 0.34%             |
| Mean entries removed, per affected VID | 431.2             |
| Most entries removed from a single VID | 535,659           |

Among the VIDs that changed, the ratio of old count to new count has a median of 1.02, a 99th percentile of 7.6, and a maximum of about 535,000. So the typical affected VID barely moves, and the effect is concentrated in a small number of VIDs rather than spread evenly.

The row count goes from 1,601,242,198 to 1,601,242,026. The 172 missing rows are VIDs that had no remaining carriers once the above was applied, so they no longer appear in the table.

## Checking it yourself

`vds_carriers_report.py` in this directory takes a list of VIDs, makes one pass over the VDS, and reports per VID both the count that ignores `FT` and the count that applies it, along with the reasons for the difference. Pass `--mapping-tsv` with the per-VID array lengths from the file you received and it will check the delivered table against the VDS directly:

```shell
python3 vds_carriers_report.py \
    --vds-path gs://.../foxtrot.vds \
    --vids-file my_vids.tsv \
    --mapping-tsv mapping_lengths.tsv \
    --output report.tsv
```

The `agrees_mapping` column is the answer to "is the file I received right for this VID". `--vat-ac-tsv` adds a check against the VAT's `gvs_all_ac`, which is the stronger one — it is built by a separate pipeline, so agreement means the `FT` and GQ 0 handling reproduces what the VAT actually did rather than merely being self-consistent.

Run it on a Hail cluster; it reads the VDS.
