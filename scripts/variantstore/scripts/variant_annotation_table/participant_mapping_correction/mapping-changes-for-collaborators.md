# Changes to the VID to participant mapping table in v9_r2_p4

This note describes what has changed in these v9_r2_p4 mappings since the previous v9_r2_p3 version, as well as what to expect when comparing them.

## The short version

For each VID the mappings now list only the participants whose genotypes contribute to that VID's `gvs_all_ac` in the VAT. Previously the mappings listed every participant with any record of that allele in the callset's source tables, including records the callset itself does not count toward `gvs_all_ac`.

**The v9_r2_p4 mappings only remove participants from the previous v9_r2_p3 mappings.** No VID gained a participant. Any count derived from the participant mappings should only go down or stay the same in the move from v9_r2_p3 to v9_r2_p4.

## What was removed

Two classes of genotype were incorrectly included in the previous mappings even though the callset does not count them:

**GQ 0 no-calls.** In the v9 srWGS deliverables a genotype with `GQ` 0 is a *no-call*. The participant has no genotype at that site and contributes nothing to `AC`. However the underlying variant record still exists in the GVS source tables and the v9_r2_p3 mappings erroneously included those participants as carriers. This is the source of the discrepancy seen when comparing the v9_r2_p3 mappings against the VDS: the mappings named participants the VDS does not call.

**Genotypes failing `FT`.** A genotype that fails the filter does not contribute to `AC`. This did not surface in the v9_r2_p3 participant count comparisons because neither the mappings nor a straightforward VDS carrier count considered `FT`; both sides of the comparison made the same error of ignoring `FT`. Consideration of `FT` was always required for a tieout against `gvs_all_ac`, and will now be required for correct tie outs of participant counts to the v9_r2_p4 mappings.

## `FT` is a property of the genotype, not of the allele

In the srWGS v9 deliverables a genotype is judged by its **worst** non-reference allele. So a `1/2` het-var genotype fails filtering when *either* of its two alleles fails; in prior srWGS versions a `1/2` het-var would only fail filtering if *both* alleles failed filtering. So in srWGS v9 deliverables a participant will stop being mapped to a VID if its partner allele in a `1/2` het-var fails filtering.

Please note: If an analysis ignores `FT` it will disagree with the v9_r2_p4 mapping table and overcount associated participants.

## How much to expect

Measured against the v9_r2_p3 mappings:

| Quantity                               | Value             |
|----------------------------------------|-------------------|
| Entries in the v9_r2_p3 mappings       | 2,559,382,818,198 |
| Entries removed for v9_r2_p4           | 2,341,606,369     |
| Fraction of entries removed            | 0.0915%           |
| VIDs whose entry changed               | 5,437,323         |
| Fraction of VIDs changed               | 0.34%             |
| Mean entries removed, per affected VID | 431.2             |
| Most entries removed from a single VID | 535,659           |

Among the VIDs that changed, the median lost 2.4% of its participants, and the 99th percentile lost 86.8%. So the typical affected VID barely moves, and the effect is concentrated in a small number of VIDs rather than spread evenly.

The row count goes from 1,601,242,198 to 1,601,242,026. The 172 missing rows are VIDs that had no remaining carriers once the corrections described above were applied, so they no longer appear in the mapping table.

