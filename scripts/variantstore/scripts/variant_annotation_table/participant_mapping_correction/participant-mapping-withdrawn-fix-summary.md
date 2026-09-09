# Withdrawn participants in the Foxtrot participant mapping table

## Symptom

A collaborator reported unexpected person IDs in the participant mapping table that were not represented in the corresponding VDS. Two sampled IDs both turned out to be withdrawn samples, and both were present in the most recent mapping table. The issue was caught during pre-release validation, so no delivered artifact was affected.

## Root cause

All three mapping-table workflows derive `person_ids` by joining `alt_allele` to `sample_info` on `sample_id` with no filter on `withdrawn`: `GvsCreateParticipantMappingTable.wdl`, `GvsMapUnmappedVIDs.wdl`, and `GvsMapDroppedDuplicateVIDs.wdl`.

The leak is `alt_allele`. `populate_alt_allele_table.py` does filter `withdrawn is NULL`, but only for rows it is inserting at that moment, and `GvsPopulateAltAllele.wdl` computes its watermark as `SELECT IFNULL(MAX(sample_id), 0) FROM alt_allele`. The table is strictly append-only for new sample IDs, so nothing ever revisits rows for a sample withdrawn after those rows were written. Echo samples withdrawn during the Foxtrot cycle therefore kept their `alt_allele` rows and flowed into the mapping table.

This is not confined to the newer unmapped-VID and dropped-duplicate steps, which was the initial hypothesis. `GvsCreateParticipantMappingTable` has never had the filter in any revision, and it is the larger source of bad rows. Participant mappings were produced manually before Foxtrot, and it was initially believed that the previous release had a manual corrective step for this which was lost in the WDLization. The data does not support that. The Echo dataset from which the previous VAT was delivered contains zero withdrawn samples — 414,838 rows comprising 414,830 active non-control, 8 controls, and none withdrawn — and its mapping table covers exactly those 414,830 active participants with the 8 controls correctly absent. That release never exercised this path, and a search of the codebase turned up no prior handling either. The accurate statement is therefore that Foxtrot is the first release in which samples were withdrawn, so it is the first to exercise this path; the WDLs were written without the filter and no earlier release had the opportunity to reveal it. This is an unanticipated case in new automation rather than a regression, and no institutional knowledge was lost.

## Fix

Commit `533078da8` adds `AND si.withdrawn IS NULL AND si.is_control = false` to the `ON` clause of all three joins. `ON` rather than `WHERE`, because two of the three already append a conditional `WHERE` for `range_filter`. All three pass `womtool validate`.

`is_control` is included so the predicate matches exactly what `GvsExtractAvroFilesForHail.wdl` puts in the VDS. Controls were not leaking, but only incidentally: `SAFE_CAST(sample_name AS INT64)` yields NULL for a non-numeric control name and `IGNORE NULLS` drops it. That silent guard is precisely what fails for withdrawn AoU samples, whose names are numeric person IDs.

## Remediation of the delivered table

The WDL fix only affects future runs, so the existing table was scrubbed in place. Re-running the pipeline would have been far more expensive.

1. Verified the scrub logic against a synthetic fixture covering all-good, mixed, all-withdrawn, numeric-named control, re-ingested person, duplicate entries, and empty-array cases. Mutation-tested it: four deliberate breaks of the logic were all caught.
2. Confirmed `vid` is unique in the mapping table, since the rebuild groups by `vid`.
3. First attempt used an explode-and-regroup formulation. It consumed 83 slot-days across 82 minutes at the full 2000-slot on-demand cap and never finished. Diagnosis from the query plan: the scan, join, and repartition stages completed in 18 minutes; the output stage then wrote 11,120 of 20,008 units in nine minutes and collapsed from ~1,233 units/min to ~0.4, flat for the following 44 minutes. Cause is shuffle volume, not partition skew — exploding the arrays turns 1.6 billion rows into ~2.58 trillion, each carrying a repeated `vid` string.
4. Switched to a row-local formulation that filters each array in place with no shuffle and no `GROUP BY`. The natural form tests membership against the ~535k good person IDs as a correlated array, which was still too slow. Inverting to test against the 4,875 withdrawn/control IDs, inlined as a constant array literal that BigQuery can hash at plan time, made it near-instantaneous **on the 1% sample**. Do not expect that at full scale: the full-table run of the same query took 18:33. Anyone re-running something of this shape — which is likely, since the other gates below need similar treatment — should not read a query that fails to return promptly as a sign it has gone wrong.
5. Confirmed the inverted form is equivalent: it differs from the good-set form only for person IDs absent from `sample_info` entirely, and a check on a 1% sample found zero such IDs.
6. Profiled the table from a 1% `TABLESAMPLE`: ~1.6 billion VIDs, ~2.58 trillion person-ID entries, mean array 1647, median 2, max 540,539.
7. Built and validated the whole approach on the 1% sample before committing to the full run.
8. Full run completed in 18:33. Swapped via `ALTER TABLE ... RENAME TO`, retaining the pre-scrub table as a backup.

## Validation

The check that gates correctness is whether `sample_info.withdrawn` actually implies absence from the delivered VDS. Nothing in the code enforces this — `GvsWithdrawSamples` stamps a timestamp in BigQuery while `GvsMergeAndRescoreVDSes` removes samples off a separately maintained TSV. Exported the VDS column list and compared it against `sample_info`.

- VDS sample count 535,662, exactly equal to the active non-control count, and every VDS column joined to a `sample_info` row.
- Withdrawn samples still present in the VDS: 0. Active samples missing from the VDS: 0. Controls in the VDS: 0.
- `sample_info` ledger reconciles exactly: 535,662 active non-control, plus 4,875 numeric all-withdrawn names, plus 8 non-numeric control names, equals 540,545 rows.
- Entries removed: predicted 0.902% from the withdrawn/control share of samples, observed 0.903% from the byte delta. Agreement in both directions rules out over-removal and under-removal.
- Post-scrub sample check: zero withdrawn or control person IDs across 25.65 billion sampled entries covering 16,017,178 VIDs. An unscrubbed person would be expected to appear ~48,000 times in a sample that size.
- Coverage: zero active non-control samples lost all mapping entries. Conclusive for the whole table, since the scrub is row-local and the sampled scrubbed table is a true subset of the full one.
- Row count 1,601,242,198 before, 1,601,242,196 after.

## Secondary finding: two VIDs removed entirely

`1-143186828-T-TG` and `17-72146919-GAGAAG...-G` each had exactly one carrier, withdrawn and confirmed absent from the VDS, so both rows emptied and were dropped.

Tracing the first one found a second defect, caused by two GVS artifacts using different normalization conventions. Both minimize allele representations; only the VAT left-aligns. `alt_allele` holds minimized representations, so an indel is trimmed but stays wherever the caller placed it. VAT VIDs are minimized AND left-aligned, because `GvsCreateVATFromVDS.wdl:731` runs `bcftools norm -f <ref>`, which does both. The two therefore coincide only when the caller happened to emit the leftmost equivalent position already.

That is what happened here. The withdrawn sample's caller emitted the deletion at the leftmost position, so its minimized `alt_allele` row coincides with the VAT VID exactly, while the active carrier's row sits elsewhere. For `1-143186828-T-TG` the withdrawn sample is at `chr1:143186828 T>TG` and the active carrier at `chr1:143186829 G>GG`. For `17-72146919-GAGAAG...-G` the active carrier's row is a 31bp-to-1bp deletion at `chr17:72146925`, six bases right, which is two periods of the surrounding `AAG`/`GAG` repeat. The base table joined on the VID's own coordinates, matched the withdrawn sample, and never saw the active carrier. The delivered mapping named the wrong person, not merely an extra one.

`GvsMapUnmappedVIDs` could not repair this because it selects VIDs with `vid NOT IN (SELECT vid FROM mapping)` — a partial match is still a match, so the VID never looked unmapped.

This is narrower than it first appears and is not an independent defect. Write R_LA for the left-aligned form, which is always the VAT VID. The base table join matches every `alt_allele` row at R_LA regardless of which sample contributed it, which gives four cases.

| Case | Condition                                                                   | Base join finds                 | Outcome                                                                                                                 |
|------|-----------------------------------------------------------------------------|---------------------------------|-------------------------------------------------------------------------------------------------------------------------|
| 1    | Active carriers all use R_LA                                                | all active carriers             | Correct.                                                                                                                |
| 2    | Active carriers use R_LA and other forms                                    | the R_LA subset only            | The VDS holds a record per form, so normalization collides them onto one VID and GvsMapDroppedDuplicateVIDs repairs it. |
| 3    | No active carrier uses R_LA, and no non-VDS sample sits at R_LA             | nothing, so the VID is unmapped | GvsMapUnmappedVIDs repairs it.                                                                                          |
| 4    | No active carrier uses R_LA, but a withdrawn or control sample sits at R_LA | only that non-VDS sample        | Broken. A row is written, so the VID never looks unmapped and nothing repairs it.                                       |

The stale `alt_allele` row is the cause rather than an incidental detail — it manufactures a spurious match that suppresses the detection which would otherwise have repaired the VID. Commit `533078da8` therefore fixes this instance: with the filter, case 4 yields no rows, no mapping row is written, and `GvsMapUnmappedVIDs` selects the VID and maps the real carriers. The `is_control` half matters as well, since a control-only match produces a row with an empty `person_ids` array that still counts as mapped.

Case 4 has a SECOND trigger that the commit does NOT fix. The spurious R_LA row need not come from a withdrawn or control sample; it is enough that the R_LA representation was excluded from the VAT while R_x was kept. VETS scores representations independently, because `filter_set_info` is keyed by location/ref/alt and nothing in the filtering path knows two representations are synonyms, so R_LA can fail calibration sensitivity while R_x passes. The VAT is hard filtered and `GvsCreateVATFromVDS.wdl` normalizes `filtered_sites_only.bcf`, so filtering happens before normalization: only R_x survives to be normalized, one record results, no duplicate is detected, and `GvsMapDroppedDuplicateVIDs` does not fire. `alt_allele` still holds R_LA regardless of filtering, and its carrier is an ordinary active sample, so the withdrawn/control predicate passes them through and the VID is mapped to a carrier whose variant failed filtering while the passing carrier is missed.

The general statement is that the base join must be restricted to exactly the population the VAT represents, on two axes: sample eligibility (active and non-control, which the commit adds) and variant quality (passing the callset's calibration-sensitivity threshold, which nothing does today). A mismatch on either axis produces the same spurious match and the same suppression of unmapped detection. Unlike the withdrawn trigger, this one is detectable by joining mapping rows back to `filter_set_info` and looking for source alleles that fail; a bounded sample would establish whether it occurs in Foxtrot in practice.

### Restoration of the two rows

Both VIDs were restored by hand rather than by re-running the pipeline. For each, the synonymous representation was resolved by left-aligning the candidate VDS alleles with `bcftools norm -f <ref>` and taking the one that lands on the VID's position, then looking up that representation's active carriers in `alt_allele`. Each returned exactly one active carrier, matching the carrier count the VDS reports for the allele, and in both cases a different person from the withdrawn one the row previously named.

```sql
INSERT INTO `foxtrot.vid_to_participant_mapping_2026_08_28` (vid, person_ids)
VALUES ('17-72146919-GAGAAGAAGAAGAAGAAGAAGAAGAAGAAGA-G', [<PERSON_A>]);

INSERT INTO `foxtrot.vid_to_participant_mapping_2026_08_28` (vid, person_ids)
VALUES ('1-143186828-T-TG', [<PERSON_B>]);
```

The row count is therefore back to 1,601,242,198, matching the pre-scrub table. Note that these two rows are correct only for the representations found; a pipeline re-run with commit `533078da8` would derive them automatically, since with the withdrawn filter in place both VIDs come out unmapped and `GvsMapUnmappedVIDs` maps them.

## The general defect: the base join applies none of the VDS's inclusion criteria

Everything above is one root cause. `GvsCreateParticipantMappingTable` derives `person_ids` by joining `alt_allele` to `sample_info`, and `alt_allele` holds every call ever ingested. Three separate gates determine what actually reaches the VDS and therefore the VAT, and the join applied none of them.

1. **Sample eligibility.** Withdrawn and control samples are excluded from the VDS by `GvsExtractAvroFilesForHail.wdl`. This is the FIRST defect the collaborator reported, and the one fixed by `533078da8`.
2. **Genotype confidence.** `import_gvs.py:262` nulls the genotype when `GQ` is 0: `LGT = hl.parse_call(hl.or_missing(hl.is_missing(var_ht.GQ) | (var_ht.GQ != 0), var_ht.GT))`. A GQ0 call exists in `vet` and `alt_allele`, and `LA` is even populated in the VDS, but the sample is not a carrier there. Not applied by the join. This is the SECOND defect the collaborator reported, raised separately as a participant count mismatch at a specific VID and traced to this cause; see the worked example below.
3. **Variant and genotype quality.** `FT` is set per genotype at `import_gvs.py:380` as `~any_no & (any_yes | all_ok)`, from `yng_status` and calibration sensitivity. Filtered genotypes stop contributing to `AC`, and `GvsCreateVATFromVDS.wdl:726` drops `AC=0` records. Not applied by the join.

A collaborator's spot check of `13-32368001-C-CTT` measured all three at one VID, and the chain reconciles exactly:

- 1,067 participants in the mapping table
- 1,063 carriers in the VDS. The gap of 4 is exactly the GQ0 calls at that allele, confirmed as `call_GT = '0/1'` with `call_GQ = 0`, and there are exactly 4 such calls.
- 937 of those pass `FT`. The gap of 126 is genotypes failing `FT`, all of them het-var `1/2` calls filtered on the strength of their *other* allele, since a genotype is judged by its worst allele.

The collaborator reported the 4. They did not report the 126, because their count did not consider `FT`.

**Which of those 130 belong in the mapping table is not a judgment call, because the VAT settles it.** `gvs_all_ac` for this VID is 1,289, which is the allele count computed from the 937 FT-passing genotypes and not from the 1,063 VDS carriers (1,415) or the 1,067 `alt_allele` carriers (1,419). The mapping table ships alongside the VAT, so the population it should list is the one the VAT counts. At this VID that is 937, and the mapping table over-reports by 130 participants, or 12%.

That gives the remaining work a concrete target rather than an argument about intent: `person_ids` should contain exactly the participants whose genotypes contribute to the VAT's `AC`.

Gate 2 is a one-line fix mirroring the import logic: `AND (aa.call_GQ IS NULL OR aa.call_GQ != 0)`. Gate 3 is larger — it needs a join to `filter_set_info`, the callset's calibration-sensitivity thresholds, and per-sample aggregation over each sample's alleles at the locus, since `FT` is `any_yes | all_ok` across the whole genotype rather than a property of one allele.

Two incidental facts worth recording, both of which cost time here. The VAT is keyed by `(vid, transcript)`, so a VID returns one row per transcript and any per-VID count needs a `DISTINCT`. And `gvs_all_an` at this VID is 453,164, implying roughly 226,582 samples with a confident passing genotype out of 535,662 in the callset — the gates bite hard at a poly-T tract, which is likely why a discrepancy surfaced on an early spot check.

### Negative result: normalization synonyms are rare

The two VIDs this scrub emptied were misattributed through a different mechanism, described above: `alt_allele` and the VAT normalize differently, so a vid can match `alt_allele` rows belonging to a representation the VAT excluded. That mechanism is real and confirmed for those two. We looked for further instances and found none, but the search was **not exhaustive** and should not be read as evidence of absence. The methodology proved cumbersome and slow, and the effort was set aside when the collaborator reported the second defect above.

Sampling the mapping table for VIDs whose own allele cannot pass calibration sensitivity gave ~508 candidates per 1% sample, extrapolating to roughly 54,000. That figure did not survive scrutiny. 190 of the 508 were same-length SNPs and MNPs, which have no alignment ambiguity and cannot exhibit the mechanism at all. Of the indels examined individually, every one resolved to het-var rescue — a failing allele riding into the VDS on a `1/2` genotype paired with a `Y` allele at the same site — or to unrelated alternates at the same position. None had a passing synonym at a shifted position.

So the prevalence of the normalization mechanism remains unquantified. What can be said is that the GQ0 and `FT` gates are demonstrated to account for real discrepancies with measured counts, while normalization has two confirmed instances and no established rate. If the question matters, the search is worth resuming with a better method than the one we abandoned.


## Follow-ups

- PR commit `533078da8`.
- Verify whether `dropped_duplicate_mappings.tsv` contains every synonym in a cluster or only the dropped ones. `GvsMapDroppedDuplicateVIDs` deletes the base row and re-inserts carriers joined on `dup.input_*` only, so if it carries only the dropped synonyms it discards carriers found at the surviving representation. Case two above rests entirely on that workflow and the commit does not touch it.
- Consider adding the post-scrub assertion as a pipeline task: no `person_id` in the mapping table belongs to a withdrawn or control sample. About a dollar per run against a 1% sample. This matters more than it would if a safeguard had merely been lost: there was never prior knowledge to carry forward, so nothing but an explicit check will catch the next case that has not yet arisen.
- `GvsCreateParticipantMappingTable` still uses the explode-and-regroup shape that failed at this scale. The commit adds the filter without changing that structure, so the workflow remains on the failing path for a table this size. Folding in the row-local formulation is worth doing before the next callset.
- Drop `_prescrub_backup` once the collaborator confirms.
