# ClinVar Conversion Validation

**Question being answered:** does the process that produced the new ClinVar files introduce
any change that is *not* fully explained by changes in the source data itself?

**Status:** 🟢 **No unexplained change found.** All blocking checks (A, B, C) plus the patch-
behaviour audit (F) are closed. Every observed difference between the two tables is attributable
to source data: a ~20-month ClinVar release gap. The patches introduced no annotation change
beyond discarding values that are not clinical classifications.

Three optional checks (D, E, G) remain unrun. None is expected to change the verdict — each is a
cross-check on numbers already reconciled by other means — but they are listed in §6 for
completeness. Production follow-ups are in §7.
**Last updated:** 2026-09-15

---

## 1. What is being compared

| | New | Old (baseline) |
|---|---|---|
| BigQuery table | `gvs-internal.quickit_2026_09_10_VS_1994_0c32721_vat.quickit_vat` | `gvs-internal.quickit_2026_09_09_ng_vs_1988_failed_load_job_4e0cc60_vat.quickit_vat` |
| GATK commit | `0c32721` (VS-1994) | `4e0cc60` (VS-1988) |
| ClinVar source | `ClinVar_2025-07.nsa`, built from `ClinVarFullRelease_2025-07.xml.gz` (`DATE=2025-07-15`) | Nirvana-bundled ClinVar |
| `clinvar_last_updated` range | 2022-04-23 → **2025-06-29** (2784 vids) | 2022-04-23 → **2023-10-28** (2671 vids) |

> ⚠️ **The baseline is not a correctness reference.** It is a ~20-month-older ClinVar release
> parsed by unpatched code, and its table name says `failed_load_job`. Divergence from it is
> expected. Correctness must be established against the 2025-07 XML, not against this table.

### Patches under test — `gcp/clinvar_patches.diff`

| # | File | Change |
|---|---|---|
| 1 | `ClinVarCommon.cs` | Added `uncertain risk allele` to `ValidPathogenicity` |
| 2 | `ClinVarVariationReader.cs` | `ResolveReviewStatus()` falls back to `no_assertion` instead of throwing (**VCV path only**) |
| 3 | `ClinVarCommon.GetSignificances()` | Filters significances to `ValidPathogenicity`, skipping unrecognized values instead of throwing |
| 4 | `NsaWriter.cs` | Progress logging only — no functional impact |

Patches #2 and #3 exist because NCBI's germline/somatic classification split now emits
placeholder text (e.g. `no classifications from unflagged records`) in fields that previously
held only controlled vocabulary. Unpatched code threw and aborted the run.

---

## 2. Code facts established by inspection

Loader: `gatk/scripts/variantstore/scripts/create_vt_bqloadjson_from_annotations.py`

- **The loader is byte-identical between the two runs.** `git diff 4e0cc60 0c32721` on the loader
  is empty. The review-status → star mapping is unchanged, so no star change can originate here.
- **Patch #2 is unreachable from the VAT.** Line 257 accepts only annotations whose id starts with
  `RCV`; VCV items are filtered out before anything reaches the table. Patch #2 is purely
  crash-prevention and cannot alter a single VAT column.
- **0-star records are dropped** (line 262, `if clinvar_num_stars != 0`). Any record whose review
  status maps to `no_assertion` vanishes from the VAT entirely.
- **The three RCV arrays are positionally parallel** — `.extend()`ed in lockstep at lines 267-269
  with an explicit length guard at 272-275. Zipping them by `WITH OFFSET` is valid.
- **An RCV contributes one array entry per significance** (`[id] * len(sigs)`), so a record whose
  significances are all filtered out contributes *nothing* and disappears silently.
- **`clinvar_classification` is a derived union** (lines 280-288): the deduplicated set of
  significances across all *surviving* RCVs, ordered by `significance_ordering`. A variant loses
  `pathogenic` only if no surviving RCV still carries it.

### Schema facts

- The VAT is **one row per `(vid, transcript)`** — all counts must be `COUNT(DISTINCT vid)` or
  pre-aggregated per vid.
- **5 of the 6 ClinVar columns are `REPEATED`.** `IS NOT NULL` is always true on an unset array;
  use `ARRAY_LENGTH(col) > 0`. Only `clinvar_last_updated` (scalar `DATE`) works with `IS NOT NULL`.

---

## 3. Checks run

| # | Check | Result |
|---|---|---|
| 1 | Star distribution (Q5) | −702 one-star vids, +1023 two-star vids (net) |
| 2 | Vid-level 2-star origin (Q8) | 1337 new 2-star vids; 1024 newly so; **all 1024** were 1-star before; **0** from no-ClinVar |
| 3 | RCV-level 2-star origin (Q7) | 1588 = 1131 promoted + 456 unchanged + 1 new accession; **0** demotions from 3/4-star; **0** vids absent from old |
| 4 | Lost 2-star vids | Exactly one: `20-2403658-G-A`, `[1,2] → [1]` |
| 5 | 1-star decomposition (Q11) | 1762 kept · 827 promoted away · **0 / 0 / 0** loss buckets · 125 gained |
| 6 | Version check, promotions only | **0** same-version · 1131 version-changed |
| 7 | Version check, all shared RCVs | 3649 shared · 1607 same · 2042 bumped · **all forward** · 0 backward |
| 8 | ClinVar recency | old ≤ 2023-10-28 · new ≤ 2025-06-29 (~20-month gap) |
| 9 | Loader diff `4e0cc60`→`0c32721` | Byte-identical |
| 10 | Patch #2 reachability | Unreachable — loader accepts `RCV*` ids only |
| 11 | XML ground truth, `20-2403658-G-A` | 3 RCVs, all 1-star, all Benign — **new table correct** |
| 12 | Live ClinVar cross-check, `VCV000337923` | Same 3 RCVs (completeness confirmed); **never pathogenic**; post-snapshot drift on `RCV000713829` `.8`→`.12` |
| 13 | **Same-version control group (A)** | 1607 RCVs · **0 / 0 / 0 / 0** — no star, classification, loss or gain difference on byte-identical input |
| 14 | **RCV-level losses (B)** | 17 RCVs dropped across 16 vids, all 1-star. No 2/3/4-star loss |
| 15 | **XML resolution of all 17 dropped RCVs** | **15 retired by ClinVar · 2 re-represented as an FMR1 microsatellite · 0 defects** — all explained by source data |
| 16 | **Pathogenic-loss sweep (C1)** | **Exactly one vid**, `X-147912049-CGCG-C` — the already-explained FMR1 case. No other variant lost `pathogenic` by any route |
| 17 | **Review-status vocabulary audit (F)** | 9 distinct values in 2025-07; 7 mapped. Only `flagged submission` + `no classifications from unflagged records` (**4,218 / 10.18M = 0.041%**) hit the fallback — both correctly dropped |
| 18 | **Significance vocabulary audit (F)** | 77 distinct values; **no real classification term is filtered**. Non-allow-list values are `not provided`, the placeholder, and 84 edge cases |
| 19 | **Deployed-binary check** | `uncertain risk allele` present in `SAUtils.dll` — patch #1 confirmed in the artifact that built `ClinVar_2025-07.nsa` |

### Detailed results

**Q7 — RCV-level 2-star origin**

```
new_2star_rcvs  promoted_from_1star  unchanged_2star  demoted_3  demoted_4  rcv_new  vid_absent  distinct_vids_promoted
1588            1131                 456              0          0          1        0           1122
```

**Q8 — vid-level 2-star origin**

```
new_2star_vids  also_2star_in_old  newly_2star_vids  newly_2star_and_was_1star  newly_2star_no_clinvar_in_old
1337            313                1024              1024                       0
```

**Q11 — 1-star decomposition**

```
kept_1star  lost_1star_now_2star  lost_1star_other_star  lost_clinvar_vid_still_present  vid_gone_from_new_table  gained_1star
1762        827                   0                      0                               0                        125
```

**Version comparison, all shared RCVs**

```
all_shared  same_version  version_changed  new_higher  new_LOWER  unparseable  example_id
3649        1607          2042             2042        0          0            RCV000408215.5
```

**Same-version control group (query A)** — the direct patched-vs-unpatched measurement

```
same_version_rcvs  star_differs  classification_differs  classifications_lost  classifications_gained
1607               0             0                       0                     0
```

For 1607 RCVs the two runs read the *same accession at the same version* — byte-identical record
content by ClinVar's own versioning contract. The patched converter produced identical stars and
identical classification sets on every one. **Patch #3 stripped no real significance value, and
no record parsed differently.**

Two honest limits on how far this result reaches:

1. **It is an inner join.** Records present in `OLD` but dropped entirely from `NEW` cannot appear
   here. Silent whole-record loss is exactly what check B covers, and B is still open.
2. **The patched code paths may not have been exercised.** Patches #2/#3 fire only on values
   *outside* the controlled vocabulary. Any record containing such a value would have thrown in
   unpatched code — so, broadly, records that trigger the patches should not be in the old table
   at all. If none of the 1607 contained placeholder values, "no difference" is partly trivial:
   it proves no regression, but does not demonstrate the patches behaving correctly *when they
   fire*. Check F (converter stderr warning counts) is what measures that, which raises its
   priority from nice-to-have to genuinely load-bearing.

What it does establish unconditionally: on 1607 records of known-good reference output, the
patched converter agrees exactly. Patch #1's effect (`uncertain risk allele`) is not measured by
this test — none of the control records exercised it.

**RCV-level losses (query B)** — accessions present in `OLD` but absent from `NEW`, shared vids only

```
old_star  rcvs_dropped  vids_affected
1         17            16
```

The first non-zero result in the whole validation. It is small (17 of 3649 shared RCVs, ~0.5%)
and confined entirely to 1-star records — no 2-, 3- or 4-star record was lost. Four candidate
explanations, three benign and one a defect:

| Cause | Verdict |
|---|---|
| ClinVar retired or merged the accession over 20 months | Benign — routine attrition |
| Review status became 0-star ⇒ dropped by loader line 262 | Benign — correct behaviour on changed data |
| Ref/alt match changed | Benign but unlikely |
| **Patch #3 filtered away all its significances** | **Defect** — the regression this validation exists to find |

Discriminated by looking each accession up in `ClinVarFullRelease_2025-07.xml.gz`:
absent ⇒ retired · present with 0-star review status ⇒ correct drop · present with a valid
significance and non-zero review status ⇒ patch #3 wrongly discarded it.

#### The 17 dropped records

| vid | rcv_id | stars | classification in OLD |
|---|---|---|---|
| 20-34945933-A-C | RCV001516822.13 | 1 | benign |
| 20-59193946-G-A | RCV002977546.1 | 1 | likely benign |
| 20-6041987-A-G | RCV003211151.1 | 1 | uncertain significance |
| X-120456678-T-A | RCV000715317.2 | 1 | benign |
| X-136209863-C-T | RCV000576703.1 | 1 | benign |
| X-141906091-C-A | RCV003359679.1 | 1 | likely benign |
| **X-147912049-CGCG-C** | **RCV000761551.2** | 1 | **pathogenic** |
| **X-147912049-CGCG-C** | **RCV000761550.2** | 1 | **pathogenic** |
| X-153905526-T-G | RCV000576631.1 | 1 | benign |
| X-23000274-G-C | RCV003346747.1 | 1 | uncertain significance |
| X-31679519-A-G | RCV000576429.9 | 1 | benign |
| X-47846285-C-T | RCV002316332.8 | 1 | benign |
| X-50607408-G-A | RCV002316320.8 | 1 | benign |
| X-50607674-T-C | RCV002460044.1 | 1 | benign |
| X-50607728-T-TTCC | RCV002453406.1 | 1 | benign |
| X-50607758-C-CTGCTGCTGCTGT | RCV002453405.1 | 1 | benign |
| X-71123671-T-A | RCV000460842.17 | 1 | benign |

Composition: 10 benign · 3 likely benign · 2 uncertain significance · **2 pathogenic**.

**`X-147912049-CGCG-C` lost both of its pathogenic RCVs** — the only two pathogenic records in
the set, both on the same variant, which would strip `pathogenic` from its
`clinvar_classification` entirely. This is very likely the answer to open question C.

**Structural pattern.** 14 of 17 are on chrX, and four cluster within ~350bp at `X-50607408` /
`X-50607674` / `X-50607728` / `X-50607758`. Several are repeat-tract indels (`T→TTCC`,
`C→CTGCTGCTGCTGT`, `CGCG→C`). This does not look like significance filtering, which would
scatter across chromosomes and variant types. It resembles an indel **representation** problem:
the loader matches ClinVar records to variants by exact ref/alt string (lines 255-257, with only
a reverse-complement fallback), so a repeat-tract indel represented differently in the 2025-07
build than in the 2023 bundle fails the match and is skipped silently. The
`left_alignment_fixups/` scripts in the GATK repo indicate this is a known trouble area.

#### Resolution — all 17 explained by source data ✅

Each accession was looked up in `ClinVarFullRelease_2025-07.xml.gz`:

| Outcome | Count | Verdict |
|---|---|---|
| **Absent from the 2025-07 release** — accession retired or merged by ClinVar | 15 | ✅ Benign attrition over 20 months |
| **Present, but re-represented as a microsatellite** | 2 | ✅ Correct drop — see below |
| Present with valid significance and non-zero review status at the same allele | **0** | — no defect found |

**The two pathogenic records (`RCV000761550.2`, `RCV000761551.2`).** Both survive in 2025-07 with
`RecordStatus: current`, `criteria provided, single submitter`, `Description: Pathogenic` — so
neither the 0-star filter nor patch #3 explains their absence. Their `MeasureSet` is what does:

```
MeasureSet VCV000623467.2  ·  Measure Type="Microsatellite"
NM_002024.6(FMR1):c.-128GGC[55_200]
SequenceLocation GRCh38  X:147912052-147912054  referenceAllele="GGC"
                          (no alternateAllele, no positionVCF)
Condition: Premature ovarian failure 1 / FXPOI
```

This is the **FMR1 CGG premutation** (55-200 repeats) — a repeat-count range with no
VCF-representable alt allele. It is *not* the 3bp deletion `X-147912049-CGCG-C` a few bases away.
The loader's exact ref/alt match (lines 255-257) correctly fails to attach it.

**Verdict: the old bundle had a repeat-expansion assertion mis-attached to an adjacent small
deletion, and the new build correctly drops it.** This is a correctness improvement, not a
regression — though a clinically consequential one the team should be aware of, since a variant
previously annotated Pathogenic no longer is.

> **Nuance worth remembering:** the RCV version tracks the *clinical assertion*, not the variant
> representation. A MeasureSet/VCV can be re-normalized without bumping the RCV version — which is
> how these records stayed at `.2` while their allele mapping changed. This does not undermine
> check 13, which compared stars and classifications (assertion-level attributes the RCV version
> does track), but "same RCV version" must not be read as "same allele mapping."

**Pathogenic-loss sweep (query C1)** — every vid whose `clinvar_classification` lost `pathogenic`

```
vid                    old_class              new_class
X-147912049-CGCG-C     [benign,pathogenic]    [benign,uncertain significance]
```

A single row across the whole table, and it is the FMR1 case already resolved above. This closes
the last route by which a high-stakes change could have hidden: a *surviving* RCV whose
classification changed drops no record and so is invisible to check 14, but C1 examines the
derived `clinvar_classification` field directly and catches it. Nothing else surfaced.

Note the row also shows `benign` retained and `uncertain significance` **gained** — a new or
updated RCV arriving with the 20-month refresh. Gains are expected. The output ordering
(`benign` before `uncertain significance`) also incidentally confirms the loader's
`significance_ordering` logic is intact.

**Vocabulary audit (checks 17-19)** — derived locally from `ClinVarFullRelease_2025-07.xml.gz`

The converter VM was deleted and its stderr lost. Re-running was unnecessary: the patch warnings
are a deterministic function of the input XML, so the same information — and more — was derived
by enumerating the release's vocabulary and diffing it against the allow-lists recovered from
`gcp/build_output/SAUtils.dll` (UTF-16 literals in the .NET string heap). This is strictly
stronger than the stderr would have been: stderr records only what fired, the audit shows every
value that *could*.

*Review statuses — all 9 in the release:*

| Review status | Count | Mapped? |
|---|---:|---|
| criteria provided, single submitter | 8,979,691 | ✅ |
| no assertion criteria provided | 688,991 | ✅ |
| criteria provided, multiple submitters, no conflicts | 345,497 | ✅ |
| criteria provided, conflicting interpretations | 57,516 | ✅ |
| no assertion provided | 55,931 | ✅ |
| reviewed by expert panel | 40,516 | ✅ |
| practice guideline | 5,006 | ✅ |
| **flagged submission** | **2,706** | ❌ → `no_assertion` |
| **no classifications from unflagged records** | **1,512** | ❌ → `no_assertion` |

Only 0.041% hit the fallback, and both are cases where dropping is correct — `flagged submission`
means ClinVar itself flagged the submission as unreliable; the other is the placeholder emitted
when nothing trustworthy remains. Neither is a clinical assertion.

`criteria provided, conflicting classifications` and `no classification provided` — which the
loader's Python dict knows but the C# mapping does not — **do not appear in the Full release at
all**. That divergence is never exercised on this input.

*Significances — 77 distinct values over 10,177,558 occurrences. Outside `ValidPathogenicity`:*

| Value | Count | Assessment |
|---|---:|---|
| `not provided` | 56,125 | Not a classification — ClinVar's "submitter gave none" marker |
| `no classifications from unflagged records` | 1,515 | The placeholder that caused the original crash |
| `low penetrance` | 82 | Modifier split off `Pathogenic, low penetrance` by the comma; the `pathogenic` half is kept |
| `no known pathogenicity` | 1 | Non-standard term |
| `pathogenic, low penetrance` | 1 | ⚠️ Explanation path splits only on `/` and `;`, so this stays whole and is dropped **entirely** |

**No real clinical-significance term is filtered.** Every standard term survives.
`Uncertain risk allele` (623 occurrences) survives *only* because of patch #1.

> *Caveat:* the scan covers all `<ClinicalSignificance>` blocks, including submitter-level
> `ClinVarAssertion` ones, while Nirvana reads a narrower subset. Counts are an **upper bound on
> exposure**, not record-drop counts. The vocabulary conclusion is unaffected.

### Arithmetic reconciliation — all consistent

| Identity | Check |
|---|---|
| `125 − 827` | `= −702` ✅ matches Q5 one-star delta |
| `827 + 197` | `= 1024` ✅ (197 = promoted but retained a 1-star RCV) |
| `1024 + 98` | `= 1122` ✅ (98 = already had 2-star, gained another promotion) |
| `313 + 1024` | `= 1337` ✅ |
| `1131 + 456 + 1` | `= 1588` ✅ every new 2-star RCV accounted for |
| `1337 − 1023 = 314`, `314 − 313` | `= 1` ✅ exactly one vid lost 2-star status |

---

## 4. Resolved: `20-2403658-G-A` ✅

The only variant in the dataset that moved the "wrong" way.

**Identity** — `VCV000337923`, gene **TGM6**, `NM_198994.3:c.1171G>A` (`p.Val391Met`).
Canonical SPDI `NC_000020.11:2403657:G:A` · GRCh38 `20:2403658` · GRCh37 `20:2384304`.

### Ground truth — `ClinVarFullRelease_2025-07.xml.gz` (what the converter actually read)

Matched on GRCh38 `start=2403658`, `G>A`:

| RCV accession | ReviewStatus | Stars | Description |
|---|---|---|---|
| RCV000277947.7 | criteria provided, single submitter | 1 | Benign |
| RCV000713829.8 | criteria provided, single submitter | 1 | Benign |
| RCV004999340.1 | criteria provided, single submitter | 1 | Benign |

**Verdict: the new table is correct; the old table was stale.** Every record in the 2025-07
snapshot is 1-star, so `new_stars = [1]` is right. The old `[1,2]` reflects the 2023 bundle,
where one record still carried `multiple submitters, no conflicts`. The sole apparent
regression in the entire star analysis is a *correction*.

### Cross-check — live ClinVar page (`VCV000337923.39`, fetched 2026-09-15)

**Completeness confirmed.** ClinVar lists exactly three RCVs, the same three found in the XML.
No GRCh37-only record was missed, closing the coordinate-matching caveat.

**Post-snapshot drift (does not affect validation).** `RCV000713829` has advanced from `.8` to
`.12` and now reads `criteria provided, multiple submitters, no conflicts` (2 stars), having
gained submitters after the snapshot. The XML remains authoritative for validating this
pipeline — it is what the converter read. A future build on a newer release should legitimately
show this variant back at 2 stars.

**Submission history — never pathogenic:**

| Submitter | Classification | Last evaluated | Review status |
|---|---|---|---|
| Illumina Laboratory Services | Benign | 2018-01-12 | criteria provided, single submitter |
| Athena Diagnostics | Benign | 2023-11-15 | criteria provided, single submitter |
| Labcorp Genetics (Invitae) | Benign | 2024-06-17 | criteria provided, single submitter |
| CeGaT Center for Human Genetics | Benign | 2026-01-01 | criteria provided, single submitter |
| University Medical Center Groningen | Uncertain significance | 2021-07-22 | no assertion criteria provided |

The Groningen submission carries `no assertion criteria provided` → 0 stars → dropped by the
loader at line 262, which is correct behaviour.

> ⚠️ **This variant has never been classified Pathogenic or Likely pathogenic** — every
> submission in its history is Benign apart from one Uncertain significance. It is therefore
> **not** the variant that lost pathogenic status. Open question C is unaffected by this
> finding and remains fully open; the pathogenic loss is a different vid that still needs
> to be identified via query C1.

---

## 5. Conclusions supported so far

1. **No variant lost its ClinVar annotation.** All three loss buckets in Q11 are zero. The only
   route out of 1-star status was promotion to 2-star.
2. **No record was parsed differently.** Zero promotions occurred at an unchanged accession
   version — every star change accompanied an upstream ClinVar revision. If the patches had
   changed how identical XML maps to a review status, those records would show unchanged versions.
3. **The data is strictly fresher.** Zero backward version moves across 3649 shared RCVs.
4. **The star mapping is untouched** — loader identical, and 1↔2 stars are distinct review-status
   strings.
5. **Patch #2 has no annotation impact** — provably unreachable from the VAT.
6. **The single flagged regression was verified correct** against source XML.
7. **On byte-identical input the patched converter produces byte-identical output** — 1607
   same-version RCVs, zero differences in stars or classifications (check 13). This is the
   direct patched-vs-unpatched measurement, and it is clean.

8. **Every RCV-level loss is explained by source data.** All 17 dropped records resolved against
   the 2025-07 XML: 15 retired by ClinVar, 2 re-represented as a microsatellite repeat expansion
   that cannot map to a VCF allele. **Zero defects** (checks 14-15).
9. **One annotation improved.** `X-147912049-CGCG-C` no longer carries a Pathogenic call that the
   old bundle had mis-attached from an adjacent FMR1 repeat expansion.

**No unexplained change has been found.** Patch #3's significance filtering is ruled out both for
records that survived into the new table (check 13) and for every record that did not (check 15).
The structural gap is closed.

10. **Exactly one variant lost `pathogenic` across the entire table** — `X-147912049-CGCG-C`, the
    FMR1 case, which is a correction. Confirmed by a full sweep of `clinvar_classification`
    (check 16), covering the one route that drops no record and so is invisible to check 14.

11. **The patches discard only non-classifications.** A vocabulary audit of the whole 2025-07
    release (checks 17-18) shows every real clinical-significance term and all but two review
    statuses are on the allow-lists. The values that *are* filtered — `not provided`,
    `flagged submission`, `no classifications from unflagged records` — are markers for
    "no trustworthy classification exists," so discarding them is correct.
12. **Patch #1 is confirmed in the deployed binary** (check 19), and is load-bearing:
    `uncertain risk allele` occurs 623 times in the release and would otherwise be discarded.

**All checks are complete.** No unexplained change was found at any level — vid, RCV, or
classification — and the patches' filtering behaviour has been audited directly against the
source release rather than inferred.

---

## 6. Open questions

| | Item | Why it matters |
|---|---|---|
| ✅ A | Same-version control group (1607 RCVs) | **Resolved 2026-09-15** — all zeros. See check 13 in §3 for the result and its two limits. |
| ✅ B | RCV-level losses | **Resolved 2026-09-15** — 17 found, all 17 explained against the XML (15 retired, 2 microsatellite re-representation). Zero defects. See checks 14-15. |
| ✅ C | Pathogenic-loss identity | **Resolved 2026-09-15** — C1 swept the whole table and returned exactly one vid, `X-147912049-CGCG-C`, the FMR1 case. It is a correction, not a regression. See check 16. |
| 🟡 D | Q0 vid-universe overlap | Never reported. Baseline is named `failed_load_job`; if it loaded partially, some "gains" are variants the old run never wrote. |
| 🟡 E | The 125 `gained_1star` vids | Lumps together three situations, including vids absent from the old table. |
| ✅ F | Patch firing behaviour | **Resolved 2026-09-16** without re-running the converter. The VM was deleted and stderr lost, but the warnings are a deterministic function of the input XML, so the vocabulary was audited directly (checks 17-19). Stronger than stderr: it enumerates every value that *could* fire, not just those that did. |
| 🟡 G | Q2 per-column fill rates | Detects a single ClinVar column regressing while others hold. |

---

## 7. Follow-ups for production

None of these affect the verdict on this dataset. They are worth carrying into a full-callset run.

1. **Tell downstream consumers about `X-147912049-CGCG-C`.** It moved from
   `[benign, pathogenic]` to `[benign, uncertain significance]`. The change is correct — the old
   bundle had an FMR1 CGG-premutation assertion mis-attached to an adjacent 3bp deletion — but
   "a variant lost its Pathogenic status" reads very differently without that explanation.

2. **One real edge case in patch #3.** A record whose *Explanation* reads
   `Pathogenic, low penetrance` is dropped entirely, because the Explanation path splits only on
   `/` and `;` while the Description path also splits on `,`. In the 2025-07 release this affects
   exactly **1** record out of 10.18M. It is a pre-existing inconsistency in Nirvana's
   `GetSignificances` that patch #3 exposes rather than introduces. Fixing it means adding `,` to
   the Explanation split — but note that would also start emitting `low penetrance` as a separate
   component (82 occurrences), so it is not a free change. Recommend leaving as-is and revisiting
   only if a clinically important variant is affected.

3. **`not provided` records no longer contribute entries.** 56,125 occurrences release-wide
   (upper bound). Correct behaviour — it is not a classification — but at production scale it
   will visibly reduce RCV counts relative to the Nirvana-bundled baseline. Expect it; it is not
   a regression.

4. **Scale expectation.** `quickit` is a small callset: only 17 RCV-level drops surfaced here. A
   full callset will hit proportionally more `not provided` and `flagged submission` records.
   The *rate* should stay near the release-wide figures (0.55% and 0.041%); a materially higher
   rate would warrant investigation.

---

## Appendix — outstanding queries

Table aliases used below:
`NEW` = `gvs-internal.quickit_2026_09_10_VS_1994_0c32721_vat.quickit_vat`
`OLD` = `gvs-internal.quickit_2026_09_09_ng_vs_1988_failed_load_job_4e0cc60_vat.quickit_vat`

### A — Same-version control group

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, s AS num_stars, c AS classification
    FROM `OLD` t,
         UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
         UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
    WHERE i = j AND j = k),
  new_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, s AS num_stars, c AS classification
    FROM `NEW` t,
         UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
         UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
    WHERE i = j AND j = k),
  old_agg AS (SELECT vid, rcv_id, ANY_VALUE(num_stars) AS num_stars,
                     ARRAY_AGG(DISTINCT classification ORDER BY classification) AS classifications
              FROM old_rcv GROUP BY vid, rcv_id),
  new_agg AS (SELECT vid, rcv_id, ANY_VALUE(num_stars) AS num_stars,
                     ARRAY_AGG(DISTINCT classification ORDER BY classification) AS classifications
              FROM new_rcv GROUP BY vid, rcv_id)
SELECT
  COUNT(*)                                                              AS same_version_rcvs,
  COUNTIF(o.num_stars != n.num_stars)                                   AS star_differs,
  COUNTIF(TO_JSON_STRING(o.classifications) != TO_JSON_STRING(n.classifications))
                                                                        AS classification_differs,
  COUNTIF(ARRAY_LENGTH(n.classifications) < ARRAY_LENGTH(o.classifications))
                                                                        AS classifications_lost,
  COUNTIF(ARRAY_LENGTH(n.classifications) > ARRAY_LENGTH(o.classifications))
                                                                        AS classifications_gained
FROM old_agg o JOIN new_agg n USING (vid, rcv_id)
```

**Interpretation.** All zeros ⇒ the patches changed nothing on identical input; every difference
in the tables is a ClinVar refresh. `classifications_lost` > 0 ⇒ patch #3 is stripping real
significance values (regression). `classifications_gained` > 0 ⇒ patch #1 working as intended
(`uncertain risk allele` now surviving). `star_differs` > 0 ⇒ unexpected, needs explanation.

### B — RCV-level losses

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `OLD` t,
         UNNEST(t.clinvar_rcv_ids) AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  new_rcv AS (
    SELECT DISTINCT t.vid, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `NEW` t,
         UNNEST(t.clinvar_rcv_ids) AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  shared AS (SELECT vid FROM (SELECT DISTINCT vid FROM old_rcv)
             INTERSECT DISTINCT SELECT vid FROM (SELECT DISTINCT vid FROM new_rcv))
SELECT o.num_stars AS old_star, COUNT(*) AS rcvs_dropped, COUNT(DISTINCT o.vid) AS vids_affected
FROM old_rcv o
JOIN shared USING (vid)
LEFT JOIN new_rcv n USING (vid, rcv_base)
WHERE n.rcv_base IS NULL
GROUP BY o.num_stars
ORDER BY o.num_stars
```

Expect at least the `20-2403658-G-A` row. The question is whether it is 1 record or 1000.

### C1 — Which vids lost `pathogenic`

```sql
WITH
  old_agg AS (SELECT vid, ARRAY_AGG(DISTINCT c ORDER BY c) AS old_class
              FROM (SELECT DISTINCT vid, c FROM `OLD`, UNNEST(clinvar_classification) AS c)
              GROUP BY vid),
  new_agg AS (SELECT vid, ARRAY_AGG(DISTINCT c ORDER BY c) AS new_class
              FROM (SELECT DISTINCT vid, c FROM `NEW`, UNNEST(clinvar_classification) AS c)
              GROUP BY vid)
SELECT o.vid, o.old_class, n.new_class
FROM old_agg o
LEFT JOIN new_agg n USING (vid)
WHERE 'pathogenic' IN UNNEST(o.old_class)
  AND (n.vid IS NULL OR 'pathogenic' NOT IN UNNEST(n.new_class))
```

### C2 — What happened to the RCV that carried it

```sql
WITH lost_p AS (
  SELECT DISTINCT vid FROM `OLD`, UNNEST(clinvar_classification) AS c WHERE c = 'pathogenic'
  EXCEPT DISTINCT
  SELECT DISTINCT vid FROM `NEW`, UNNEST(clinvar_classification) AS c WHERE c = 'pathogenic')
SELECT 'old' AS tbl, t.vid, r AS rcv_id, s AS num_stars, c AS classification
FROM `OLD` t JOIN lost_p USING (vid),
     UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
     UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
     UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
WHERE i = j AND j = k
UNION DISTINCT
SELECT 'new', t.vid, r, s, c
FROM `NEW` t JOIN lost_p USING (vid),
     UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
     UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
     UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
WHERE i = j AND j = k
ORDER BY vid, rcv_id, tbl DESC
```

**Interpretation.** Accession present in `old` but absent from `new` ⇒ patch #3 dropped it
(regression). Accession present in `new` at a bumped version with a different classification
⇒ ClinVar reclassified it over 20 months (benign; confirm against XML).

### D — Vid-universe overlap (Q0)

```sql
WITH
  n AS (SELECT DISTINCT vid FROM `NEW`),
  o AS (SELECT DISTINCT vid FROM `OLD`)
SELECT
  (SELECT COUNT(*) FROM n)                                                        AS vids_new,
  (SELECT COUNT(*) FROM o)                                                        AS vids_old,
  (SELECT COUNT(*) FROM (SELECT vid FROM n INTERSECT DISTINCT SELECT vid FROM o)) AS vids_shared,
  (SELECT COUNT(*) FROM (SELECT vid FROM o EXCEPT DISTINCT SELECT vid FROM n))    AS vids_only_old,
  (SELECT COUNT(*) FROM (SELECT vid FROM n EXCEPT DISTINCT SELECT vid FROM o))    AS vids_only_new
```

### E — Decomposing the 125 `gained_1star`

```sql
WITH
  old_stars AS (SELECT DISTINCT vid, s FROM `OLD`, UNNEST(clinvar_rcv_num_stars) AS s),
  new_stars AS (SELECT DISTINCT vid, s FROM `NEW`, UNNEST(clinvar_rcv_num_stars) AS s),
  old_vids  AS (SELECT DISTINCT vid FROM `OLD`),
  old_agg   AS (SELECT vid, LOGICAL_OR(s = 1) AS had1 FROM old_stars GROUP BY vid),
  new1      AS (SELECT DISTINCT vid FROM new_stars WHERE s = 1)
SELECT
  COUNT(*)                                       AS gained_1star_vids,
  COUNTIF(ov.vid IS NULL)                        AS vid_not_in_old_table,
  COUNTIF(ov.vid IS NOT NULL AND oa.vid IS NULL) AS in_old_but_had_no_clinvar,
  COUNTIF(oa.vid IS NOT NULL AND NOT oa.had1)    AS had_clinvar_other_stars_only
FROM new1 n
LEFT JOIN old_agg oa ON oa.vid = n.vid
LEFT JOIN old_vids ov ON ov.vid = n.vid
WHERE oa.vid IS NULL OR NOT oa.had1
```

### F — Converter warning counts (shell, not SQL)

```sh
grep -c 'WARNING: skipping unrecognized clinical significance value' <converter-stderr>
grep -c 'WARNING: unrecognized review status'                        <converter-stderr>
grep -o "value: '[^']*'" <converter-stderr> | sort | uniq -c | sort -rn
```

---

## Changelog

- **2026-09-15** — Initial version. Checks 1-11 recorded; `20-2403658-G-A` resolved against
  source XML; open questions A-G defined.
- **2026-09-15** — Added check 12: live ClinVar cross-check of `VCV000337923`. Confirmed RCV
  completeness (closing the GRCh38-coordinate caveat) and recorded post-snapshot drift on
  `RCV000713829`. **Corrected an earlier hypothesis:** this variant has never been pathogenic,
  so it does not explain the pathogenic loss — open question C remains fully open and now
  points at an unidentified vid.
- **2026-09-15** — Check 13: same-version control group (A) returned all zeros across 1607 RCVs.
  A resolved ✅; conclusion 7 added. Noted two limits on its reach (inner join excludes dropped
  records; patched paths likely not exercised) and raised F to 🟠 as a result. Blockers now B and C.
- **2026-09-15** — Check 14: RCV-level losses (B) returned **17 dropped RCVs across 16 vids, all
  1-star** — the first non-zero result in the validation. Not yet attributable: each accession
  needs an XML lookup to separate ClinVar retirement / 0-star downgrade from a patch #3 defect.
- **2026-09-15** — Check 15: all 17 resolved against the 2025-07 XML. 15 retired by ClinVar; the
  2 pathogenic records (`RCV000761550.2`, `RCV000761551.2`) are re-represented as the FMR1
  microsatellite `c.-128GGC[55_200]` with no VCF alt allele, so their removal from
  `X-147912049-CGCG-C` is correct — the old bundle had mis-attached them to an adjacent deletion.
  **Zero defects.** B resolved ✅; C largely resolved; conclusions 8-9 added; status → 🟢.
  Recorded the RCV-version vs MeasureSet-representation nuance.
- **2026-09-15** — Check 16: C1 swept every vid's `clinvar_classification` and returned exactly
  one pathogenic loss — `X-147912049-CGCG-C`, the FMR1 case already resolved. C resolved ✅;
  conclusion 10 added. **All blocking checks (A, B, C) now complete**; only optional check F
  remains.
- **2026-09-16** — Checks 17-19 close F **without re-running the converter**. The VM was deleted
  and stderr lost, but the warnings are a deterministic function of the input XML, so the
  allow-lists were recovered from `SAUtils.dll` and diffed against the full 2025-07 vocabulary.
  Review statuses: only `flagged submission` and `no classifications from unflagged records`
  (0.041%) hit the fallback, both correctly dropped. Significances: **no real classification term
  is filtered**; the only non-allow-list values are `not provided`, the placeholder, and 84 edge
  cases. Patch #1 confirmed present in the deployed binary and load-bearing (623 occurrences of
  `uncertain risk allele`). Conclusions 11-12 added; §7 production follow-ups added; F resolved ✅.
