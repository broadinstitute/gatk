# ClinVar Conversion Validation

**Question being answered:** does the process that produced the new ClinVar files introduce any
change that is *not* fully explained by changes in the source data itself?

**Status:** 🟢 **No unexplained change found.** All blocking checks are closed. Every observed
difference between the two VAT tables is attributable to a ~20-month ClinVar release gap. The
patches introduced no annotation change beyond discarding values that were never clinical
classifications.

**Last updated:** 2026-10-02 · **Reusable:** §8 and Appendix A are written to be re-run for the
next ClinVar upgrade.

---

## 1. What is being compared

|                              | New                                                                   | Old (baseline)                                         |
|------------------------------|-----------------------------------------------------------------------|--------------------------------------------------------|
| Terra workspace              | `gvs-dev/GVS Integration mcovarr`                                     | `gp-dsp-gvs-operations-terra/GVS Integration - Nalini` |
| Terra submission             | `ab719aab-c2fb-43ea-ac1e-6abe31810a4f`                                | `552694c3-807b-44dd-bb2d-57eddbf15b4d`                 |
| GATK commit                  | `0c32721` (VS-1994)                                                   | `4e0cc60` (VS-1988)                                    |
| ClinVar source               | `ClinVar_2025-07.nsa`, built from `ClinVarFullRelease_2025-07.xml.gz` | Nirvana-bundled ClinVar                                |
| `clinvar_last_updated` range | 2022-04-23 → **2025-06-29** (2784 vids)                               | 2022-04-23 → **2023-10-28** (2671 vids)                |

> ⚠️ **The baseline is not a correctness reference.** It is a ~20-month-older ClinVar release
> parsed by unpatched code, built from the development branch `ng_vs_1988_failed_load_job`.
> Divergence from it is expected. Correctness must be established against the source XML, not
> against this table.

Both VATs were loaded into integration-test BigQuery datasets that have since expired. Each
submission's VAT survives as the TSV it wrote, `vat_complete.bgz.tsv.gz` from `MergeVatTSVs`. The
baseline submission is marked failed, but only because of a VCF-side cost assertion; its VAT
steps completed. Both TSVs were reproduced byte for byte by the later runs in §9.

```
# New
gs://fc-358fc8a2-f84e-4d45-846a-a533d08f6103/submissions/ab719aab-c2fb-43ea-ac1e-6abe31810a4f/GvsQuickstartIntegration/d0ccc0c1-2d65-4914-8ea6-c63f4bd3d5e0/call-GvsQuickstartVATIntegration/GvsQuickstartVATIntegration/29d6f1d9-b8bd-4a99-8473-e34110fb84c7/call-CreateVATFromVDS/GvsCreateVATfromVDS/ccb39e44-1d27-4939-8856-7694777a23c6/call-GvsCreateVATFilesFromBigQuery/GvsCreateVATFilesFromBigQuery/9036451e-e0a4-441d-8336-85531b2aee24/call-MergeVatTSVs/vat_complete.bgz.tsv.gz
# Old (baseline)
gs://fc-22843ae9-4dbd-4751-ba8a-10fccdb797cd/submissions/552694c3-807b-44dd-bb2d-57eddbf15b4d/GvsQuickstartIntegration/79c6b053-c05d-4d13-a174-21671b84f327/call-GvsQuickstartVATIntegration/GvsQuickstartVATIntegration/5c25560f-f074-448d-9b95-75e5d21af22a/call-CreateVATFromVDS/GvsCreateVATfromVDS/9cd05997-641e-4781-9df7-2473ba552d1f/call-GvsCreateVATFilesFromBigQuery/GvsCreateVATFilesFromBigQuery/0ac5d4b1-469a-487c-810d-5a652e1c2c8b/call-MergeVatTSVs/vat_complete.bgz.tsv.gz
```

The TSVs are gzipped and unsorted; compare them with `gunzip -c <tsv> | LC_ALL=C sort`.

### Patches under test — [`clinvar_patches.diff`](clinvar_patches.diff)

Build and apply instructions are in [`clinvar_patches.README.md`](clinvar_patches.README.md).

| # | File                               | Change                                                                                          |
|---|------------------------------------|-------------------------------------------------------------------------------------------------|
| 1 | `ClinVarCommon.cs`                 | Added `uncertain risk allele` to `ValidPathogenicity`                                           |
| 2 | `ClinVarVariationReader.cs`        | `ResolveReviewStatus()` falls back to `no_assertion` instead of throwing (**VCV path only**)    |
| 3 | `ClinVarCommon.GetSignificances()` | Filters significances to `ValidPathogenicity`, skipping unrecognized values instead of throwing |
| 4 | `NsaWriter.cs`                     | Progress logging only — no functional impact                                                    |

Patches #2 and #3 exist because NCBI's germline/somatic classification split now emits placeholder
text (e.g. `no classifications from unflagged records`) in fields that previously held only
controlled vocabulary. Unpatched code threw and aborted the run.

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
- **All six ClinVar columns are written together** inside one `if len(clinvar_rcv_ids) > 0:` block
  (lines 271-294), so they are populated atomically per variant. A single column cannot regress
  independently of the others.

### Schema facts

- The VAT is **one row per `(vid, transcript)`** — all counts must be `COUNT(DISTINCT vid)` or
  pre-aggregated per vid.
- **5 of the 6 ClinVar columns are `REPEATED`.** `IS NOT NULL` is always true on an unset array;
  use `ARRAY_LENGTH(col) > 0`. Only `clinvar_last_updated` (scalar `DATE`) works with `IS NOT NULL`.
- Significance values are stored **lowercased**. Match exactly — `LIKE '%pathogenic%'` wrongly
  catches `conflicting interpretations of pathogenicity`.

---

## 3. Checks run

Every query and command is reproduced in **Appendix A**, keyed by the check number below.

| #  | Check                                    | Result                                                                                                                                           |
|----|------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| 1  | Star distribution                        | −702 one-star vids, +1023 two-star vids (net)                                                                                                    |
| 2  | Vid-level 2-star origin                  | 1337 new 2-star vids; 1024 newly so; **all 1024** were 1-star before; **0** from no-ClinVar                                                      |
| 3  | RCV-level 2-star origin                  | 1588 = 1131 promoted + 456 unchanged + 1 new accession; **0** demotions from 3/4-star; **0** vids absent from old                                |
| 4  | Lost 2-star vids                         | Exactly one: `20-2403658-G-A`, `[1,2] → [1]`                                                                                                     |
| 5  | 1-star decomposition                     | 1762 kept · 827 promoted away · **0 / 0 / 0** loss buckets · 125 gained                                                                          |
| 6  | Version check, promotions only           | **0** same-version · 1131 version-changed                                                                                                        |
| 7  | Version check, all shared RCVs           | 3649 shared · 1607 same · 2042 bumped · **all forward** · 0 backward                                                                             |
| 8  | ClinVar recency                          | old ≤ 2023-10-28 · new ≤ 2025-06-29 (~20-month gap)                                                                                              |
| 9  | Loader diff `4e0cc60`→`0c32721`          | Byte-identical                                                                                                                                   |
| 10 | Patch #2 reachability                    | Unreachable — loader accepts `RCV*` ids only                                                                                                     |
| 11 | XML ground truth, `20-2403658-G-A`       | 3 RCVs, all 1-star, all Benign — **new table correct**                                                                                           |
| 12 | Live ClinVar cross-check, `VCV000337923` | Same 3 RCVs (completeness confirmed); **never pathogenic**; post-snapshot drift `.8`→`.12`                                                       |
| 13 | **Same-version control group**           | 1607 RCVs · **0 / 0 / 0 / 0** — no star, classification, loss or gain difference on byte-identical input                                         |
| 14 | **RCV-level losses**                     | 17 RCVs dropped across 16 vids, all 1-star. No 2/3/4-star loss                                                                                   |
| 15 | **XML resolution of all 17**             | **15 retired by ClinVar · 2 re-represented as an FMR1 microsatellite · 0 defects**                                                               |
| 16 | **Pathogenic-loss sweep**                | **Exactly one vid**, `X-147912049-CGCG-C` — the FMR1 case. No other variant lost `pathogenic` by any route                                       |
| 17 | **Review-status vocabulary audit**       | 9 distinct values; 7 mapped. Only `flagged submission` + the placeholder (**4,218 / 10.18M = 0.041%**) hit the fallback — both correctly dropped |
| 18 | **Significance vocabulary audit**        | 77 distinct values; **no real classification term is filtered**                                                                                  |
| 19 | **Deployed-binary check**                | `uncertain risk allele` present in `SAUtils.dll` — patch #1 confirmed in the artifact that built the `.nsa`                                      |

### Detailed results

**Check 3 — RCV-level 2-star origin**

```
new_2star_rcvs  promoted_from_1star  unchanged_2star  demoted_3  demoted_4  rcv_new  vid_absent  distinct_vids_promoted
1588            1131                 456              0          0          1        0           1122
```

**Check 2 — vid-level 2-star origin**

```
new_2star_vids  also_2star_in_old  newly_2star_vids  newly_2star_and_was_1star  newly_2star_no_clinvar_in_old
1337            313                1024              1024                       0
```

**Check 5 — 1-star decomposition**

```
kept_1star  lost_1star_now_2star  lost_1star_other_star  lost_clinvar_vid_still_present  vid_gone_from_new_table  gained_1star
1762        827                   0                      0                               0                        125
```

**Check 7 — version comparison, all shared RCVs**

```
all_shared  same_version  version_changed  new_higher  new_LOWER  unparseable  example_id
3649        1607          2042             2042        0          0            RCV000408215.5
```

**Check 13 — same-version control group.** The direct patched-vs-unpatched measurement.

```
same_version_rcvs  star_differs  classification_differs  classifications_lost  classifications_gained
1607               0             0                       0                     0
```

For 1607 RCVs the two runs read the *same accession at the same version* — byte-identical record
content by ClinVar's own versioning contract. The patched converter produced identical stars and
identical classification sets on every one. **Patch #3 stripped no real significance value, and no
record parsed differently.**

Two limits on how far this result reaches, both since closed by other checks:

1. **It is an inner join.** Records present in `OLD` but dropped entirely from `NEW` cannot appear
   here — covered by checks 14-15, which found 17 such records and explained all of them.
2. **The patched code paths were probably not exercised here.** Patches #2/#3 fire only on values
   *outside* the controlled vocabulary, and a record containing such a value would have thrown in
   unpatched code, so it largely cannot be in the old table at all. "No difference" therefore
   proves no regression but does not by itself show the patches behaving correctly *when they
   fire* — covered by checks 17-18.

Patch #1's effect is not measured by this test; none of the control records carried
`uncertain risk allele`.

**Check 14 — RCV-level losses.** Accessions present in `OLD` but absent from `NEW`, shared vids only.

```
old_star  rcvs_dropped  vids_affected
1         17            16
```

The only non-zero loss result in the validation: 17 of 3649 shared RCVs (~0.5%), confined entirely
to 1-star records. The 17:

| vid                        | rcv_id             | stars | classification in OLD  |
|----------------------------|--------------------|-------|------------------------|
| 20-34945933-A-C            | RCV001516822.13    | 1     | benign                 |
| 20-59193946-G-A            | RCV002977546.1     | 1     | likely benign          |
| 20-6041987-A-G             | RCV003211151.1     | 1     | uncertain significance |
| X-120456678-T-A            | RCV000715317.2     | 1     | benign                 |
| X-136209863-C-T            | RCV000576703.1     | 1     | benign                 |
| X-141906091-C-A            | RCV003359679.1     | 1     | likely benign          |
| **X-147912049-CGCG-C**     | **RCV000761551.2** | 1     | **pathogenic**         |
| **X-147912049-CGCG-C**     | **RCV000761550.2** | 1     | **pathogenic**         |
| X-153905526-T-G            | RCV000576631.1     | 1     | benign                 |
| X-23000274-G-C             | RCV003346747.1     | 1     | uncertain significance |
| X-31679519-A-G             | RCV000576429.9     | 1     | benign                 |
| X-47846285-C-T             | RCV002316332.8     | 1     | benign                 |
| X-50607408-G-A             | RCV002316320.8     | 1     | benign                 |
| X-50607674-T-C             | RCV002460044.1     | 1     | benign                 |
| X-50607728-T-TTCC          | RCV002453406.1     | 1     | benign                 |
| X-50607758-C-CTGCTGCTGCTGT | RCV002453405.1     | 1     | benign                 |
| X-71123671-T-A             | RCV000460842.17    | 1     | benign                 |

Composition: 10 benign · 3 likely benign · 2 uncertain significance · **2 pathogenic**.

**Check 15 — resolution of all 17 against `ClinVarFullRelease_2025-07.xml.gz`**

| Outcome                                                                       | Count | Verdict                           |
|-------------------------------------------------------------------------------|-------|-----------------------------------|
| **Absent from the 2025-07 release** — accession retired or merged by ClinVar  | 15    | ✅ Benign attrition over 20 months |
| **Present, but re-represented as a microsatellite**                           | 2     | ✅ Correct drop — see §4           |
| Present with valid significance and non-zero review status at the same allele | **0** | ✅ No defect found                 |

**Check 16 — pathogenic-loss sweep.** Every vid whose `clinvar_classification` lost `pathogenic`:

```
vid                    old_class              new_class
X-147912049-CGCG-C     [benign,pathogenic]    [benign,uncertain significance]
```

A single row across the whole table, and it is the FMR1 case resolved in §4. This closes the one
route by which a high-stakes change could hide: a *surviving* RCV whose classification changed
drops no record and so is invisible to check 14, but this examines the derived
`clinvar_classification` field directly. The row also shows `uncertain significance` **gained** —
a new or updated RCV from the refresh — and the output ordering incidentally confirms the loader's
`significance_ordering` logic is intact.

**Checks 17-19 — vocabulary audit.** The converter VM was deleted and its stderr lost. Re-running
was unnecessary: the patch warnings are a deterministic function of the input XML, so the same
information — and more — was derived by enumerating the release's vocabulary and diffing it against
the allow-lists recovered from `build_output/SAUtils.dll` (UTF-16 literals in the .NET string
heap). This is strictly stronger than stderr: stderr records only what fired, the audit shows every
value that *could*.

*Review statuses — all 9 in the release:*

| Review status                                        | Count     | Mapped?            |
|------------------------------------------------------|----------:|--------------------|
| criteria provided, single submitter                  | 8,979,691 | ✅                  |
| no assertion criteria provided                       | 688,991   | ✅                  |
| criteria provided, multiple submitters, no conflicts | 345,497   | ✅                  |
| criteria provided, conflicting interpretations       | 57,516    | ✅                  |
| no assertion provided                                | 55,931    | ✅                  |
| reviewed by expert panel                             | 40,516    | ✅                  |
| practice guideline                                   | 5,006     | ✅                  |
| **flagged submission**                               | **2,706** | ❌ → `no_assertion` |
| **no classifications from unflagged records**        | **1,512** | ❌ → `no_assertion` |

Only 0.041% hit the fallback, and both are cases where dropping is correct — `flagged submission`
means ClinVar itself flagged the submission as unreliable; the other is the placeholder emitted
when nothing trustworthy remains. Neither is a clinical assertion.

`criteria provided, conflicting classifications` and `no classification provided` — which the
loader's Python dict knows but the C# mapping does not — **do not appear in the Full release at
all**, so that divergence is never exercised.

*Significances — 77 distinct values over 10,177,558 occurrences. Outside `ValidPathogenicity`:*

| Value                                       | Count  | Assessment                                                                                      |
|---------------------------------------------|-------:|-------------------------------------------------------------------------------------------------|
| `not provided`                              | 56,125 | Not a classification — ClinVar's "submitter gave none" marker                                   |
| `no classifications from unflagged records` | 1,515  | The placeholder that caused the original crash                                                  |
| `low penetrance`                            | 82     | Modifier split off `Pathogenic, low penetrance` by the comma; the `pathogenic` half is kept     |
| `no known pathogenicity`                    | 1      | Non-standard term                                                                               |
| `pathogenic, low penetrance`                | 1      | ⚠️ Explanation path splits only on `/` and `;`, so this stays whole and is dropped **entirely** |

**No real clinical-significance term is filtered.** `Uncertain risk allele` (623 occurrences)
survives *only* because of patch #1.

> *Caveat:* the scan covers all `<ClinicalSignificance>` blocks, including submitter-level
> `ClinVarAssertion` ones, while Nirvana reads a narrower subset. Counts are an **upper bound on
> exposure**, not record-drop counts. The vocabulary conclusion is unaffected.

### Arithmetic reconciliation — all consistent

| Identity                         | Check                                                          |
|----------------------------------|----------------------------------------------------------------|
| `125 − 827`                      | `= −702` ✅ matches the check-1 one-star delta                  |
| `827 + 197`                      | `= 1024` ✅ (197 = promoted but retained a 1-star RCV)          |
| `1024 + 98`                      | `= 1122` ✅ (98 = already had 2-star, gained another promotion) |
| `313 + 1024`                     | `= 1337` ✅                                                     |
| `1131 + 456 + 1`                 | `= 1588` ✅ every new 2-star RCV accounted for                  |
| `1337 − 1023 = 314`, `314 − 313` | `= 1` ✅ exactly one vid lost 2-star status                     |

---

## 4. The two anomalies, both resolved ✅

### 4.1 `20-2403658-G-A` — lost a 2-star rating

**Identity** — `VCV000337923`, gene **TGM6**, `NM_198994.3:c.1171G>A` (`p.Val391Met`).
Canonical SPDI `NC_000020.11:2403657:G:A` · GRCh38 `20:2403658` · GRCh37 `20:2384304`.

Ground truth from `ClinVarFullRelease_2025-07.xml.gz`, matched on GRCh38 `start=2403658`, `G>A`:

| RCV accession  | ReviewStatus                        | Stars | Description |
|----------------|-------------------------------------|-------|-------------|
| RCV000277947.7 | criteria provided, single submitter | 1     | Benign      |
| RCV000713829.8 | criteria provided, single submitter | 1     | Benign      |
| RCV004999340.1 | criteria provided, single submitter | 1     | Benign      |

**Verdict: the new table is correct; the old table was stale.** Every record in the snapshot is
1-star, so `new_stars = [1]` is right. The old `[1,2]` reflects the 2023 bundle, where one record
still carried `multiple submitters, no conflicts`. The sole apparent regression in the entire star
analysis is a *correction*.

Cross-checked against the live ClinVar page (fetched 2026-09-15): the same three RCVs, confirming
completeness and closing the GRCh38-coordinate-matching caveat. The variant has **never** been
classified Pathogenic — four Benign submissions plus one Uncertain significance, the latter
carrying `no assertion criteria provided` → 0 stars → correctly dropped by the loader. Note
`RCV000713829` has since advanced to `.12` at 2 stars; that is post-snapshot drift and a future
build should legitimately show this variant back at 2 stars.

### 4.2 `X-147912049-CGCG-C` — lost its Pathogenic call

Both pathogenic records survive in 2025-07 with `RecordStatus: current`,
`criteria provided, single submitter`, `Description: Pathogenic` — so neither the 0-star filter nor
patch #3 explains their absence. Their `MeasureSet` is what does:

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
deletion, and the new build correctly drops it.** A correctness improvement — though a clinically
consequential one, since a variant previously annotated Pathogenic no longer is.

> **Nuance worth remembering:** the RCV version tracks the *clinical assertion*, not the variant
> representation. A MeasureSet/VCV can be re-normalized without bumping the RCV version — which is
> how these records stayed at `.2` while their allele mapping changed. This does not undermine
> check 13, which compared stars and classifications (assertion-level attributes the RCV version
> does track), but **"same RCV version" must not be read as "same allele mapping."**

---

## 5. Conclusions

1. **No variant lost its ClinVar annotation.** All three loss buckets in check 5 are zero; the only
   route out of 1-star status was promotion to 2-star.
2. **No record was parsed differently.** Zero promotions occurred at an unchanged accession version
   — every star change accompanied an upstream ClinVar revision.
3. **The data is strictly fresher.** Zero backward version moves across 3649 shared RCVs.
4. **The star mapping is untouched** — loader byte-identical, and 1↔2 stars are distinct
   review-status strings.
5. **Patch #2 has no annotation impact** — provably unreachable from the VAT.
6. **On byte-identical input the patched converter produces identical output** — 1607 same-version
   RCVs, zero differences (check 13). This is the direct patched-vs-unpatched measurement.
7. **Every RCV-level loss is explained by source data** — all 17 resolved against the XML: 15
   retired by ClinVar, 2 re-represented as a microsatellite. Zero defects.
8. **Exactly one variant lost `pathogenic` table-wide**, and it is a correction (§4.2).
9. **The patches discard only non-classifications** — `not provided`, `flagged submission`, and
   NCBI's placeholder are all markers for "no trustworthy classification exists."
10. **Patch #1 is confirmed in the deployed binary and is load-bearing** — 623 occurrences of
    `uncertain risk allele` would otherwise be discarded.
11. **Both flagged anomalies were verified against source XML** and both turned out to be
    corrections, not regressions.
12. **The arithmetic closes** — independent queries reconcile exactly, with no remainder.

**No unexplained change was found at any level — vid, RCV, or classification.** Patch #3's
filtering is ruled out both for records that survived into the new table (check 13) and for every
record that did not (check 15), and the patches' behavior was audited directly against the source
release rather than inferred.

---

## 6. Optional checks not run

All blocking questions are closed. These three were defined early and never run; each is a
cross-check on numbers already reconciled by other means. Queries are in **Appendix B**.

|   | Item                                         | Why it was not needed                                                                                                                                                                                                                                                                                                               |
|---|----------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| D | Vid-universe overlap                         | Largely answered indirectly: check 3 returned `vid_absent_from_old = 0`, check 2 returned `newly_2star_no_clinvar_in_old = 0`, and check 5 returned `vid_gone_from_new_table = 0`. No variant appeared or vanished in any population examined. Worth running on a production-scale comparison, where partial loads are a real risk. |
| E | Decomposition of the 125 `gained_1star` vids | Gains are the expected direction for a newer release. Splits the 125 into "new to the table" vs "new ClinVar data for an existing variant" — informative, not load-bearing.                                                                                                                                                         |
| G | Per-column fill rates                        | **Structurally near-moot.** The loader writes all six ClinVar columns inside one `if len(clinvar_rcv_ids) > 0:` block, so they populate atomically per variant and cannot diverge independently. Only `clinvar_phenotype` can legitimately be empty while the others are set.                                                       |

---

## 7. Follow-ups for production

None affect the verdict on this dataset. All are worth carrying into a full-callset run.

1. **Tell downstream consumers about `X-147912049-CGCG-C`.** It moved from `[benign, pathogenic]`
   to `[benign, uncertain significance]`. The change is correct, but "a variant lost its Pathogenic
   status" reads very differently without the explanation.
2. **One real edge case in patch #3.** A record whose *Explanation* reads
   `Pathogenic, low penetrance` is dropped entirely, because that path splits on `/` and `;` but
   not `,`. Exactly **1** record in 10.18M. It is a pre-existing inconsistency in Nirvana's
   `GetSignificances`, not something the patches introduced. The obvious fix — adding `,` to the
   Explanation split — would also start emitting `low penetrance` as a separate component (82
   occurrences), so it is not free. Recommend leaving as-is unless a clinically important variant
   is affected.
3. **`not provided` records no longer contribute entries.** 56,125 occurrences release-wide (upper
   bound). Correct behavior, but at production scale it will visibly reduce RCV counts relative to
   the Nirvana-bundled baseline. Expect it; it is not a regression.
4. **Scale expectation.** `quickit` is a small callset — only 17 RCV-level drops surfaced here. A
   full callset will hit proportionally more `not provided` and `flagged submission` records. The
   *rate* should stay near the release-wide figures (0.55% and 0.041%); a materially higher rate
   warrants investigation.

---

## 8. How to reuse this for the next upgrade

### Suggested order

Each step either frames the comparison or closes a specific hiding place for an unexplained change.

| Order | Check(s)                       | Purpose                                                                                                                                                     |
|-------|--------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1     | **D** (Appendix B), then **8** | Establish comparability and size the release gap *before* interpreting anything. If the variant universes differ materially, every later number is suspect. |
| 2     | **1**                          | Frame the question — where did stars move?                                                                                                                  |
| 3     | **2, 3**                       | Attribute the gains. Vid-level first, then record-level; only the record-level query distinguishes "existing record re-graded" from "new record appeared."  |
| 4     | **6, 7**                       | **The pivotal step.** Compare accession *versions*. A change at an unchanged version implicates your code; a change with a version bump implicates ClinVar. |
| 5     | **5**                          | Decompose the losses and confirm the arithmetic closes against step 2.                                                                                      |
| 6     | **13**                         | The control group — the only true A/B test available.                                                                                                       |
| 7     | **14**, then **15**            | Find record-level losses, then resolve every one against the source XML. Do not stop at the count.                                                          |
| 8     | **16**                         | Sweep the derived `clinvar_classification` for high-stakes losses that drop no record.                                                                      |
| 9     | **17, 18, 19**                 | Audit the release vocabulary against the allow-lists in the deployed binary.                                                                                |

### Principles that made this work

- **The old table is never ground truth.** It is a different release parsed by different code. Only
  the source XML can settle correctness. Both anomalies here looked like regressions until the XML
  showed they were corrections.
- **Same accession + same version is the only clean A/B.** ClinVar's versioning contract makes those
  records byte-identical input, so any output difference is attributable to your code alone. This is
  the single most valuable check in the set — but remember §4.2's nuance: the version tracks the
  *assertion*, not the allele representation.
- **Resolve every non-zero individually.** 17 losses sounds benign and 2 of them were pathogenic.
  Counts hide the cases that matter.
- **Derive, don't re-run.** Patch warnings are a deterministic function of the input XML. When the
  run logs are gone, auditing the release vocabulary against the binary's allow-lists is cheaper
  *and* more complete than reproducing the run.
- **Make the arithmetic close.** Independent queries that reconcile to the same totals are strong
  evidence against a missed population. A remainder means something is unaccounted for.

### Traps specific to this schema

- One row per `(vid, transcript)` — always `COUNT(DISTINCT vid)` or collapse with `QUALIFY`.
- `SELECT DISTINCT` fails on `ARRAY` columns; use `QUALIFY ROW_NUMBER() OVER (PARTITION BY vid ...)`.
- Five of six ClinVar columns are `REPEATED` — `ARRAY_LENGTH(col) > 0`, never `IS NOT NULL`.
- Zip parallel arrays with `UNNEST(...) WITH OFFSET i` / `WITH OFFSET j` and `WHERE i = j`.
- Join RCVs on the **base accession** (`SPLIT(r,'.')[OFFSET(0)]`), keeping the full id to compare
  versions separately.
- Significances are lowercase; match exactly.

---

## 9. Reproduced through the VAT workflow

The checks above compare tables from the VS-1994 build. On 2026-09-30, two `GvsQuickstartIntegration`
runs of `vs_2029_clinvar_clean` used reference disks and differed only in
`use_manual_clinvar_update`. Checks 1, 7, 8, 13, 14 and 16 were re-run with the flag-off VAT as
`{{OLD}}` and the flag-on VAT as `{{NEW}}`, and they reproduced every number in §3. That covers the
star deltas, version moves, `clinvar_last_updated` range, control group, the 17 dropped RCVs and
the single pathogenic loss. All non-ClinVar columns were identical between the two runs. The
flag-off run passed the integration test's exactness assertions, and the flag-on run failed them,
as expected until the truth data is regenerated.

The two runs are submissions `6ee83e18-2181-4e54-ba0b-22b9c47d0da0` (flag off) and
`af0f9e5d-dd99-499f-b8ae-77ba3ae260f0` (flag on) in `gvs-dev/GVS Integration mcovarr`.

A third run on 2026-10-01, `c8d31878-9245-43e8-986c-1e3c47f0c7e9`, set the flag and `use_reference_disk = false`, the configuration of the
VS-1994 build. Its sorted VAT TSV was byte-for-byte identical to the new table's TSV (§1). The
flag-off run's sorted TSV was likewise byte-for-byte identical to the baseline's, so both sides of
the original comparison have been reproduced exactly through the VS-2029 WDL.

The three runs' TSVs:

```
# Flag off, reference disks (6ee83e18)
gs://fc-358fc8a2-f84e-4d45-846a-a533d08f6103/submissions/6ee83e18-2181-4e54-ba0b-22b9c47d0da0/GvsQuickstartIntegration/9f1fda0d-4886-4e25-8fa8-ea22852466c4/call-GvsQuickstartVATIntegration/GvsQuickstartVATIntegration/e3ca8dd0-b741-405e-878b-5c56ca33ced2/call-CreateVATFromVDS/GvsCreateVATfromVDS/325d1f80-614d-4499-a1d4-f48aacc87c4e/call-GvsCreateVATFilesFromBigQuery/GvsCreateVATFilesFromBigQuery/a618183b-7469-49b3-b8e9-85ca8a7beeb9/call-MergeVatTSVs/vat_complete.bgz.tsv.gz
# Flag on, reference disks (af0f9e5d)
gs://fc-358fc8a2-f84e-4d45-846a-a533d08f6103/submissions/af0f9e5d-dd99-499f-b8ae-77ba3ae260f0/GvsQuickstartIntegration/65c8b0c2-7249-4069-928c-a058229c4b9d/call-GvsQuickstartVATIntegration/GvsQuickstartVATIntegration/d7b85b65-aea2-45ba-b693-293762594609/call-CreateVATFromVDS/GvsCreateVATfromVDS/37fb09e7-eb94-4dc9-80a7-6e43f6bfebf5/call-GvsCreateVATFilesFromBigQuery/GvsCreateVATFilesFromBigQuery/59dd24bd-21e3-4921-8065-dba030017baa/call-MergeVatTSVs/vat_complete.bgz.tsv.gz
# Flag on, no reference disks (c8d31878)
gs://fc-358fc8a2-f84e-4d45-846a-a533d08f6103/submissions/c8d31878-9245-43e8-986c-1e3c47f0c7e9/GvsQuickstartIntegration/0fcad72f-b28d-4820-aba5-95d6e3c89c9b/call-GvsQuickstartVATIntegration/GvsQuickstartVATIntegration/fc56b43c-3335-44dd-b777-61ebe26aac61/call-CreateVATFromVDS/GvsCreateVATfromVDS/83932f26-ca8d-4ad4-bc67-454a2223fbd6/call-GvsCreateVATFilesFromBigQuery/GvsCreateVATFilesFromBigQuery/a9a1f60a-ac94-45a0-8da4-307f136b8e55/call-MergeVatTSVs/vat_complete.bgz.tsv.gz
```

---

## Appendix A — queries and commands, by check

Substitute your two VAT tables throughout:

```
{{NEW}} = <project>.<dataset>.<table>   -- VAT built with the new ClinVar
{{OLD}} = <project>.<dataset>.<table>   -- VAT built with the previous ClinVar
```

Several queries share the same `old_rcv` / `new_rcv` CTE pair that zips the parallel RCV arrays;
it is repeated in each query so every block is independently runnable.

### Check 1 — Star distribution

```sql
WITH
  old_s AS (SELECT s AS num_stars, COUNT(DISTINCT vid) AS old_vids
            FROM `{{OLD}}`, UNNEST(clinvar_rcv_num_stars) AS s GROUP BY s),
  new_s AS (SELECT s AS num_stars, COUNT(DISTINCT vid) AS new_vids
            FROM `{{NEW}}`, UNNEST(clinvar_rcv_num_stars) AS s GROUP BY s)
SELECT num_stars,
       IFNULL(old_vids, 0) AS old_vids,
       IFNULL(new_vids, 0) AS new_vids,
       IFNULL(new_vids, 0) - IFNULL(old_vids, 0) AS delta
FROM old_s FULL OUTER JOIN new_s USING (num_stars)
ORDER BY num_stars
```

### Check 2 — Vid-level 2-star origin

```sql
WITH
  old_v AS (SELECT DISTINCT vid, s AS num_stars
            FROM `{{OLD}}`, UNNEST(clinvar_rcv_num_stars) AS s),
  new_v AS (SELECT DISTINCT vid, s AS num_stars
            FROM `{{NEW}}`, UNNEST(clinvar_rcv_num_stars) AS s),
  old_agg AS (SELECT vid,
                     LOGICAL_OR(num_stars = 1) AS had_1star,
                     LOGICAL_OR(num_stars = 2) AS had_2star
              FROM old_v GROUP BY vid),
  new_2 AS (SELECT DISTINCT vid FROM new_v WHERE num_stars = 2)
SELECT
  COUNT(*)                                                               AS new_2star_vids,
  COUNTIF(IFNULL(o.had_2star, FALSE))                                    AS also_2star_in_old,
  COUNTIF(NOT IFNULL(o.had_2star, FALSE))                                AS newly_2star_vids,
  COUNTIF(NOT IFNULL(o.had_2star, FALSE) AND IFNULL(o.had_1star, FALSE)) AS newly_2star_and_was_1star,
  COUNTIF(NOT IFNULL(o.had_2star, FALSE) AND o.vid IS NULL)              AS newly_2star_no_clinvar_in_old
FROM new_2 n
LEFT JOIN old_agg o USING (vid)
```

### Check 3 — RCV-level 2-star origin

Distinguishes "the same accession was re-graded" from "a new accession appeared" — check 2 cannot.

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{OLD}}` t,
         UNNEST(t.clinvar_rcv_ids)       AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  new_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{NEW}}` t,
         UNNEST(t.clinvar_rcv_ids)       AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  old_vids AS (SELECT DISTINCT vid FROM `{{OLD}}`)
SELECT
  COUNT(*)                                           AS new_2star_rcvs,
  COUNTIF(o.num_stars = 1)                           AS promoted_from_1star,
  COUNTIF(o.num_stars = 2)                           AS unchanged_2star,
  COUNTIF(o.num_stars = 3)                           AS demoted_from_3star,
  COUNTIF(o.num_stars = 4)                           AS demoted_from_4star,
  COUNTIF(o.rcv_base IS NULL AND ov.vid IS NOT NULL) AS rcv_new_but_vid_in_old,
  COUNTIF(ov.vid IS NULL)                            AS vid_absent_from_old,
  COUNT(DISTINCT IF(o.num_stars = 1, n.vid, NULL))   AS distinct_vids_promoted_1to2
FROM new_rcv n
LEFT JOIN old_rcv o USING (vid, rcv_base)
LEFT JOIN old_vids ov ON ov.vid = n.vid
WHERE n.num_stars = 2
```

Optional guard — should return 0, confirming one star value per accession per vid:

```sql
SELECT COUNT(*) AS dup_star_conflicts FROM (
  SELECT vid, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, COUNT(DISTINCT s) AS star_values
  FROM `{{NEW}}`,
       UNNEST(clinvar_rcv_ids) AS r WITH OFFSET i,
       UNNEST(clinvar_rcv_num_stars) AS s WITH OFFSET j
  WHERE i = j
  GROUP BY 1, 2 HAVING COUNT(DISTINCT s) > 1)
```

### Check 4 — Vids that lost a star level

```sql
WITH
  old_agg AS (SELECT vid, ARRAY_AGG(DISTINCT s ORDER BY s) AS old_stars
              FROM (SELECT DISTINCT vid, s FROM `{{OLD}}`, UNNEST(clinvar_rcv_num_stars) AS s)
              GROUP BY vid),
  new_agg AS (SELECT vid, ARRAY_AGG(DISTINCT s ORDER BY s) AS new_stars
              FROM (SELECT DISTINCT vid, s FROM `{{NEW}}`, UNNEST(clinvar_rcv_num_stars) AS s)
              GROUP BY vid)
SELECT o.vid, o.old_stars, n.new_stars
FROM old_agg o
LEFT JOIN new_agg n USING (vid)
WHERE 2 IN UNNEST(o.old_stars)
  AND (n.vid IS NULL OR 2 NOT IN UNNEST(n.new_stars))
```

### Check 5 — 1-star decomposition

The result must reconcile: `gained_1star − (the four loss buckets)` equals the check-1 delta.

```sql
WITH
  old_agg AS (SELECT vid, LOGICAL_OR(s = 1) AS had1, LOGICAL_OR(s = 2) AS had2
              FROM (SELECT DISTINCT vid, s FROM `{{OLD}}`, UNNEST(clinvar_rcv_num_stars) AS s)
              GROUP BY vid),
  new_agg AS (SELECT vid, LOGICAL_OR(s = 1) AS has1, LOGICAL_OR(s = 2) AS has2
              FROM (SELECT DISTINCT vid, s FROM `{{NEW}}`, UNNEST(clinvar_rcv_num_stars) AS s)
              GROUP BY vid),
  new_vids AS (SELECT DISTINCT vid FROM `{{NEW}}`)
SELECT
  COUNTIF(o.had1 AND IFNULL(n.has1, FALSE))                             AS kept_1star,
  COUNTIF(o.had1 AND NOT IFNULL(n.has1, FALSE)
                 AND IFNULL(n.has2, FALSE))                             AS lost_1star_now_2star,
  COUNTIF(o.had1 AND NOT IFNULL(n.has1, FALSE)
                 AND NOT IFNULL(n.has2, FALSE) AND n.vid IS NOT NULL)   AS lost_1star_other_star,
  COUNTIF(o.had1 AND n.vid IS NULL AND nv.vid IS NOT NULL)              AS lost_clinvar_vid_still_present,
  COUNTIF(o.had1 AND nv.vid IS NULL)                                    AS vid_gone_from_new_table,
  COUNTIF(NOT IFNULL(o.had1, FALSE) AND IFNULL(n.has1, FALSE))          AS gained_1star
FROM old_agg o
FULL OUTER JOIN new_agg n USING (vid)
LEFT JOIN new_vids nv ON nv.vid = IFNULL(o.vid, n.vid)
```

### Check 6 — Accession versions on re-graded records

**The pivotal query.** `same_accession_version > 0` means identical input produced different output
— that implicates your code and must be explained.

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{OLD}}` t,
         UNNEST(t.clinvar_rcv_ids)       AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  new_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{NEW}}` t,
         UNNEST(t.clinvar_rcv_ids)       AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j)
SELECT
  COUNTIF(o.rcv_id =  n.rcv_id) AS same_accession_version,
  COUNTIF(o.rcv_id != n.rcv_id) AS version_changed
FROM new_rcv n JOIN old_rcv o USING (vid, rcv_base)
WHERE n.num_stars = 2 AND o.num_stars = 1
```

### Check 7 — Accession versions, all shared records

`new_version_LOWER > 0` would mean the new build read *older* data — alarming, and it inverts the
meaning of every "gain."

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{OLD}}` t,
         UNNEST(t.clinvar_rcv_ids)       AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  new_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{NEW}}` t,
         UNNEST(t.clinvar_rcv_ids)       AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  pairs AS (
    SELECT o.rcv_id AS old_rcv_id, n.rcv_id AS new_rcv_id,
           SAFE_CAST(SPLIT(o.rcv_id, '.')[SAFE_OFFSET(1)] AS INT64) AS old_v,
           SAFE_CAST(SPLIT(n.rcv_id, '.')[SAFE_OFFSET(1)] AS INT64) AS new_v
    FROM new_rcv n JOIN old_rcv o USING (vid, rcv_base))
SELECT
  COUNT(*)                                AS all_shared_rcvs,
  COUNTIF(old_rcv_id  = new_rcv_id)       AS same_version_all,
  COUNTIF(old_rcv_id != new_rcv_id)       AS version_changed_all,
  COUNTIF(new_v > old_v)                  AS new_version_higher,
  COUNTIF(new_v < old_v)                  AS new_version_LOWER,
  COUNTIF(new_v IS NULL OR old_v IS NULL) AS unparseable_version,
  ANY_VALUE(old_rcv_id)                   AS example_old_id,
  ANY_VALUE(new_rcv_id)                   AS example_new_id
FROM pairs
```

### Check 8 — ClinVar recency

Dates the two releases and sizes the gap.

```sql
SELECT 'old' AS tbl, MIN(clinvar_last_updated) AS earliest,
       MAX(clinvar_last_updated) AS latest, COUNT(DISTINCT vid) AS vids_with_date
FROM `{{OLD}}` WHERE clinvar_last_updated IS NOT NULL
UNION ALL
SELECT 'new', MIN(clinvar_last_updated), MAX(clinvar_last_updated), COUNT(DISTINCT vid)
FROM `{{NEW}}` WHERE clinvar_last_updated IS NOT NULL
```

### Checks 9-10 — Code inspection

```sh
# Check 9 — is the loader identical between the two runs? Empty output ⇒ yes.
cd /path/to/gatk
git diff <OLD_COMMIT> <NEW_COMMIT> \
  -- scripts/variantstore/scripts/create_vt_bqloadjson_from_annotations.py

# Check 10 — which annotation ids does the loader accept, and what is dropped?
grep -n 'RCV\|num_stars != 0\|ReviewStatusNameMapping\|vat_clinvar_review_status_dictionary' \
  scripts/variantstore/scripts/create_vt_bqloadjson_from_annotations.py
```

### Check 11 — XML ground truth for one variant

Two passes: locate the variant, then extract the enclosing records. Each pass decompresses the
whole file (several minutes); run in the background.

```sh
XML=ClinVarFullRelease_2025-07.xml.gz
POS=2403658

# Pass 1 — find the GRCh38 SequenceLocation lines and note their line numbers
gzcat $XML | grep -n "start=\"$POS\"" > pos_hits.txt

# Pass 2 — dump a window around each hit (ReviewStatus precedes SequenceLocation,
# so weight the window backwards), then read off accession / status / description
gzcat $XML | awk '
  (NR>=683134264 && NR<=683136064) { print NR": "$0 }
  NR>683136064 { exit }' > window.txt
grep -E 'ClinVarAccession Acc="RCV|<ReviewStatus>|<Description>|SequenceLocation' window.txt
```

### Check 12 — Live ClinVar cross-check

Confirms RCV completeness and surfaces submission history, including any past Pathogenic call:

```
https://www.ncbi.nlm.nih.gov/clinvar/variation/<VCV numeric id>/
```

The live page reflects ClinVar *today*, which may post-date your snapshot. Where they disagree the
XML is authoritative for validating the pipeline — it is what the converter read.

### Check 13 — Same-version control group

**The most valuable single check.** Joining on the full `rcv_id` (version included) restricts to
records whose content is identical by ClinVar's versioning contract.

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, s AS num_stars, c AS classification
    FROM `{{OLD}}` t,
         UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
         UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
    WHERE i = j AND j = k),
  new_rcv AS (
    SELECT DISTINCT t.vid, r AS rcv_id, s AS num_stars, c AS classification
    FROM `{{NEW}}` t,
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

To list offenders instead of counting them, swap the `SELECT` for:

```sql
SELECT o.vid, o.rcv_id, o.num_stars AS old_star, n.num_stars AS new_star,
       o.classifications AS old_class, n.classifications AS new_class
FROM old_agg o JOIN new_agg n USING (vid, rcv_id)
WHERE o.num_stars != n.num_stars
   OR TO_JSON_STRING(o.classifications) != TO_JSON_STRING(n.classifications)
LIMIT 100
```

### Check 14 — RCV-level losses

Vid-level zeros do not prove record-level zeros: a variant with four RCVs can lose two and still
register as "has ClinVar."

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{OLD}}` t,
         UNNEST(t.clinvar_rcv_ids) AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars) AS s WITH OFFSET j
    WHERE i = j),
  new_rcv AS (
    SELECT DISTINCT t.vid, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, s AS num_stars
    FROM `{{NEW}}` t,
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

Then list them with what they carried — the `old_classifications` column is what tells you whether
anything high-stakes was lost:

```sql
WITH
  old_rcv AS (
    SELECT DISTINCT t.vid, SPLIT(r, '.')[OFFSET(0)] AS rcv_base, r AS rcv_id,
           s AS num_stars, c AS classification
    FROM `{{OLD}}` t,
         UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
         UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
         UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
    WHERE i = j AND j = k),
  new_base AS (
    SELECT DISTINCT t.vid, SPLIT(r, '.')[OFFSET(0)] AS rcv_base
    FROM `{{NEW}}` t, UNNEST(t.clinvar_rcv_ids) AS r),
  shared AS (
    SELECT vid FROM (SELECT DISTINCT vid FROM old_rcv)
    INTERSECT DISTINCT
    SELECT vid FROM (SELECT DISTINCT vid FROM new_base))
SELECT o.vid, o.rcv_id, o.num_stars,
       STRING_AGG(DISTINCT o.classification, '; ' ORDER BY o.classification) AS old_classifications
FROM old_rcv o
JOIN shared USING (vid)
LEFT JOIN new_base n USING (vid, rcv_base)
WHERE n.rcv_base IS NULL
GROUP BY o.vid, o.rcv_id, o.num_stars
ORDER BY o.vid
```

### Check 15 — Resolve dropped accessions against the XML

One pass for all accessions. **Absent from the output ⇒ ClinVar retired or merged it** (benign).

```sh
XML=ClinVarFullRelease_2025-07.xml.gz
ACC='RCV000761550|RCV000761551|RCV001516822'   # ...all dropped accessions, pipe-separated

gzcat $XML | grep -E -A 8 "ClinVarAccession Acc=\"($ACC)\"" > dropped_rcvs.txt
```

Each hit shows `RecordStatus`, `ReviewStatus` and `Description`. Interpretation:

| Finding                                                | Verdict                                                                     |
|--------------------------------------------------------|-----------------------------------------------------------------------------|
| Accession absent entirely                              | ClinVar retired/merged it — benign                                          |
| Present, review status maps to 0 stars                 | Correctly dropped by loader line 262 — benign                               |
| Present, valid significance and non-zero review status | **Investigate** — check the allele mapping before concluding it is a defect |

For the third case, pull the `MeasureSet` to see which variant the record now describes — a
re-normalized or microsatellite representation legitimately detaches it from a SNV/indel:

```sh
gzcat $XML | grep -E -A 400 'ClinVarAccession Acc="(RCV000761550|RCV000761551)" Version=' \
  | grep -E '<Measure |MeasureSet |<Name|ElementValue|SequenceLocation Assembly="GRCh38"|ClinVarAccession Acc="RCV'
```

### Check 16 — Pathogenic-loss sweep

Catches losses that drop no record — a surviving RCV whose classification changed.

```sql
WITH
  old_agg AS (SELECT vid, ARRAY_AGG(DISTINCT c ORDER BY c) AS old_class
              FROM (SELECT DISTINCT vid, c FROM `{{OLD}}`, UNNEST(clinvar_classification) AS c)
              GROUP BY vid),
  new_agg AS (SELECT vid, ARRAY_AGG(DISTINCT c ORDER BY c) AS new_class
              FROM (SELECT DISTINCT vid, c FROM `{{NEW}}`, UNNEST(clinvar_classification) AS c)
              GROUP BY vid)
SELECT o.vid, o.old_class, n.new_class
FROM old_agg o
LEFT JOIN new_agg n USING (vid)
WHERE 'pathogenic' IN UNNEST(o.old_class)
  AND (n.vid IS NULL OR 'pathogenic' NOT IN UNNEST(n.new_class))
```

Then dump the full RCV detail for whatever it returns, old and new side by side:

```sql
WITH lost_p AS (
  SELECT DISTINCT vid FROM `{{OLD}}`, UNNEST(clinvar_classification) AS c WHERE c = 'pathogenic'
  EXCEPT DISTINCT
  SELECT DISTINCT vid FROM `{{NEW}}`, UNNEST(clinvar_classification) AS c WHERE c = 'pathogenic')
SELECT 'old' AS tbl, t.vid, r AS rcv_id, s AS num_stars, c AS classification
FROM `{{OLD}}` t JOIN lost_p USING (vid),
     UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
     UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
     UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
WHERE i = j AND j = k
UNION DISTINCT
SELECT 'new', t.vid, r, s, c
FROM `{{NEW}}` t JOIN lost_p USING (vid),
     UNNEST(t.clinvar_rcv_ids)             AS r WITH OFFSET i,
     UNNEST(t.clinvar_rcv_num_stars)       AS s WITH OFFSET j,
     UNNEST(t.clinvar_rcv_classifications) AS c WITH OFFSET k
WHERE i = j AND j = k
ORDER BY vid, rcv_id, tbl DESC
```

Swap `'pathogenic'` for `'likely pathogenic'` to sweep that tier too.

### Checks 17-18 — Release vocabulary audit

One pass over the XML, scoped to `<ClinicalSignificance>` blocks. **Scoping matters:**
`<Description>` also appears in assay-method and observation elements, and an unscoped scan returns
thousands of irrelevant values.

```sh
gzcat ClinVarFullRelease_2025-07.xml.gz | awk '
 /<ClinicalSignificance/    { incs=1 }
 /<\/ClinicalSignificance>/ { incs=0 }
 incs && match($0, /<ReviewStatus>[^<]*<\/ReviewStatus>/) {
   v=substr($0,RSTART,RLENGTH); gsub(/<\/?ReviewStatus>/,"",v); rs[v]++ }
 incs && match($0, /<Description>[^<]*<\/Description>/) {
   v=substr($0,RSTART,RLENGTH); gsub(/<\/?Description>/,"",v); de[v]++ }
 incs && match($0, /<Explanation[^>]*>[^<]*<\/Explanation>/) {
   v=substr($0,RSTART,RLENGTH); sub(/<Explanation[^>]*>/,"",v); gsub(/<\/Explanation>/,"",v); ex[v]++ }
 END {
   for (k in rs) printf "REVIEWSTATUS\t%d\t%s\n", rs[k], k
   for (k in de) printf "DESCRIPTION\t%d\t%s\n",  de[k], k
   for (k in ex) printf "EXPLANATION\t%d\t%s\n",  ex[k], k
 }' > vocab.txt
```

Then diff against the allow-lists, replicating `GetSignificances`' own splitting rules — the
Description path splits on `/ , ;` while the Explanation path splits only on `/ ;` and strips
`(n)` counts:

```python
import re
VALID = {  # ValidPathogenicity, recovered from the binary — see check 19
 'benign','likely benign','pathogenic','likely pathogenic','uncertain significance',
 'affects','association','association not found','confers sensitivity',
 'conflicting data from submitters','conflicting interpretations of pathogenicity',
 'drug response','established risk allele','likely risk allele','protective','risk factor',
 'uncertain risk allele','other','others'}

desc, expl = {}, {}
for line in open('vocab.txt', encoding='latin-1'):
    p = line.rstrip('\n').split('\t')
    if len(p) != 3: continue
    kind, n, v = p[0], int(p[1]), p[2]
    if   kind == 'DESCRIPTION': desc[v] = desc.get(v, 0) + n
    elif kind == 'EXPLANATION': expl[v] = expl.get(v, 0) + n

def cd(s): return [x.strip() for x in re.split(r'[/,;]', s.lower())]
def ce(s):
    out = []
    for part in re.split(r'[/;]', s.lower()):
        i = part.find('(')
        out.append((part if i < 0 else part[:i]).strip())
    return out

bad = {}
for v, n in desc.items():
    for c in cd(v):
        if c and c not in VALID: bad[c] = bad.get(c, 0) + n
for v, n in expl.items():
    for c in ce(v):
        if c and c not in VALID: bad[c] = bad.get(c, 0) + n

for k, n in sorted(bad.items(), key=lambda x: -x[1]):
    print(f'{n:>9,}  {k!r}')
```

Compare the `REVIEWSTATUS` rows against `ReviewStatusNameMapping` by hand — there are under a dozen.

### Check 19 — Recover allow-lists from the deployed binary

.NET string literals are UTF-16LE in the `#US` heap, so plain `strings` misses them.

```sh
python3 - <<'PY'
import re
d = open('build_output/SAUtils.dll', 'rb').read()
out = set()
for m in re.finditer(rb'(?:[\x20-\x7e]\x00){4,}', d):
    out.add(m.group(0).decode('utf-16-le'))
pat = r'^(benign|likely |pathogenic|uncertain|risk |protective|affects|association|drug |' \
      r'conflicting|no assertion|no classification|criteria provided|reviewed by|' \
      r'practice guideline|established|low penetrance|confers)'
for t in sorted(x for x in out if re.match(pat, x, re.I)):
    print(' ', repr(t))
PY
```

This both recovers the allow-lists for checks 17-18 *and* verifies that the patches are present in
the artifact that actually built the `.nsa` — `uncertain risk allele` appearing confirms patch #1.

---

## Appendix B — optional checks (not run)

### D — Vid-universe overlap

Run this **first** in any future comparison. If the variant universes differ materially, "lost
ClinVar" and "this variant wasn't called in that run" get conflated and every later number is
suspect.

```sql
WITH
  n AS (SELECT DISTINCT vid FROM `{{NEW}}`),
  o AS (SELECT DISTINCT vid FROM `{{OLD}}`)
SELECT
  (SELECT COUNT(*) FROM n)                                                        AS vids_new,
  (SELECT COUNT(*) FROM o)                                                        AS vids_old,
  (SELECT COUNT(*) FROM (SELECT vid FROM n INTERSECT DISTINCT SELECT vid FROM o)) AS vids_shared,
  (SELECT COUNT(*) FROM (SELECT vid FROM o EXCEPT DISTINCT SELECT vid FROM n))    AS vids_only_old,
  (SELECT COUNT(*) FROM (SELECT vid FROM n EXCEPT DISTINCT SELECT vid FROM o))    AS vids_only_new
```

### E — Decompose newly 1-star vids

```sql
WITH
  old_stars AS (SELECT DISTINCT vid, s FROM `{{OLD}}`, UNNEST(clinvar_rcv_num_stars) AS s),
  new_stars AS (SELECT DISTINCT vid, s FROM `{{NEW}}`, UNNEST(clinvar_rcv_num_stars) AS s),
  old_vids  AS (SELECT DISTINCT vid FROM `{{OLD}}`),
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

### G — Per-column fill rates

Near-moot for this schema (see §6), but useful if the loader ever stops writing the six ClinVar
columns atomically.

```sql
WITH
  shared AS (
    SELECT vid FROM `{{NEW}}` INTERSECT DISTINCT SELECT vid FROM `{{OLD}}`),
  per_vid AS (
    SELECT 'old' AS tbl, vid,
           LOGICAL_OR(ARRAY_LENGTH(clinvar_classification)      > 0) AS classification,
           LOGICAL_OR(clinvar_last_updated IS NOT NULL)              AS last_updated,
           LOGICAL_OR(ARRAY_LENGTH(clinvar_phenotype)           > 0) AS phenotype,
           LOGICAL_OR(ARRAY_LENGTH(clinvar_rcv_ids)             > 0) AS rcv_ids,
           LOGICAL_OR(ARRAY_LENGTH(clinvar_rcv_classifications) > 0) AS rcv_classifications,
           LOGICAL_OR(ARRAY_LENGTH(clinvar_rcv_num_stars)       > 0) AS rcv_num_stars
    FROM `{{OLD}}` WHERE vid IN (SELECT vid FROM shared) GROUP BY vid
    UNION ALL
    SELECT 'new', vid,
           LOGICAL_OR(ARRAY_LENGTH(clinvar_classification)      > 0),
           LOGICAL_OR(clinvar_last_updated IS NOT NULL),
           LOGICAL_OR(ARRAY_LENGTH(clinvar_phenotype)           > 0),
           LOGICAL_OR(ARRAY_LENGTH(clinvar_rcv_ids)             > 0),
           LOGICAL_OR(ARRAY_LENGTH(clinvar_rcv_classifications) > 0),
           LOGICAL_OR(ARRAY_LENGTH(clinvar_rcv_num_stars)       > 0)
    FROM `{{NEW}}` WHERE vid IN (SELECT vid FROM shared) GROUP BY vid)
SELECT col AS clinvar_column,
       COUNTIF(tbl = 'old' AND populated) AS old_vids,
       COUNTIF(tbl = 'new' AND populated) AS new_vids,
       COUNTIF(tbl = 'new' AND populated) - COUNTIF(tbl = 'old' AND populated) AS delta,
       ROUND(SAFE_DIVIDE(COUNTIF(tbl = 'new' AND populated),
                         COUNTIF(tbl = 'old' AND populated)), 4) AS new_over_old
FROM per_vid
UNPIVOT(populated FOR col IN (classification, last_updated, phenotype,
                              rcv_ids, rcv_classifications, rcv_num_stars))
GROUP BY col
ORDER BY new_over_old
```

### Classification vocabulary shift

Not formally part of the check list, but a cheap overview of how the classification vocabulary moved
between the two tables:

```sql
WITH
  old_c AS (SELECT c AS classification, COUNT(DISTINCT vid) AS old_vids
            FROM `{{OLD}}`, UNNEST(clinvar_classification) AS c GROUP BY c),
  new_c AS (SELECT c AS classification, COUNT(DISTINCT vid) AS new_vids
            FROM `{{NEW}}`, UNNEST(clinvar_classification) AS c GROUP BY c)
SELECT classification, IFNULL(old_vids, 0) AS old_vids, IFNULL(new_vids, 0) AS new_vids,
       IFNULL(new_vids, 0) - IFNULL(old_vids, 0) AS delta
FROM old_c FULL OUTER JOIN new_c USING (classification)
ORDER BY delta
```

---

## Changelog

- **2026-09-15** — Initial version. Checks 1-11 recorded; `20-2403658-G-A` resolved against source
  XML; open questions A-G defined.
- **2026-09-15** — Check 12: live ClinVar cross-check of `VCV000337923`. Confirmed RCV completeness
  and recorded post-snapshot drift. **Corrected an earlier hypothesis:** that variant has never been
  pathogenic, so it does not explain the pathogenic loss.
- **2026-09-15** — Check 13: same-version control group returned all zeros across 1607 RCVs.
- **2026-09-15** — Check 14: RCV-level losses returned 17 dropped RCVs across 16 vids, all 1-star —
  the first non-zero result in the validation.
- **2026-09-15** — Check 15: all 17 resolved against the XML. 15 retired by ClinVar; the 2 pathogenic
  records are the FMR1 microsatellite `c.-128GGC[55_200]`, correctly detached from an adjacent
  deletion. **Zero defects.**
- **2026-09-15** — Check 16: full `clinvar_classification` sweep returned exactly one pathogenic
  loss, the FMR1 case already resolved.
- **2026-09-16** — Checks 17-19 close the patch-behavior question **without re-running the
  converter**: allow-lists recovered from `SAUtils.dll` and diffed against the full 2025-07
  vocabulary. No real classification term is filtered; patch #1 confirmed present and load-bearing.
- **2026-10-01** — Restructured as a reusable model. Added **Appendix A** (every query and command,
  keyed by check), **Appendix B** (the three optional checks), and **§8** (suggested run order,
  principles, schema traps). Renumbered conclusions, folded both resolved anomalies into §4, and
  cleared stale cross-references left by incremental editing.
- **2026-10-02** — Added §9, reproducing the results through the VS-2029 WDL with and without
  reference disks. Pointed patch and binary references at the checked-in diff and its build.
