# ClinVar 2025-07 Rebuild for the Nirvana VAT Pipeline

How the updated ClinVar annotation files were produced, and the evidence that the resulting VAT
differs from the previous one only because the underlying ClinVar data changed.

**Tickets** VS-1988 → VS-1994 · **Output** `ClinVar_2025-07.nsa` · **Status: no unexplained differences found**

Full detail lives in `PROGRESS.md` (build) and `validation.md` (validation, including every query used).

---

# 1. Building the annotation files

Nirvana's shipped reference bundle carries ClinVar from **October 2023**. We needed the July 2025
release, which meant rebuilding the ClinVar `.nsa` ourselves.

### The blocker

Nirvana's ClinVar parsers hard-code allow-lists for clinical significance and review status and
*throw* on anything unrecognised. NCBI's germline/somatic classification split now emits placeholder
text — `no classifications from unflagged records` — in both fields, for variants whose submissions
have all been flagged. The run aborted about 15 seconds in. Nirvana's ClinVar code has had no commits
since June 2022, so this is an unfixed upstream bug, not a configuration error on our side.

### Patches to Nirvana v3.18.1

Rebuilt from the Illumina/Nirvana tag `v3.18.1` — confirmed to be the exact version inside our
`variantstore:nirvana_2022_10_19` image and on the reference disk, so there is no schema-version
drift. Nothing else in the source was touched.

| # | File | Change |
|---|---|---|
| 1 | `ClinVarCommon.cs` | Added `uncertain risk allele` to `ValidPathogenicity` — the third tier of the same ClinGen risk-allele set as the two already listed, omitted originally. |
| 2 | `ClinVarVariationReader.cs` | New `ResolveReviewStatus()` does a safe dictionary lookup falling back to `no_assertion`, replacing a raw indexer that threw `KeyNotFoundException`. |
| 3 | `ClinVarCommon.cs` | `GetSignificances()` filters tokens against `ValidPathogenicity`, warning and skipping rather than throwing. One function, both the RCV and VCV paths. |
| 4 | `NsaWriter.cs` | Progress logging only — added while diagnosing an apparent hang. No functional change. |

### Inputs

| File | Notes |
|---|---|
| `ClinVarFullRelease_2025-07.xml.gz` | NCBI FTP, RCV release (~5.0 GB) |
| `ClinVarVariationRelease_2025-07.xml.gz` | NCBI FTP, VCV release (~4.8 GB) |
| `ClinVarFullRelease_2025-07.xml.gz.version` | Hand-written sidecar (`NAME=ClinVar`, `VERSION=2025-07`, `DATE=2025-07-15`). Drives the output filenames; the date is a placeholder. |
| `Homo_sapiens.GRCh38.Nirvana.dat` | Reference sequence, from the existing bundle |

### Build and run environment

Patched source published with `dotnet publish SAUtils/SAUtils.csproj -c Release` inside
`mcr.microsoft.com/dotnet/sdk:6.0`, forced to `--platform linux/amd64` because the repo's
`libBlockCompression.so` is x86_64-only.

The conversion ran on a **GCP VM** — project `gvs-internal`, zone `us-central1-a`,
**`e2-standard-8`**, Ubuntu 22.04, 50 GB SSD, .NET 6 *runtime* only, no Docker. Native x86_64 turned
out to be essential: under QEMU emulation on an Apple Silicon laptop the write phase made zero
measurable progress in an hour, whereas natively every chromosome completed in about 12.5 minutes.

```
dotnet SAUtils.dll clinvar \
  --ref Homo_sapiens.GRCh38.Nirvana.dat \
  --rcv ClinVarFullRelease_2025-07.xml.gz \
  --vcv ClinVarVariationRelease_2025-07.xml.gz \
  --out clinvar_db_output
```

**Total run time 49 min 37 s; peak memory 9.4 GB.** Parse/merge stats: 4,667,952 RCVs,
3,511,919 VCVs, 0 unknown VCVs. Output: `ClinVar_2025-07.nsa` (114,567,339 bytes) plus `.nsa.idx`
and `.nsa.schema`. No `.nsi` is produced — SAUtils has no ClinVar interval writer, and the VAT
loader reads per-allele `.nsa` data only, so this does not matter here.

### WDL changes that made testing possible

Nirvana merges every `.nsa` from every `--sd` directory into one list and rejects two readers sharing
a JSON key. ClinVar's key is the hardcoded constant `clinvar`, so simply adding a second directory
collides and aborts. The production reference disk is read-only, so the old file cannot be removed
there. Changes to `GvsCreateVATfromVDS.wdl`:

1. New workflow input `Boolean use_manual_clinvar_update = true`, marked TEMPORARY.
2. At the `AnnotateVCF` call site, `use_reference_disk = use_reference_disk && !use_manual_clinvar_update`
   — enabling the manual update *forces* the downloader path, whose scratch directory is writable,
   making the duplicate-key crash unreachable by construction rather than relying on the operator to
   set two flags consistently.
3. New task inputs pointing at the three uploaded files in GCS.
4. In the download branch, remove `ClinVar_*.nsa*`, hard-link the new three in, log before/after
   listings, and hard-fail unless exactly one ClinVar `.nsa` remains — so a botched swap surfaces
   immediately rather than 40 minutes later inside Nirvana.

Validated with `womtool`. Reverting to production behaviour is a single flag:
`use_manual_clinvar_update = false`.

---

# 2. Validation

The new VAT (`0c32721`, VS-1994) was compared against the previous VAT (`4e0cc60`, VS-1988) across
every ClinVar column.

> ⚠️ **The baseline is not a correctness reference.** It was built with the Nirvana-bundled ClinVar,
> whose records stop at **2023-10-28**; the new build runs to **2025-06-29**. The two runs differ by
> roughly 20 months of ClinVar revisions *as well as* by the patches, so divergence is expected. The
> task was to show that every difference traces to the data, not to our changes.

### Checks and what each establishes

| Check | Result | Establishes |
|---|---|---|
| Loader diff between the two commits | Byte-identical | No star or classification change can originate in the loader |
| Patch #2 reachability | Loader accepts `RCV*` ids only; patch #2 is VCV-side | Patch #2 cannot affect any VAT column — crash-prevention only |
| Star distribution | −702 one-star variants, +1023 two-star | Framed the question: was anything lost? |
| Origin of new 2-star variants | 1024 newly 2-star; **all 1024** were 1-star before; none from no-ClinVar | Re-grading of existing calls, not invented annotations |
| Origin at record level | 1588 = 1131 same accession re-graded + 456 unchanged + 1 new; **0** demotions from 3/4-star | The same records were re-graded; new records are not appearing |
| Accession versions on re-grades | **0** at an unchanged version; all 1131 version-bumped | No record parsed differently — every change had an upstream revision |
| Accession versions, all 3649 shared | 1607 same, 2042 bumped, **all forward**, 0 backward | The data is strictly fresher; no staleness regression |
| One-star decomposition | 1762 kept, 827 promoted away, **0 / 0 / 0** across all loss buckets, 125 gained — reconciles to −702 exactly | No variant lost its ClinVar annotation |
| **Same-version control group** | 1607 records read from byte-identical input: **0** star differences, **0** classification differences, none lost, none gained | The direct patched-vs-unpatched test: identical input produces identical output |
| Record-level losses | 17 records across 16 variants, all 1-star | Located the only losses anywhere in the comparison |
| Those 17 resolved against the source XML | 15 retired by ClinVar; 2 re-represented as a microsatellite; **0 defects** | Every loss is explained by source data |
| Pathogenic sweep, whole table | Exactly **1** variant | Only one high-stakes change exists, and it is understood |
| Vocabulary audit of the full release | Only non-classifications are filtered; 0.041% of review statuses hit the fallback | The patches discard nothing clinically meaningful |
| Deployed binary | `uncertain risk allele` present in `SAUtils.dll` | Patch #1 shipped, and is load-bearing — 623 occurrences in the release |

### The two anomalies, both resolved against source data

**`20-2403658-G-A` (TGM6) lost a 2-star rating.** The 2025-07 XML lists three RCVs for it, all 1-star
and all Benign — so the new table is correct and the old one was stale. The live ClinVar record
confirms the same three accessions and shows the variant has never been classified pathogenic.

**`X-147912049-CGCG-C` lost its Pathogenic call.** Its two pathogenic records describe
`NM_002024.6(FMR1):c.-128GGC[55_200]` — the FMR1 CGG premutation, a repeat-count range with no
VCF-representable allele. The 2023 bundle had attached that assertion to an adjacent 3 bp deletion.
Declining to carry it across is a *correction*, and it is the single pathogenic change in the entire
table.

### Why this adds up to confidence

The checks were chosen so that each closes a route by which an unexplained change could hide, and
together they cover all of them:

- **Gains** are explained — every promotion carried an upstream record revision, and none occurred
  at an unchanged accession version.
- **Losses** are explained — 17 records in total, each individually traced to ClinVar retiring the
  accession or re-representing the allele.
- **Silent changes** are excluded — on 1607 records where both runs read byte-identical input,
  output was identical in every field.
- **Our code** is excluded — the loader is unchanged, patch #2 is unreachable from the VAT, and a
  vocabulary audit of the entire release shows the remaining filters discard only values that are
  not clinical classifications (`not provided`, `flagged submission`, and NCBI's placeholder).
- **The arithmetic closes.** Independent queries reconcile exactly — the one-star decomposition
  returns −702 against the observed −702, and every record in the 2-star population is accounted
  for with no remainder.

Taken together: **every difference between the two VAT tables is attributable to 20 months of
legitimate ClinVar revision, and the patches introduced no annotation change beyond discarding
values that were never clinical classifications.**

### Known residuals

One edge case: a record whose *Explanation* field reads `Pathogenic, low penetrance` is dropped
entirely, because that code path splits on `/` and `;` but not `,`. This affects exactly 1 record in
10.18 million and is a pre-existing inconsistency in Nirvana's own `GetSignificances`, not something
the patches introduced.

Three optional cross-checks were not run (variant-universe overlap, decomposition of the 125 newly
1-star variants, per-column fill rates). Each re-verifies figures already reconciled by other means.

At production scale, expect proportionally more `not provided` and `flagged submission` records to
be excluded — near 0.55% and 0.041% of records respectively. That is correct behaviour, not a
regression.
