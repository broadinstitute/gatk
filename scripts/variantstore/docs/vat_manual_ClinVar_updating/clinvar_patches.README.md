# `clinvar_patches.diff` — what it is and how to use it

Self-contained instructions for rebuilding Nirvana's `SAUtils` with ClinVar support for
ClinVar releases from 2024 through 2025-07, the last release NCBI published in the XML
format SAUtils reads (see VS-2038). Assume no other context.

## Why this patch exists

Nirvana's ClinVar parsers hard-code allow-lists for clinical significance and review status and
**throw** on anything unrecognized. Since the germline/somatic classification split, NCBI emits
placeholder text — `no classifications from unflagged records` — in both fields for variants whose
submissions have all been flagged. Stock Nirvana aborts about 15 seconds into the run:

```
ERROR: Invalid clinical significance found. Observed: no classifications from unflagged records
  at SAUtils.InputFileParsers.ClinVar.ClinVarVariationReader.GetSignificances(...)
```

Nirvana's ClinVar code has had no upstream commits since June 2022, so this is an unfixed upstream
bug, not a misconfiguration. Without this patch you cannot convert any recent ClinVar release.

## What the patch changes

| # | File                        | Change                                                                                                                                                       |
|---|-----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1 | `ClinVarCommon.cs`          | Adds `uncertain risk allele` to `ValidPathogenicity`. **Not cosmetic** — 623 occurrences in the 2025-07 release would otherwise be silently discarded.       |
| 2 | `ClinVarVariationReader.cs` | `ResolveReviewStatus()`: safe lookup falling back to `no_assertion` instead of a raw indexer that threw `KeyNotFoundException`.                              |
| 3 | `ClinVarCommon.cs`          | `GetSignificances()` filters tokens against `ValidPathogenicity`, warning and skipping instead of throwing. One function, covers both the RCV and VCV paths. |
| 4 | `NsaWriter.cs`              | **Diagnostic progress logging only.** Added while chasing an apparent hang (see step 5). No functional change — safe to drop if you want a minimal patch.    |

---

## 1. Clone the exact tag

```sh
git clone --branch v3.18.1 --depth 1 https://github.com/Illumina/Nirvana.git
cd Nirvana
```

**`v3.18.1` is not arbitrary.** It is the exact version inside the
`us.gcr.io/broad-dsde-methods/variantstore:nirvana_2022_10_19` image and on the Terra Nirvana
reference disk. Building from any other tag risks a `.nsa` schema version that the production
pipeline will not read. Verify before continuing if you are unsure.

## 2. Apply

```sh
git apply /path/to/clinvar_patches.diff      # or: patch -p1 < clinvar_patches.diff
```

Paths in the diff are relative to the repository root.

## 3. Build

The project targets `net6.0`. If you have no local .NET 6 SDK, build in a container:

```sh
docker run --rm --platform linux/amd64 \
  -v "$PWD":/src -w /src \
  mcr.microsoft.com/dotnet/sdk:6.0 \
  dotnet publish SAUtils/SAUtils.csproj -c Release -o /src/build_output
```

> ⚠️ **`--platform linux/amd64` is required on Apple Silicon.** The repo ships
> `libBlockCompression.so` as x86_64-only. Without the flag you get the arm64 SDK image and a
> misleading failure: `ERROR: Unable to find the block GZip compression library
> (BlockCompression)`. That is an architecture problem, not a code problem.

Output includes `SAUtils.dll`, every dependency DLL, and the native `libBlockCompression.so` —
nothing needs copying in by hand. The DLLs are architecture-agnostic MSIL, so you can build here
and run elsewhere.

## 4. Gather inputs

| File                                      | Where from                                           |
|-------------------------------------------|------------------------------------------------------|
| `ClinVarFullRelease_<ver>.xml.gz`         | NCBI FTP — the RCV release (~5 GB)                   |
| `ClinVarVariationRelease_<ver>.xml.gz`    | NCBI FTP — the VCV release (~4.8 GB)                 |
| `Homo_sapiens.GRCh38.Nirvana.dat`         | Nirvana reference sequence, from the existing bundle |
| `ClinVarFullRelease_<ver>.xml.gz.version` | **You must create this by hand** — see below         |

NCBI does not ship the `.version` sidecar; `SAUtils` requires it and derives the output filenames
from it. Create it next to the RCV XML:

```
NAME=ClinVar
VERSION=2025-07
DATE=2025-07-15
DESCRIPTION=
```

`NAME` and `VERSION` determine the output name (`ClinVar_2025-07.nsa`). `DATE` is cosmetic.

## 5. Run — on native x86_64, not emulated

```sh
dotnet build_output/SAUtils.dll clinvar \
  --ref Homo_sapiens.GRCh38.Nirvana.dat \
  --rcv ClinVarFullRelease_2025-07.xml.gz \
  --vcv ClinVarVariationRelease_2025-07.xml.gz \
  --out clinvar_db_output
```

> 🛑 **Do not run the conversion under QEMU emulation.** On an Apple Silicon Mac the parse and merge
> phases complete normally, then the **write phase appears to hang**: CPU pegged near 100%, memory
> climbing, and *zero* output for over an hour. It is not deadlocked and it is not a bug in the
> patch — it is an emulation artifact. Patch #4's progress logging exists solely because this was so
> hard to diagnose.
>
> On native x86_64 the same write phase finishes **all** chromosomes in about 12.5 minutes.

Reference environment that worked: GCP VM, `e2-standard-8`, Ubuntu 22.04, 50 GB SSD, .NET 6
**runtime** only (no SDK, no Docker). Install via Microsoft's `dotnet-install.sh --channel 6.0
--runtime dotnet`. Budget **~50 minutes wall clock and ~10 GB peak RAM**; the two XMLs are ~10 GB,
so size the disk accordingly.

## 6. Verify it worked

Expected stderr warnings — these confirm the patches are firing, they are **not** errors:

```
WARNING: skipping unrecognized clinical significance value: 'no classifications from unflagged records'
WARNING: unrecognized review status '...', defaulting to no_assertion
WARNING: skipping unrecognized clinical significance value: 'low penetrance'
```

Expected parse stats for 2025-07: `rcvCount: 4,667,952`, `vcvCount: 3,511,919`, `0 unknown VCVs`.

Expected output — exactly three files:

```
ClinVar_2025-07.nsa          114,567,339 bytes
ClinVar_2025-07.nsa.idx            4,342 bytes
ClinVar_2025-07.nsa.schema           902 bytes
```

**No `.nsi` is produced, and that is correct.** `SAUtils` has no ClinVar interval writer — Illumina
builds the `.nsi` with unreleased internal tooling. It does not matter for the VAT: `.nsa` and
`.nsi` readers are duplicate-key-checked in separate namespaces, and the VAT loader reads per-allele
`.nsa` data only. A stale `.nsi` on the reference disk coexists with a new `.nsa` without conflict.

## 7. Before you deliver it downstream

Nirvana merges every `.nsa` from every `--sd` directory into one list and rejects two readers
sharing a JSON key. ClinVar's key is the hardcoded constant `clinvar`, baked in at creation time
regardless of filename or version.

**You cannot simply add a second `--sd` directory containing the new ClinVar** — it collides with
the bundled one and aborts with `Duplicate variant-level JSON keys found for: clinvar`. Nor can the
old `.nsa` be deleted in place, because the Terra reference disk is mounted read-only.

`GvsCreateVATfromVDS.wdl` handles this when `use_manual_clinvar_update` is true. `AnnotateVCF`
builds a writable directory of symlinks to every bundled supplementary annotation file except
`ClinVar_*.nsa*`, links in the three manually built ClinVar files, and points Nirvana's `--sd`
at that directory. This works the same with or without reference disks. A guard fails the task
unless exactly one ClinVar `.nsa` ends up in the directory. To deliver a new build, upload its three
files to GCS and set `manual_clinvar_path_prefix` to their common path prefix (the path without the
`.nsa` extension). The default is `gs://gvs_quickstart_storage/Nirvana/ClinVar/ClinVar_2025-07`.

## Known behavior worth knowing

- The patches discard only values that are **not** clinical classifications — `not provided`,
  `flagged submission`, and NCBI's placeholder. A vocabulary audit of the full 2025-07 release
  confirmed no real significance term is filtered.
- One edge case: a record whose `<Explanation>` reads `Pathogenic, low penetrance` is dropped
  entirely, because that code path splits on `/` and `;` but not `,`. Exactly 1 record in 10.18
  million, and a pre-existing inconsistency in Nirvana's `GetSignificances` rather than something
  these patches introduce.

## Further reading

These sit alongside this file in `scripts/variantstore/docs/vat_manual_ClinVar_updating/`:

- [`summary.md`](summary.md) — two-page overview of the build and its validation
- [`validation.md`](validation.md) — the evidence that the resulting VAT differs from its
  predecessor only because ClinVar itself changed
