# Loading VCF Headers in the Parquet Code Path

**Tickets:** VS-1968 — *[Spike] How to load headers in Parquet code path*

**Motivation:** To support GVS on TSPS and the next All of Us callset. Today VCF header
loading works only on the legacy BigQuery Write API ingest path; the Parquet ingest
path (now the default for bulk/AoU) cannot load headers and actively blocks the attempt.

---

## 1. Current state

### 1.1 The header data model (unchanged, works on the BQ path)

Three tables, all created *conditionally* on `load_vcf_headers` in
`GvsAssignIds.wdl:96-129` (schemas at `GvsAssignIds.wdl:30-32`):

| Table | Columns | Role |
|---|---|---|
| `vcf_header_lines_scratch` | `sample_id`, `vcf_header_lines` (NULLABLE), `vcf_header_lines_hash`, `is_expected_unique` | Transient staging during ingest |
| `vcf_header_lines` | `vcf_header_lines_hash`, `vcf_header_lines`, `is_expected_unique` | Deduplicated header text, one row per distinct hash |
| `sample_vcf_header` | `sample_id`, `vcf_header_lines_hash` | Many-to-many map: which chunks each sample has |

A sample's header is split into **chunks** by
`CreateVariantIngestFiles.buildAllVcfLineHeaders()` (`CreateVariantIngestFiles.java:524-537`):

- each `CommandLine` header line individually → `is_expected_unique = true` (varies per sample);
- all remaining header lines joined into one blob → `is_expected_unique = false`
  (near-identical across *all* samples — this blob is the entire dedup cost story).

Each chunk is MD5-hashed; the hash is the dedup key.

### 1.2 Where dedup actually happens today

Dedup is done **entirely at write time** in `VcfHeaderLineScratchCreator.java:93-108`
(BQ path): for each chunk it runs existence checks against both the scratch and the
already-loaded tables (`VcfHeaderLineScratchCreator.java:43-52`, backed by
`BigQueryUtils.doRowsExistFor` at `BigQueryUtils.java:455`). If the hash already exists,
it writes a row with **NULL** text (association only); otherwise it writes the full text.

The scratch→final merge in `process_sample_vcf_headers.py` does **NOT** deduplicate:

```sql
-- process_sample_vcf_headers.py:21  (no DISTINCT / GROUP BY)
INSERT INTO ...vcf_header_lines (vcf_header_lines_hash, vcf_header_lines, is_expected_unique)
SELECT vcf_header_lines_hash, vcf_header_lines, is_expected_unique
FROM ...vcf_header_lines_scratch
WHERE vcf_header_lines IS NOT NULL;   -- relies on write-time NULLing for uniqueness
```

**Consequence for Parquet:** any design that drops the write-time existence checks
*must* add deduplication somewhere else, or `vcf_header_lines` gets N copies of the
shared blob.

### 1.3 What the Parquet path does today

- `VcfHeaderLineScratchCreator.java:109-113` — the `PARQUET` branch is a **stub**: it
  writes only `sample_id` + `headerLineHash`. No header text, no `is_expected_unique`,
  no dedup.
- The Parquet header schema `headersRowSchema` (`CreateVariantIngestFiles.java:226`)
  declares only 2 fields vs the 4-column scratch table.
- `HeaderParquetFileWriter.writeJson` (`HeaderParquetFileWriter.java:38`) mirrors that
  2-field shape.
- The header Parquet file (`header_file.parquet`, `VcfHeaderLineScratchCreator.java:76`)
  is written locally but **never uploaded** — the upload block globs only
  `vet_*`/`ref_*`/`sample_chromosome_ploidy_*` (`GvsImportGenomes.wdl:586-596`) and then
  `rm *.parquet` (`GvsImportGenomes.wdl:599`) deletes it.
- The load chain never sees headers: `DiscoverParquetFiles` prefixes are only
  `sample_chromosome_ploidy` / `vet` / `ref_ranges` (`GvsImportGenomes.wdl:278-279`).
- The combination is explicitly blocked: `GvsImportGenomes.wdl:97-100`
  ("*The combination of Parquet ingest and VCF header loading is not currently supported*").
- A gate placeholder is left at `GvsImportGenomes.wdl:252`
  ("*add a gate for Parquet header loading here once that's implemented*").

### 1.4 TODO in Parquet

The Parquet-generating task already localizes each VCF once
(`GvsImportGenomes.wdl:551-553`) and runs a **single** `CreateVariantIngestFiles` call
(`GvsImportGenomes.wdl:560`) that already passes `--enable-vcf-headers`
(`GvsImportGenomes.wdl:575`). Header extraction piggybacks on the existing VCF read —
**no second localization is required**. The work is finishing the plumbing, not adding a
task.

---

## 2. Key findings / constraints

1. **The merge does not dedup** (§1.2). Naive Parquet writing requires a dedup step to be
   added.
2. **`is_expected_unique` is an unused lever.** It cleanly separates the one shared blob
   (dedup-worthy) from per-sample command lines (not dedup-worthy). Designs can route on it.
3. **Per-sample invocations are isolated.** Each `CreateVariantIngestFiles` run sees one
   gVCF and cannot see other samples' chunks without a shared store — so cross-sample
   write-time dedup (as the BQ path does via BQ queries) is not naturally available offline.
4. **Incremental ingest must be handled.** Callsets add samples in batches; the BQ path
   checks *already-loaded* tables (`VcfHeaderLineScratchCreator.java:43-52`). Any Parquet
   design must avoid re-inserting a hash already present in `vcf_header_lines` from a prior
   batch.
5. **The current BQ dedup is itself expensive** — see §5.

---

## 3. Design options

### Option 1 — Naive write, dedup in the merge SQL (low risk)
- Writer: emit full text + `is_expected_unique` for every chunk, every sample; drop the
  existence queries.
- Merge: change `process_sample_vcf_headers.py:21` from a blind INSERT to a **MERGE** into
  `vcf_header_lines` on `vcf_header_lines_hash` (or `SELECT ... GROUP BY hash`), so it
  dedups *and* is idempotent across batches (satisfies constraint #4).
- **Cost:** max (transient) scratch bloat; one large dedup/MERGE query scanning the text
  column; zero per-sample BQ round-trips during ingest.
- **Pros:** simplest; keeps existing scratch + load machinery; ingest stays fully offline
  and parallel (the point of the Parquet path).
- **Cons:** transient N× storage; large dedup query.

### Option 2 — Content-addressed header text (larger change)
- Split writer output: (a) tiny per-sample **associations** (`sample_id`, `hash`,
  `is_expected_unique`) → `sample_vcf_header`; (b) header **text** written to a GCS path
  keyed by hash, e.g. `.../headers/text/<hash>.parquet`, with an `ifGenerationMatch=0`
  precondition so identical concurrent writes are skipped.
- Load: distinct hash-named objects → `vcf_header_lines` (one object = one row = inherently
  deduped). Incremental ingest is free (object/hash already exists → skip).
- **Pros:** shared blob stored once; no dedup query; **scratch table can be eliminated**.
- **Cons:** `is_expected_unique = true` lines can explode into ~N tiny objects (GCS
  object-count/listing cost); needs a bespoke loader (doesn't fit the superpartition FOFN
  grouping in `DiscoverParquetFiles`).

### Option 3 — Hybrid, routed by `is_expected_unique` (probably optimal for implementation)
- `false` (shared blob): content-address it (Option 2) → stored once. The big win.
- `true` (command lines): write naively inline with associations → no small-object
  explosion, no dedup benefit lost (they were never going to dedup).
- **Pros:** cost-optimal; uses the schema flag as intended.
- **Cons:** two code paths in the writer.

### Option 4 — Pre-load local consolidation
- Dedup header Parquet by hash in a consolidation step before load (git precedent:
  "Experimenting with local consolidation and deduping of parquet data").
- **Cons:** already paid N× GCS writes; only saves BQ-side cost. Weaker than 2/3.

**Recommendation:** implement **Option 1 first** (low risk, unblocks TSPS/AoU, benchmarks
the "naive" cost for AC1), then evaluate **Option 3** if the naive cost is unacceptable.
Options 1 and 3 share the same writer-schema and WDL wiring work; only dedup differs.

---

## 4. Implementation work breakdown (Option 1, with Option 3 notes)

### 4.1 Java — writer + schema
- `CreateVariantIngestFiles.java:226` — expand `headersRowSchema` to the full 4 columns
  (`sample_id`, `vcf_header_lines` text, `vcf_header_lines_hash`, `is_expected_unique`).
- `VcfHeaderLineScratchCreator.java:109-113` — flesh out the `PARQUET` branch to write the
  full record (text + hash + `is_expected_unique`); drop the write-time existence checks on
  this path.
- `HeaderParquetFileWriter.java:38` (`writeJson`) — extend to the full record shape.
- *(Option 3)* route on `is_expected_unique`: content-address the shared blob, inline the
  unique lines.
- Reference template: commits `e852dcf27` / `2d48b4a84` (ref & ploidy Parquet writers) and
  `gvs_ingest_refactoring.md` (note VS-1803 is called out there as the pending header work).

### 4.2 WDL — `GvsImportGenomes.wdl`
- Remove the block at `:97-100` (`CannotLoadHeadersWithParquetIngest`).
- Upload the header Parquet in the generate task: extend the glob/copy block
  (`:586-596`) to pick up `header_file.parquet` and copy it to a new
  `.../vcf_header_lines_scratch/` GCS dir (before the `rm *.parquet` at `:599`).
- `DiscoverParquetFiles` (`:273`, prefixes at `:278-279`): add `vcf_header_lines_scratch`
  as a `regular_table_prefixes` entry so the FOFN grouping (`parse_and_group_files.py`)
  picks it up. The loader derives the table name from the FOFN filename
  (`load_parquet_to_bq.py:46-52`), so a `vcf_header_lines_scratch.fofn` loads into the
  scratch table with no loader change.
- Gate: wire `ProcessVCFHeaders` (`:249`, gate at `:252`) to also depend on
  `LoadParquetFilesToBQ.done` when `use_parquet_ingest` (currently only
  `LoadDataViaBigQueryWriteAPI.done`).

### 4.3 Python — merge
- `process_sample_vcf_headers.py:21` — change the `vcf_header_lines` INSERT to a
  **MERGE**/dedup on `vcf_header_lines_hash` (dedups the naive-loaded scratch rows *and*
  makes it idempotent across incremental batches). The `sample_vcf_header` INSERT
  (`:24`) stays as-is (naturally unique per sample+hash).
- `load_parquet_to_bq.py` — no change needed for Option 1 (raw load into scratch);
  Option 2/3 would need a content-addressed loader.

### 4.4 Table schema
- No schema change to the three tables. `vcf_header_lines_scratch` stays as the load
  target for Option 1. *(Option 2/3 could drop scratch for the Parquet path.)*

---

## 5. AC1 — cost of a naive header load into BQ

"Naive" = Option 1. Measure and compare against the current BQ Write API path.

**Option 1 (Parquet) cost components:**
- transient scratch storage ≈ N samples × (blob + command-line text) bytes;
- the BQ Parquet load-job bytes;
- the dedup/MERGE query bytes scanned (dominated by the large text column).

**Current BQ Write API path baseline (for comparison):**
- `doRowsExistFor` runs a `SELECT COUNT(*)` per lookup (`BigQueryUtils.java:455`), called
  **twice per chunk** (scratch + non-scratch) × ~2 chunks/sample ≈ **~4 BQ query jobs per
  sample** — on the order of **~1.6M query jobs for a 400k-sample callset**, serialized
  inside each ingest task, plus the streaming inserts.

Escaping that per-sample query storm is a primary motivation; the naive Parquet approach
trades it for transient storage + one dedup query.

**Suggested measurement harness:** `GvsQuickstartVcfIntegration.wdl` (smallest end-to-end
path with `load_vcf_headers = true`); capture bytes/slot-ms from the `cost_observability`
table and BQ job stats for both paths at a fixed sample count, then extrapolate.

---

## 6. Comparative cost of the design options

AC1 (§5) covers the naive load (Option 1). This section puts Options 2–4 and the current
BQ path on the same cost model so they can be compared directly. **Bottom line up front:**
at header data volumes every Parquet option costs pennies in BQ/GCS dollars — the real
differentiators are *transient storage footprint*, *GCS object count / small-file
handling*, and *engineering complexity*, and the current BQ path's cost is an *operational*
query storm rather than dollars.

### 6.1 Cost model & parameters

Per ingest batch of **N** samples (AoU ≈ 400k):

- **B** — size of the shared non-command-line blob (`is_expected_unique = false`),
  ~10–100 KB (working figure: 50 KB).
- **k** — command-line chunks per sample (`is_expected_unique = true`), ~1–3, each ~**C**
  bytes (small, ~200 B).
- **d_B** — distinct shared blobs across the batch (usually 1; a few if samples span
  pipeline/reference versions).
- **d_C** — distinct command-line chunks (up to ≈N if genuinely per-sample; fewer if
  pipelines share command lines).

Cost dimensions (the BQ/GCS pricing levers):

1. **Ingest-time BQ queries** — query jobs fired *inside* each `CreateVariantIngestFiles` run.
2. **GCS write ops + bytes** uploaded.
3. **GCS transient storage** — held until the 14-day lifecycle delete.
4. **GCS object count** — drives discovery `ls` cost and BQ small-file load overhead.
5. **BQ load jobs** — free in dollars, but subject to quota + small-file overhead.
6. **BQ query bytes scanned for dedup** — the on-demand dollar driver.
7. **Final BQ storage** — deduped; ~identical and negligible across all options.

Representative unit prices (US, *verify current*): on-demand query **$6.25 / TiB scanned**;
GCS Standard **~$0.02 / GB-month**; BQ Storage Write API **~$0.025 / GB** (streaming); BQ
batch **load jobs are free**.

### 6.2 Comparison table (dominant term per cell)

| Dimension | Current BQ Write API | Opt 1 Naive+MERGE | Opt 2 Content-addr (all) | Opt 3 Hybrid | Opt 4 Pre-load consolidate |
|---|---|---|---|---|---|
| Ingest BQ queries | **~2(k+1)·N ≈ 1.6M** | 0 | 0 | 0 | 0 |
| GCS bytes written | 0 | **N·B** | d_B·B + d_C·C | N·k·C + d_B·B | **N·B** |
| GCS transient storage | 0 | **N·B** | d_B·B + d_C·C | N·k·C + d_B·B | N·B (then small) |
| GCS object count | 0 | N | **d_B + d_C (→ ~N tiny)** | N (small) + d_B | N → small |
| BQ load jobs | 0 (stream) | N/10k | ~N/10k (tiny files) | N/10k (small) | few |
| BQ dedup scan | 0 | **N·B** | 0 | N·k·C | 0 (local) |
| Streaming $ | d_B·B + assoc | 0 | 0 | 0 | 0 |
| Extra compute | — | — | — | — | **consolidation VM reads N·B** |
| Complexity / risk | (baseline) | Low | High | Medium | Med-High |

Bold marks the dominant cost for each option. Note `C ≪ B`, so any `N·k·C` term is orders
of magnitude below an `N·B` term.

### 6.3 Per-option read

- **Current BQ path** — negligible dollars (small bytes, cheap streaming) but its cost is
  *operational*: ~1.6M query jobs against per-project query rate limits / concurrency,
  serialized inside ingest tasks (`VcfHeaderLineScratchCreator.java:93-108`,
  `BigQueryUtils.java:455`). This is the fragility the Parquet path exists to escape, not a
  dollar problem.
- **Option 1 (naive)** — the only nontrivial costs are transient GCS storage of **N·B**
  (every sample's copy of the blob) for up to 14 days, plus one dedup MERGE that scans
  **N·B** once. Both are ~N·B — cheap in absolute dollars (§6.4) but the largest
  storage/scan footprint of the Parquet options. Zero ingest queries; lowest complexity.
- **Option 2 (content-addressed, all chunks)** — eliminates N·B: the blob is stored once,
  dedup is free via the object namespace, no MERGE scan. But `is_expected_unique = true`
  lines can explode to **~N tiny objects**, which hits (a) `ls` discovery cost/time, (b) BQ
  small-file load overhead, and (c) GCS Class-A op counts from N precondition writes.
  Lowest bytes, worst object-count profile, highest complexity (bespoke loader).
- **Option 3 (hybrid)** — keeps the blob content-addressed (stored once, no MERGE scan for
  it) while writing the small command-line/association data the naive way. Object count ≈ N
  (same as naive) but each object is tiny (no B), so transient storage and any residual
  dedup scan drop from **N·B** to **N·k·C** — orders of magnitude smaller. Best
  cost/footprint overall; moderate complexity (two writer paths + one bespoke blob loader).
- **Option 4 (pre-load consolidation)** — still pays the full **N·B** GCS write, then adds
  a consolidation VM that reads **N·B** to dedup locally before a tiny load. Saves the BQ
  dedup scan but not the write, and adds compute. Dominated by Option 3 on every axis.

### 6.4 Worked example (plausible numbers)

N = 400,000; B = 50 KB; k = 1; C = 200 B; d_B = 1; d_C ≈ 400,000.

| Quantity | Opt 1 | Opt 2 | Opt 3 |
|---|---|---|---|
| GCS bytes written | ~20 GB | ~80 MB | ~80 MB + 50 KB |
| Transient storage (14 d) | ~20 GB → ~$0.14 | ~80 MB → <$0.01 | ~80 MB → <$0.01 |
| GCS objects | ~400k | ~400k tiny | ~400k tiny + 1 |
| BQ dedup scan | ~20 GB → ~$0.13 | 0 | ~80 MB → <$0.01 |

**Takeaway:** every Parquet option is dollar-trivial. The decision is *operational*, not
financial: Option 1's downside is a ~20 GB transient footprint (already auto-cleaned by the
14-day lifecycle rule) and one ~20 GB scan; Options 2/3 shrink that to ~80 MB at the price
of complexity (Opt 2 also trades it for an ~N-tiny-object problem). This argues for
**Option 1 first** (cheapest to build, only real cost is transient storage the lifecycle
already handles) and **Option 3** only if the transient footprint or object count becomes a
problem at full AoU scale.

### 6.5 How to measure each (extends the AC1 harness)

Reuse the `GvsQuickstartVcfIntegration.wdl` harness at a fixed sample count and extrapolate
by N:

- **Ingest queries:** count BQ jobs per ingest step (`cost_observability` / INFORMATION_SCHEMA.JOBS);
  confirm 0 for all Parquet options vs ~2(k+1)·N for the BQ path.
- **GCS bytes/objects:** `gcloud storage du` and object counts on the parquet output dir
  *before* the lifecycle delete.
- **Dedup scan:** `totalBytesProcessed` / `totalSlotMs` from the MERGE job.
- **Load:** job count + per-job duration from `load_parquet_to_bq.py` `stats.json`.
- **Empirically measure the drivers** from one real gVCF header: **B** and **k** (they
  dominate every formula), and **d_B** via
  `SELECT COUNT(DISTINCT vcf_header_lines_hash) ... WHERE is_expected_unique = false`.

---

## 7. AC2 — consistency queries

Consistency queries to be supplied by Aaron Hatcher. They will presumably assert
equivalence between a BQ-loaded and a Parquet-loaded dataset. Expected checks:
- row counts and distinct-hash counts in `vcf_header_lines` match across paths;
- `sample_vcf_header` (sample_id, hash) set is identical across paths;
- every `sample_vcf_header.vcf_header_lines_hash` resolves to a row in `vcf_header_lines`
  (referential integrity; no orphan hashes);
- reconstructed per-sample header text is byte-identical across paths.

*(Placeholder — fill in once queries are provided.)*

---

## 8. Open questions / risks

- **Small-object cost (Option 2/3):** how many distinct `is_expected_unique = true` chunks
  arise at AoU scale? Determines whether content-addressing the unique lines is viable.
- **MERGE cost at scale (Option 1):** the dedup MERGE scans the text column; confirm this
  is cheaper than the query storm it replaces.
- **`billing_project_id` not threaded to the loader** (`GvsImportGenomes.wdl:1248`,
  VS-1955) — verify header loads don't need it.
- **Scratch retention:** Option 1 keeps scratch; confirm `clean_up_scratch_table`
  (`process_sample_vcf_headers.py:28-31`) still applies cleanly after a MERGE.

---

## 9. File / line reference index

| Area | Location |
|---|---|
| Parquet header writer stub | `VcfHeaderLineScratchCreator.java:109-113` |
| BQ write-time dedup | `VcfHeaderLineScratchCreator.java:93-108` |
| Already-loaded existence checks | `VcfHeaderLineScratchCreator.java:43-52`, `BigQueryUtils.java:455` |
| Header chunking | `CreateVariantIngestFiles.java:524-537` |
| Parquet header schema (2-field) | `CreateVariantIngestFiles.java:226` |
| Header Parquet file writer | `HeaderParquetFileWriter.java:38` |
| Table schemas / creation | `GvsAssignIds.wdl:30-32`, `GvsAssignIds.wdl:96-129` |
| Scratch→final merge (no dedup) | `process_sample_vcf_headers.py:21,24,28-31` |
| Blocking guard | `GvsImportGenomes.wdl:97-100` |
| Single-invocation ingest + `--enable-vcf-headers` | `GvsImportGenomes.wdl:551-575` |
| Parquet upload glob | `GvsImportGenomes.wdl:586-599` |
| Discover / load / gate | `GvsImportGenomes.wdl:278-279`, `:288`, `:252` |
| Loader (table name from FOFN) | `load_parquet_to_bq.py:16,46-52` |
