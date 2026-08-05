# Loading VCF Headers in the Parquet Code Path

**Tickets:**
- **VS-1968** — *[Spike] How to load headers in Parquet code path* — deliverables are this
  design doc and the consistency checks (§7).
- **VS-1803** — *Header loading not supported on Parquet branch* — the **implementation** (the
  Java/Python/WDL changes and the `is_loaded` view fix described here belong to this ticket;
  PR 9403).

**Motivation:** To support GVS on TSPS and the next All of Us callset. Today VCF header
loading works only on the legacy BigQuery Write API ingest path; the Parquet ingest
path (now the default for bulk/AoU) cannot load headers and actively blocks the attempt.

**Two strategies, one implemented.** The design settled on two strategies, named to match the
`HeaderParquetStrategy` enum in the implementation (`VcfHeaderLineScratchCreator.java`):

- **Naive** — write full header text for every chunk of every sample; deduplicate downstream in
  the scratch→final promotion. **This is what shipped in PR 9403** and the only implemented
  strategy.
- **Hybrid** — route on `is_expected_unique`: write per-sample command-line chunks inline, but
  content-address the shared blob to GCS so it is stored once. Cheaper and structurally cleaner on
  concurrency at scale, but the content-addressing half is **not yet implemented** (the enum value
  exists and throws). Documented here as the intended future optimization.

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

The scratch→final merge in `process_sample_vcf_headers.py` did **NOT** deduplicate (original form):

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

### 1.3 What the Parquet path did before this work

- `VcfHeaderLineScratchCreator.java:109-113` — the `PARQUET` branch was a **stub**: it
  wrote only `sample_id` + `headerLineHash`. No header text, no `is_expected_unique`,
  no dedup.
- The Parquet header schema `headersRowSchema` (`CreateVariantIngestFiles.java:226`)
  declared only 2 fields vs the 4-column scratch table.
- `HeaderParquetFileWriter.writeJson` (`HeaderParquetFileWriter.java:38`) mirrored that
  2-field shape.
- The header Parquet file was written locally but **never uploaded** — the upload block globbed
  only `vet_*`/`ref_*`/`sample_chromosome_ploidy_*` and then `rm *.parquet` deleted it.
- The load chain never saw headers: `DiscoverParquetFiles` prefixes were only
  `sample_chromosome_ploidy` / `vet` / `ref_ranges`.
- The combination was explicitly blocked: `GvsImportGenomes.wdl:97-100`
  ("*The combination of Parquet ingest and VCF header loading is not currently supported*").

### 1.4 Header extraction reuses the existing single VCF read

`CreateVariantIngestFiles` already localizes each VCF once (`GvsImportGenomes.wdl:551-553`)
and runs a **single** call (`GvsImportGenomes.wdl:560`) that already accepts
`--enable-vcf-headers` (`GvsImportGenomes.wdl:575`), so header extraction piggybacks on the
existing VCF read with no extra ingest tool — the remaining work was finishing the plumbing
(upload → discover → load), not adding an extraction tool. Note that the operational model
(§1.5) runs header loading as a **separate ingest** from vet/ref, so in practice each phase
localizes its VCFs independently (twice overall) rather than folding everything into one pass.

### 1.5 Operational model: phased header-then-data loading (AoU style)

**Purpose of the header phase.** The header-only ingest is an early **sanity check of the input
gVCFs themselves**. Loading just the headers is one way to inspect each input's metadata
(reference version, contig list, sample names, expected INFO/FORMAT/FILTER definitions, pipeline
provenance, DRAGEN versions, whether reblocking was done, etc.) and confirm the cohort is
well-formed and mutually consistent *before* committing to the expensive vet/ref ingest. This
purpose must be preserved in the Parquet path.

**The supported model.** Header loading runs as its **own** ingest first
(`load_vcf_headers = true`, `load_vet_and_ref_ranges = false`); the loaded headers are then
**manually inspected** with BQ queries to sanity-check the input gVCFs; and **only if that passes**
is a second ingest run started for vet / ref_ranges / sample_chromosome_ploidy
(`load_vcf_headers = false`, `load_vet_and_ref_ranges = true`). Three phases: *load headers →
manually inspect inputs → load everything else*. This is the same operational model AoU already
uses on the BQ path, and it is the **only** model supported on the Parquet path. It localizes the
VCFs **twice** (once per ingest run) — an accepted cost, because the header phase exists precisely
to catch bad inputs before paying for the expensive vet/ref ingest.

The "avoid re-localizing the VCFs in a second task" idea (from VS-1803) is therefore **not**
pursued: the phased model re-localizes by design. Validating inputs before doing expensive
variant work fundamentally requires either re-localizing or doing that expensive work first,
and this model chooses to re-localize so that **no** vet/ref compute is spent until the manual
gate passes.

**Implication:**
- The header load must complete and be **manually validated** before the vet/ref/ploidy
  ingest is started. Because these are separate workflow runs, the gate is enforced naturally
  by the operator between runs (see the WDL wiring in §4.2).
- The sanity check reads header **content**, so the Parquet header load must persist the
  **full header text**, not just hashes. Storing it *deduplicated* does not hurt inspection:
  the check reconstructs any sample's header by joining `sample_vcf_header` (sample → hash) to
  `vcf_header_lines` (hash → text). Both Naive and Hybrid (§3) end with those same two tables, so
  both keep the header lines (contigs, reference, INFO/FORMAT) fully inspectable — the only way to
  lose that is an *incomplete* implementation that persists hashes without ever loading the text
  into `vcf_header_lines` (see the Hybrid scaffold caveat in §3).

---

## 2. Key findings / constraints

1. **The merge does not dedup** (§1.2). Naive Parquet writing requires a dedup step to be
   added.
2. **`is_expected_unique` distinguishes the shared blob from per-sample command lines** —
   `false` for the joined non-command-line blob (shared, dedup-worthy), `true` for the
   command-line chunks (expected to differ per sample). This flag is **already used** by the
   AoU header sanity check: its query filters `is_expected_unique = TRUE` to isolate the
   per-sample command lines and extract the DRAGEN version for validation
   (`AOU_DELIVERABLES.md`). Hybrid (§3) *additionally* routes storage on it.
3. **Per-sample invocations are isolated.** Each `CreateVariantIngestFiles` run sees one
   gVCF and cannot see other samples' chunks without a shared store — so cross-sample
   write-time dedup (as the BQ path does via BQ queries) is not naturally available offline.
4. **Incremental ingest must be handled.** Callsets add samples in batches; the BQ path
   checks *already-loaded* tables (`VcfHeaderLineScratchCreator.java:43-52`). Any Parquet
   design must avoid re-inserting a hash already present in `vcf_header_lines` from a prior
   batch.
5. **The current BQ dedup is itself expensive** — see §5.
6. **Header loading must be idempotent** (explicit design goal — and largely already
   realized). Re-running the load — WDL task retry (`LoadParquetFilesToBQ` runs
   `preemptible: 5`, `maxRetries: 3`), Cromwell resume, or re-ingest of a batch — must
   converge to the same final tables, with no duplicate rows. The **BQ path already
   implements this deliberately** via the per-sample `HEADERS_LOADED` status in
   `sample_load_status` (§6.1); it is not accidental. The gap the Parquet work must close is
   narrower: the Parquet path bypasses that per-sample status check and the promotion step's
   blind INSERT does not dedup — both addressed by the anti-join INSERT (§6).
7. **Loading is phased and gated — as an input-gVCF sanity check** (§1.5). The header-only
   phase exists to validate the *input gVCFs* (reference, contigs, samples, expected fields)
   before the expensive vet/ref ingest, so it must remain a distinct, gated phase in the
   Parquet path too — run as a **separate ingest** from the data load. Consequences: (a) the
   header load must stand alone and be manually validated before the data load is started;
   (b) it must persist **full header text** so the check can inspect content (a hard
   requirement on Hybrid's content-addressing); (c) this separate-run pattern makes
   constraint #6 (idempotency) more pressing, not less.

---

## 3. Design: Naive (built) and Hybrid (future)

### 3.1 Naive — write full text, dedup in the load SQL (implemented, PR 9403)

- **Writer:** emit full text + `is_expected_unique` for every chunk, every sample; drop the
  existence queries. (`HeaderParquetStrategy.NAIVE`.)
- **Load:** change the `process_sample_vcf_headers.py` `vcf_header_lines` INSERT from a blind
  INSERT to an **anti-join INSERT**
  (`INSERT ... SELECT ... LEFT JOIN <target> USING(key) WHERE <target>.key IS NULL`, with
  `GROUP BY`/`DISTINCT` on the source), so it dedups *and* is idempotent across batches
  (satisfies constraint #4). Anti-join INSERT rather than `MERGE` to match GVS convention (GVS
  uses `MERGE` nowhere) and avoid mutating-DML concerns.
- **Cost:** max (transient) scratch bloat; one large dedup query scanning the text column;
  zero per-sample BQ round-trips during ingest.
- **Pros:** simplest; keeps existing scratch + load machinery; ingest stays fully offline and
  parallel (the point of the Parquet path).
- **Cons:** transient N× storage; large dedup query; the anti-join is only *sequentially*
  idempotent (safe in the supported single-promotion wiring — see §6.6).

### 3.2 Hybrid — route by `is_expected_unique` (designed; not yet implemented)

Hybrid keeps the small per-sample data naive but content-addresses the one expensive shared blob:

- **`false` (shared blob):** write the header **text** to a GCS path keyed by its hash, e.g.
  `.../headers/text/<hash>.parquet`, with an `ifGenerationMatch=0` precondition so identical
  concurrent writes are skipped. Distinct hash-named objects are then **loaded (deduplicated)
  into `vcf_header_lines`** — the same final table as Naive, so the header text stays fully
  inspectable by the sanity check (§1.5).
- **`true` (command lines):** write naively inline with the associations → no small-object
  explosion, no dedup benefit lost (they were never going to dedup).

- **Pros:** cost-optimal footprint (the shared blob is stored once, not N times).
- **Pros (concurrency — the strongest structural argument):** content-addressing gives
  **lock-free cross-shard dedup**. Generation is sharded (§1.5, `GvsImportGenomes.wdl:186`) and
  shards are isolated (constraint #3), so no shard can see another's chunks. The BQ path solved
  this with a *synchronous shared-table* existence check (`doRowsExistFor`), the source of its
  ~1.6M-query storm (§5). Hybrid instead makes the **GCS object namespace the coordination
  point**: concurrent shards racing to write the same `<hash>.parquet` are resolved by the
  `ifGenerationMatch=0` precondition (first writer wins; the rest get a 412 and skip). So the
  shared blob is written **once** across all shards with **no locks, no round-trips, and no
  central table**. This holds regardless of shard count and sidesteps the anti-join's
  concurrent-writer gap (§6.6), because dedup happens at the storage layer, not via
  read-then-INSERT.
- **Cons:** two code paths in the writer; needs a bespoke content-addressed loader (doesn't fit
  the superpartition FOFN grouping in `DiscoverParquetFiles`).
- **Scaffold caveat:** the current `HYBRID` enum value (`VcfHeaderLineScratchCreator`, dormant
  behind the hardcoded strategy switch) is a *stub* — selecting it throws, because the GCS
  content-addressing half is not built. A complete Hybrid must actually store the blob and load
  its text into `vcf_header_lines`; until then Hybrid is not functional and would leave the shared
  header lines uninspectable. **Only Naive is fully implemented.**

### 3.3 Cost & concurrency comparison

At header data volumes **every option is dollar-trivial**; the real differentiators are
*transient storage footprint*, *concurrency model*, and *engineering complexity*. The current BQ
path's cost is an *operational* query storm rather than dollars.

Parameters (per batch of **N** samples; AoU ≈ 400k): **B** = shared blob size (~50 KB); **k** =
command-line chunks per sample (~1–3); **C** = size of each command-line chunk (~200 B, so
`C ≪ B`); **d_B** = distinct shared blobs (usually 1).

| Dimension | Current BQ Write API | Naive | Hybrid |
|---|---|---|---|
| Ingest-time BQ queries | **~2(k+1)·N ≈ 1.6M** | 0 | 0 |
| GCS transient bytes | 0 | **N·B** | N·k·C + d_B·B |
| GCS object count | 0 | N | N (tiny) + d_B |
| BQ dedup scan | 0 | **N·B** | N·k·C |
| Dollars @ 400k (B=50 KB) | negligible (operational storm) | ~$0.14 storage + ~$0.13 scan | <$0.01 |
| Concurrency dedup | synchronous shared-table check | anti-join, *sequential*-only (§6.6) | lock-free via object namespace |
| Complexity / risk | (baseline) | Low | Medium (bespoke blob loader) |

**Recommendation — and what shipped.** Implement **Naive first** (low risk, unblocks TSPS/AoU,
benchmarks the "naive" cost for AC1) — this is what PR 9403 delivered. Naive and Hybrid share the
same writer-schema and WDL wiring; only the dedup locus differs. **Hybrid is the eventual target**,
not for the (trivial) dollars but for the cleaner *concurrency* story: it makes the shared blob's
*"stored once"* a structural, lock-free, shard-count-independent invariant of the write path,
rather than a property a single after-the-fact anti-join scan has to reconstruct (and that scan is
only sequentially safe, §6.6). The staging is deliberate: ship Naive to unblock and to get the AC1
measurement, then move to Hybrid once the content-addressed loader is built.

---

## 4. Implementation work breakdown (Naive, with Hybrid notes)

### 4.1 Java — writer + schema
- `CreateVariantIngestFiles.java:226` — expand `headersRowSchema` to the full 4 columns
  (`sample_id`, `vcf_header_lines` text, `vcf_header_lines_hash`, `is_expected_unique`).
- `VcfHeaderLineScratchCreator.java:109-113` — flesh out the `PARQUET` branch to write the
  full record (text + hash + `is_expected_unique`); drop the write-time existence checks on
  this path.
- `HeaderParquetFileWriter.java:38` (`writeJson`) — extend to the full record shape.
- *(Hybrid)* route on `is_expected_unique`: content-address the shared blob, inline the
  unique lines.
- Reference template: commits `e852dcf27` / `2d48b4a84` (ref & ploidy Parquet writers) and
  `gvs_ingest_refactoring.md` (note VS-1803 is called out there as the pending header work).

### 4.2 WDL — `GvsImportGenomes.wdl` (implemented)
- Removed the `CannotLoadHeadersWithParquetIngest` guard.
- **Header parquet naming:** the Java writer now names the file
  `vcf_header_lines_scratch_<sampleId>.parquet` (was `header_file.parquet`) so it matches the
  regular-table discovery regex in `parse_and_group_files.py`
  (`/(prefix)_([0-9]+)...`). Without this the file would never be discovered.
- **Upload:** the generate task's per-sample upload block now guards the vet/ref/ploidy
  copies behind `load_vet_and_ref_ranges` and adds a `load_vcf_headers`-guarded copy of the
  header parquet to `.../vcf_header_lines_scratch/`. (`rm *.parquet` → `rm -f`.)
- **One combined discover/load/verify** now runs when
  `use_parquet_ingest && (load_vet_and_ref_ranges || load_vcf_headers)` — pulled out of the
  old `load_vet_and_ref_ranges`-only block so a **headers-only** run also loads. Prefixes:
  `regular = load_vcf_headers ? [sample_chromosome_ploidy, vcf_header_lines_scratch] :
  [sample_chromosome_ploidy]`, `superpartitioned = [vet, ref_ranges]` (always non-empty to
  avoid empty-array / argparse issues; extra prefixes are harmless — no files match). The
  loader/verifier need no changes beyond receiving these prefixes.
- **VerifyParquetLoading** now takes and forwards `--regular-table-prefixes` /
  `--superpartitioned-table-prefixes` (previously defaulted, which would flag header files as
  unverified).
- **Gate:** `ProcessVCFHeaders` now gates on
  `flatten([select_all([VerifyParquetLoading.done]), select_all(LoadDataViaBigQueryWriteAPI.done)])`
  — i.e. the Parquet header load for parquet runs, or the BQ Write API for BQ runs.
  `SetIsLoadedColumn` stays under `load_vet_and_ref_ranges` (headers-only runs don't mark
  samples loaded).
- **Phasing (per §1.5):** the supported operational model runs the workflow **twice** — a
  headers-only ingest (`load_vcf_headers = true`, `load_vet_and_ref_ranges = false`), manual
  inspection, then a data-only ingest. Each run localizes its own VCFs and cleans up its own
  parquet after loading, so no intra-run phase-aware deletion is needed. (Setting both
  booleans true in one run would load headers and data together and is technically possible,
  but it bypasses the manual gate and is **not** the supported model.)
- Validated with `womtool validate`. End-to-end coverage is the Terra `quickstart*` tests.

### 4.3 Python — scratch → final load
- `process_sample_vcf_headers.py` — change the `vcf_header_lines` INSERT to an **anti-join
  INSERT**: `INSERT ... SELECT s.hash, ANY_VALUE(s.text), ANY_VALUE(s.is_expected_unique)
  FROM scratch s LEFT JOIN vcf_header_lines t USING(vcf_header_lines_hash) WHERE s.text IS
  NOT NULL AND t.vcf_header_lines_hash IS NULL GROUP BY s.hash`. This dedups the naive-loaded
  scratch rows *and* makes the step idempotent across retries and incremental batches
  (satisfies constraints #4 and #6).
- The `sample_vcf_header` INSERT **also** needs the same anti-join treatment, keyed on
  (`sample_id`, `vcf_header_lines_hash`). *(A blind INSERT is not idempotent here — a task retry
  or re-ingest would duplicate the (sample_id, hash) associations. See §6.)*
- Anti-join INSERT is used rather than `MERGE` to match GVS convention (GVS uses `MERGE`
  nowhere) and sidestep mutating-DML concerns; `INSERT` also never touches a streaming buffer.
- `load_parquet_to_bq.py` — no change needed for Naive (raw load into scratch), but the
  scratch loads rely on deterministic BQ load job IDs for retry-idempotency
  (`load_parquet_to_bq.py:166`, `_make_job_id`) — see §6.
- Hybrid would need a content-addressed loader.

### 4.4 Table schema
- No schema change to the three tables. `vcf_header_lines_scratch` stays as the load
  target for Naive. *(Hybrid could drop scratch for the shared blob on the Parquet path.)*

---

## 5. AC1 — cost of a naive header load into BQ

"Naive" = the implemented strategy (§3.1). Measure and compare against the current BQ Write API
path.

**Naive (Parquet) cost components:**
- transient scratch storage ≈ N samples × (blob + command-line text) bytes;
- the BQ Parquet load-job bytes;
- the dedup (anti-join INSERT) query bytes scanned (dominated by the large text column).

**Current BQ Write API path baseline (for comparison):**
- `doRowsExistFor` runs a `SELECT COUNT(*)` per lookup (`BigQueryUtils.java:455`), called
  **twice per chunk** (scratch + non-scratch) × ~2 chunks/sample ≈ **~4 BQ query jobs per
  sample** — on the order of **~1.6M query jobs for a 400k-sample callset**, serialized
  inside each ingest task, plus the streaming inserts.

Escaping that per-sample query storm is a primary motivation; the naive Parquet approach
trades it for transient storage + one dedup query.

**Suggested measurement harness:** `GvsQuickstartVcfIntegration.wdl` (smallest end-to-end
path with `load_vcf_headers = true`); capture bytes/slot-ms from the `cost_observability`
table and BQ job stats for both paths at a fixed sample count, then extrapolate. Empirically
measure the drivers from one real gVCF header — **B** and **k** (they dominate every formula) and
**d_B** via `SELECT COUNT(DISTINCT vcf_header_lines_hash) ... WHERE is_expected_unique = false`.

---

## 6. Idempotency

**Goal:** re-running header loading converges to the same final tables — no duplicate rows,
no partial corruption — under three triggers: (a) a WDL **task retry** (preemption /
`maxRetries`), (b) a Cromwell **workflow resume / re-run**, and (c) a deliberate
**re-ingest** of the same batch.

### 6.1 Existing idempotency, and the gap the Parquet path leaves

The system **was designed with idempotency in mind — this is not accidental.** The BQ ingest
path tracks a per-sample `HEADERS_LOADED` status in `sample_load_status` (`LoadStatus.java`):
`CreateVariantIngestFiles` checks it and **skips** header writes for a sample whose headers are
already loaded, then writes the status on success (`CreateVariantIngestFiles.java:338-354`,
`:466-468`). Re-ingesting an already-loaded sample is a deliberate no-op. The BQ path also
dedups header *text* at write time by hash (`VcfHeaderLineScratchCreator`), so the
scratch→final blind INSERT sees one non-null-text row per hash.

Two narrower gaps remain, both specific to the **Parquet** path:

- **It bypasses the per-sample status logic.** `CreateVariantIngestFiles.java:356-361`
  short-circuits the existence checks for `PARQUET` ("operate as though they don't exist"), and
  `HEADERS_LOADED` is never written on the Parquet path (`shouldWriteVCFHeadersLoadedStatusRow`
  is set only on the BQ branch, `:353`). So Parquet generation does not skip already-loaded
  samples. *(Aside: the Parquet path not writing `HEADERS_LOADED` at all is a separate gap —
  see §8.)*
- **The promotion step is not self-contained.** The original scratch→final blind INSERT relied
  on the scratch DELETE having run to avoid re-inserting; and because the Parquet path does no
  write-time hash dedup, its scratch holds full text for *every* (sample, hash), which a blind
  INSERT would duplicate in `vcf_header_lines`.

The anti-join INSERT (§6.2) closes both: it dedups by key and refuses to re-insert rows already
in the target, so the promotion is deduplicating and idempotent regardless of the scratch
DELETE or the (absent) Parquet per-sample skip.

### 6.2 Where idempotency must be enforced per layer

The Parquet path deliberately moves dedup/idempotency *off* the ingest write (no BQ
round-trips there) and onto the scratch → final load. Each layer needs a keyed, re-runnable
operation:

| Layer | Non-idempotent form | Idempotent form |
|---|---|---|
| Parquet write (offline, per sample) | (already safe — writes local files only; re-run overwrites) | deterministic file contents |
| GCS upload of header parquet | `cp` (overwrite — safe, identical bytes) | idempotent by nature; Hybrid blob uses `ifGenerationMatch=0` |
| BQ load scratch ← parquet | `WRITE_APPEND` retried → double-load | **deterministic load job ID** (`load_parquet_to_bq.py:166`) so BQ dedups the retried job |
| Load scratch → `vcf_header_lines` | blind INSERT | **anti-join INSERT** on `vcf_header_lines_hash` |
| Load scratch → `sample_vcf_header` | blind INSERT | **anti-join INSERT** on (`sample_id`, `vcf_header_lines_hash`) |
| Scratch cleanup | — | DELETE is naturally idempotent (re-delete = no-op) |

Key point: with both final-table populations as **anti-join INSERTs keyed on their natural
keys**, the whole chain becomes idempotent *independently* of the scratch-delete step —
removing the fragile dependency in §6.1.

### 6.3 Idempotency by strategy

- **Naive (naive + anti-join INSERT):** idempotent once both loads are anti-joined (§6.2)
  and scratch loads use deterministic job IDs. Note scratch can accumulate rows across
  retries; the anti-join skips keys already in the target (and `GROUP BY`/`DISTINCT` collapses
  within-batch dups), which is exactly why an anti-join INSERT — not a blind
  `INSERT ... SELECT DISTINCT`, which would still re-insert rows already loaded — is the right
  primitive here.
- **Hybrid (content-addressed blob):** strongest form — the expensive header **text** is
  idempotent at the *storage* layer (path = hash, identical bytes, `ifGenerationMatch=0`
  skips re-writes), so re-runs are no-ops for the blob. The BQ association / `vcf_header_lines`
  loads still need the keyed anti-join INSERT / skip-existing, same as Naive.

### 6.4 Recommendation

Make idempotency structural, not incidental:
1. Both final-table populations become **anti-join INSERT** statements keyed on natural keys.
2. Rely on **deterministic BQ load job IDs** for scratch-load retry safety.
3. Keep the scratch DELETE as cleanup, but the design must be correct **without** depending
   on it having run.

This is small extra work on top of Naive (one of the two merges was already going to
change for dedup) and it satisfies constraint #6 for both strategies.

### 6.5 Table persistence & incremental ingest

Callsets are ingested in batches, so it matters which header tables persist:

- **`vcf_header_lines` + `sample_vcf_header` are the incremental accumulator.** Each batch
  merges its headers into them and they grow across batches. They must persist for the
  callset's life.
- **`vcf_header_lines_scratch` is throwaway.** `clean_up_scratch_table` empties it after every
  merge, so it holds no long-term data; deleting its *rows* is normal. Dropping the *table*
  breaks the next run only until it is recreated (`DiscoverParquetFiles` queries it).
- **Subtlety:** header "skip already-loaded" detection in `DiscoverParquetFiles` reads the
  *transient scratch* table (regular prefix `vcf_header_lines_scratch`), unlike vet/ref/ploidy
  which check their persistent data tables. So cross-batch skip is effectively a no-op for
  headers — harmless only because each batch generates parquet for its new samples and old
  parquet is deleted after loading. Incremental **correctness rests on the anti-join INSERT +
  the persistent tables**, not on skip-detection.

### 6.6 Re-run idempotency vs. concurrent-writer safety (review note)

The anti-join INSERT makes the promotion idempotent across **sequential** re-runs (task retry,
Cromwell resume, deliberate re-ingest). It does **not** make two promotions writing
**concurrently** safe: two `INSERT ... SELECT ... LEFT JOIN t ... WHERE t.key IS NULL`
statements can both observe a key as absent and both insert it, because BigQuery does **not**
conflict-abort concurrent `INSERT`s the way it does `UPDATE`/`DELETE`/`MERGE` under optimistic
concurrency. The result would be a duplicate row.

This is **not** a live risk in the current wiring: `ProcessVCFHeaders`
(`GvsImportGenomes.wdl:335`) is a **single, non-scattered** task that runs both anti-join
INSERTs exactly once per workflow run. The per-shard scatter (`LoadParquetFilesToBQ`,
`GvsImportGenomes.wdl:281-290`) only *appends* raw rows into `vcf_header_lines_scratch`
(retry-idempotent via deterministic BQ load job IDs, §6.2) and never runs the anti-join — so no
two shards race on the promotion. The **BQ path uses the same single-task promotion**, so this
is not a Parquet-specific regression.

The race would only appear if **two ingest runs promoted into the same dataset concurrently**.
The WDL guarantees only the *within-run* single-task promotion above; whether two
`GvsImportGenomes` runs can overlap on one dataset is an **orchestration-layer question** the
WDL itself cannot answer. **Resolved:** the bulk-ingest / AoU owner (mcovarr) has confirmed this
must never happen going forward — the rest of the ingest system is not designed for that
concurrency and would almost certainly enter a bad state if two runs shared a dataset. So the
anti-join's *sequential* idempotency is sufficient for the supported operational model, and no
dataset-level lock is required. (Were that invariant ever relaxed, the anti-join alone would be
insufficient and the promotion would need serializing or a dedup pass; Hybrid's
content-addressed blob dedup (§3) sidesteps the concern entirely for the expensive shared blob
regardless.)
*(Raised in review by mcovarr on the VS-1968 doc; concurrency invariant confirmed by mcovarr.)*

---

## 7. AC2 — header-analysis queries WIP based on VS-1215

The consistency checks are captured in **VS-1215** (*Create WDL to analyze VCF headers and
create a small report*). Those queries — from the manual analysis of the first Echo header
ingest — read the final `vcf_header_lines` / `sample_vcf_header` / `sample_info` tables and
fall into three groups. Because those tables are identical across load paths (§3), the queries
are **path-agnostic** and validate a Parquet-loaded dataset the same as a BQ-loaded one.

**1. Integrity invariants** (should hold for any correct header load):
- distinct `sample_id` in `sample_vcf_header` == the loaded (non-control) sample count in
  `sample_info` — i.e. every sample has header data. *(Echo: 414,830 ✓.)*
- every sample has ≥1 `is_expected_unique = true` chunk (its per-sample processing info is
  present). *(Echo: 414,830 ✓.)*
- referential integrity: every `sample_vcf_header.vcf_header_lines_hash` resolves to a row in
  `vcf_header_lines` (no orphan hashes).

**2. Shared-blob distribution** — the distinct `is_expected_unique = false` blobs and how many
samples carry each. Ideally one; more than one indicates distinct delivery batches
(informational, not necessarily an error). *(Echo: 6 distinct shared blobs.)*

**3. DRAGEN-version report (the "useful info")** — group samples by the DRAGEN software-version
string in their `is_expected_unique = true` chunk (`DRAGENCommandLine=<ID=dragen,Version="SW: …"`).
The final VS-1215 query joins `vcf_header_lines` × `sample_vcf_header` × `sample_info` against a
small `dragen_versions` CTE and counts samples per version. *(This is exactly how Echo's mixed
processing was found: four versions — `05/07.021.408.3.4.12` and `05/07.021.604.3.7.8`.)*

**Relation to this spike:** these are the concrete form of the §1.5 input-gVCF sanity check,
and they double as the Parquet consistency check — run them on a Parquet-loaded dataset and
confirm the invariants hold and the version breakdown matches expectations. If an explicit
A/B check against a BQ-loaded dataset is wanted, the same tables can be compared directly
(distinct-hash sets, the (sample_id, hash) set, and reconstructed per-sample text).

---

## 8. Open questions / risks

- **`is_loaded` would not be set for Parquet callsets that load headers** (found by static
  tracing during the VS-1968 spike; fix applied under the VS-1803 implementation — no
  dedicated bug ticket yet — but **not** yet verified end-to-end, see `quickstart*` note). The
  Parquet path writes no `sample_load_status` rows: all `shouldWrite*StatusRow` flags are set
  only on the BQ branch (`CreateVariantIngestFiles.java:318/334/353`) and the Parquet branch
  bypasses them (`:356-361`). Separately, `sample_info.is_loaded` is set post-load in bulk by
  `SetIsLoadedColumn` from the `samples_with_all_data` view, which JOINs `samples_with_header_data`
  when headers are in play. That header view read the **transient** `vcf_header_lines_scratch`
  table (emptied by `ProcessVCFHeaders`) UNION the `HEADERS_LOADED` status — **neither of which
  survives for Parquet**, so the view went empty and `is_loaded` was never set.
  - **Fix applied (review-endorsed):** point `samples_with_header_data` at the durable
    `sample_vcf_header` table instead of transient scratch
    (`GvsImportGenomes.wdl:1054-1058`), `UNION DISTINCT` the `HEADERS_LOADED` status
    (`:952-960`, `:1036-1058`). This mirrors `samples_with_reference_data` /
    `samples_with_variant_data` (data-presence source `UNION` `sample_load_status`) — exactly
    the structure mcovarr recommended in review (VS-1968 doc): the Parquet side derives presence
    from `sample_vcf_header`, and the `sample_load_status` `UNION` preserves legacy BQ-loaded
    datasets. It is path-agnostic and removes reliance on `HEADERS_LOADED` for the Parquet path.
  - **Safe for the BQ path:** the view is a shared BQ/Parquet artifact, but the change doesn't
    regress BQ — the `HEADERS_LOADED` status clause (still UNION'd in, gated on the
    `sample_load_status` table existing) covers the during-ingest window before the
    scratch→final merge, and post-merge `sample_vcf_header` is populated for BQ too, so no BQ
    samples drop out of the view. No BQ code-path change is required.
  - **The Parquet path deliberately does not write `HEADERS_LOADED`** — per mcovarr's review
    (VS-1968 doc), header-load presence is determined *by the view* (as vet/ref/ploidy already do),
    not by a per-sample status row. Writing `HEADERS_LOADED` on the Parquet path was considered
    and **rejected**.
  - **Verify in the Terra `quickstart*` run** (traced statically; not executed here).
- **Small-object cost (Hybrid):** how many distinct `is_expected_unique = true` chunks
  arise at AoU scale? Determines whether content-addressing the unique lines would ever be worth
  it (Hybrid leaves them inline for exactly this reason).
- **Dedup-scan cost at scale (Naive):** the anti-join INSERT scans the text column;
  confirm this is cheaper than the query storm it replaces.
- **`billing_project_id` not threaded to the loader** (`GvsImportGenomes.wdl:1248`,
  VS-1955) — verify header loads don't need it.
- **Scratch retention:** Naive keeps scratch; confirm `clean_up_scratch_table`
  (`process_sample_vcf_headers.py`) still applies cleanly after the anti-join INSERT.
- **Concurrent promotions into one dataset (§6.6): resolved.** The anti-join is only
  *sequentially* idempotent, which is sufficient because the bulk-ingest / AoU owner (mcovarr)
  has confirmed two `GvsImportGenomes` runs must never overlap on the same dataset going forward
  (the ingest system is not designed for that concurrency). No dataset-level lock is required; if
  that invariant were ever relaxed, serialize the promotion or move to Hybrid's
  content-addressed blob dedup.

---

## 9. File / line reference index

| Area | Location |
|---|---|
| Parquet header writer + strategy enum | `VcfHeaderLineScratchCreator.java` (`HeaderParquetStrategy`, `:109-113`) |
| BQ write-time dedup | `VcfHeaderLineScratchCreator.java:93-108` |
| Already-loaded existence checks | `VcfHeaderLineScratchCreator.java:43-52`, `BigQueryUtils.java:455` |
| Header chunking | `CreateVariantIngestFiles.java:524-537` |
| Parquet header schema | `CreateVariantIngestFiles.java:226` |
| Header Parquet file writer | `HeaderParquetFileWriter.java:38` |
| Table schemas / creation | `GvsAssignIds.wdl:30-32`, `GvsAssignIds.wdl:96-129` |
| Scratch→final merge / anti-join INSERT | `process_sample_vcf_headers.py` |
| Discover / load / gate | `GvsImportGenomes.wdl:278-279`, `:288`, `:335` |
| `is_loaded` view fix | `GvsImportGenomes.wdl:952-960`, `:1036-1058` |
| Loader (table name from FOFN) | `load_parquet_to_bq.py:16,46-52,166` |
