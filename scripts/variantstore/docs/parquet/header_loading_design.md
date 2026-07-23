# Loading VCF Headers in the Parquet Code Path

**Tickets:**
- **VS-1968** — *[Spike] How to load headers in Parquet code path* — deliverables are this
  design doc and the consistency checks (§8).
- **VS-1803** — *Header loading not supported on Parquet branch* — the **implementation** (the
  Java/Python/WDL changes and the `is_loaded` view fix described here belong to this ticket).

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

### 1.4 Header extraction: single invocation vs. the gated model

`CreateVariantIngestFiles` already localizes each VCF once (`GvsImportGenomes.wdl:551-553`)
and runs a **single** call (`GvsImportGenomes.wdl:560`) that already accepts
`--enable-vcf-headers` (`GvsImportGenomes.wdl:575`), so header extraction *can* piggyback on
the existing VCF read with no extra localization. **However**, the operational model (§1.5)
deliberately loads headers as a *separate, gated* phase before vet/ref/ploidy, which is in
tension with folding everything into one pass. The Parquet architecture reconciles this
because generation and BQ loading are already decoupled — see §1.5. Either way the remaining
work is finishing the plumbing, not adding an ingest tool.

### 1.5 Operational model: phased, gated loading

**Purpose of the header phase.** The header-only ingest is not a check that header *loading*
worked — it is an early **sanity check of the input gVCFs themselves**. Loading just the
headers is a cheap way to inspect each input's metadata (reference version, contig list,
sample names, expected INFO/FORMAT/FILTER definitions, pipeline provenance) and confirm the
cohort is well-formed and mutually consistent *before* committing to the expensive vet/ref
ingest. This purpose must be preserved in the Parquet path.

**AoU today (BQ path).** Header loading is run as its **own** ingest first
(`load_vcf_headers = true`, `load_vet_and_ref_ranges = false`), the loaded headers are then
**manually inspected** to sanity-check the input gVCFs, and **only if that passes** is a
second ingest run for vet / ref_ranges / sample_chromosome_ploidy. Three phases: *load
headers → sanity check inputs → load everything else*. On the BQ path this localizes the VCFs
**twice** (once per ingest run).

**TSPS (assumption — undecided).** Likely the same three phases, but **automated** (the
manual sanity check becomes an automated gate).

**How this maps onto Parquet.** In the Parquet path, *generation* (localize VCF → write
parquet to GCS) is already decoupled from *loading* (`DiscoverParquetFiles` →
`LoadParquetFilesToBQ`). That opens two ways to honor the gate:

- **(A) Generate-all-once, phased load.** One `CreateVariantIngestFiles` pass writes *all*
  parquet (headers + vet + ref + ploidy) with a single localization; then load **only** the
  header parquet, gate on the sanity check, and load the rest afterward from the parquet
  already in GCS. Single localization *and* a real gate.
  *Cost:* vet/ref/ploidy parquet is **generated before** the gate, so a failed header check
  wastes that generation compute (and transient GCS storage).
- **(B) Phase-separated generation.** A headers-only generation+load first (cheap), gate,
  then a second generation+load for the rest — **re-localizing** the VCFs (as AoU does today
  on the BQ path). *Cost:* pays localization + generation twice, but spends **no** vet/ref
  compute until the gate passes.

**Trade-off — sharpened by the purpose.** Because the gate exists specifically to *avoid
expensive work on bad inputs*, approach (A) partly defeats it: (A) generates the vet/ref/ploidy
parquet **before** the gate, so a caught bad input still pays the full generation cost — only
the vet/ref *BQ load* is saved. In fact the "avoid re-localizing in a second task" wish
(from VS-1803, not this spike) and the header phase's purpose are **fundamentally in
conflict**: validating inputs
before doing expensive variant work requires *either* re-localizing (B) *or* doing that
expensive work first (A). We can have at most one.

Concretely, per-sample the costs are: **localize** (copy the gVCF; needed for headers *and*
variants), cheap **header parse**, and expensive **vet/ref parse**.
- **(A)** localizes once but always does the vet/ref parse up front → wastes vet/ref parse on
  a gate failure.
- **(B)** localizes twice on success but does the vet/ref parse only after the gate → wastes
  almost nothing on failure.

So (A) wins only when localization dominates *and* input failures are rare; the stated
purpose (the phase exists because bad inputs are a real risk) leans toward **(B)**, or toward
(A) only if TSPS confirms failures are rare and localization is the dominant cost.
**Recommendation:** leave open until TSPS fixes its model, but note the conflict above rather
than assuming the re-localization can simply be optimized away.

**Implication either way:**
- The header load must complete and be validated as a *distinct, gated step* before the
  vet/ref/ploidy load (see the WDL gate in §4.2); for approach (A) the parquet cleanup
  (`DeleteParquetFiles`, which runs by default right after loading — plus the off-by-default
  14-day GCS lifecycle rule as a backstop) must **not** remove the vet/ref parquet before the
  post-gate load.
- The sanity check reads header **content**, so the Parquet header load must persist the
  **full header text**, not just hashes. Storing it *deduplicated* does not hurt inspection:
  the check reconstructs any sample's header by joining `sample_vcf_header` (sample → hash) to
  `vcf_header_lines` (hash → text). Every design in §3 ends with those same two tables, so all
  of them keep the header lines (contigs, reference, INFO/FORMAT) fully inspectable — the only
  way to lose that is an *incomplete* implementation that persists hashes without ever loading
  the text into `vcf_header_lines` (see the caveat on the current content-addressing scaffold
  in §3 / §7.3).

---

## 2. Key findings / constraints

1. **The merge does not dedup** (§1.2). Naive Parquet writing requires a dedup step to be
   added.
2. **`is_expected_unique` distinguishes the shared blob from per-sample command lines** —
   `false` for the joined non-command-line blob (shared, dedup-worthy), `true` for the
   command-line chunks (expected to differ per sample). This flag is **already used** by the
   AoU header sanity check: its query filters `is_expected_unique = TRUE` to isolate the
   per-sample command lines and extract the DRAGEN version for validation
   (`AOU_DELIVERABLES.md`). Dedup designs can *additionally* route storage on it (Option 3).
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
   `sample_load_status` (§7.1); it is not accidental. The gap the Parquet work must close is
   narrower: the Parquet path bypasses that per-sample status check and the promotion step's
   blind INSERT does not dedup — both addressed by the anti-join INSERT (§7).
7. **Loading is phased and gated — as an input-gVCF sanity check** (§1.5). The header-only
   phase exists to validate the *input gVCFs* (reference, contigs, samples, expected fields)
   before the expensive vet/ref ingest, so it must remain a distinct, gated phase in the
   Parquet path too. Consequences: (a) the header load must stand alone and gate the data
   load; (b) it must persist **full header text** so the check can inspect content (a hard
   requirement on Option 3's content-addressing); (c) for approach (A) the parquet cleanup
   must be phase-aware so vet/ref parquet survives until its post-gate load; (d) this gated
   re-run pattern makes constraint #6 (idempotency) more pressing, not less.

---

## 3. Design options

### Option 1 — Naive write, dedup in the load SQL (low risk)
- Writer: emit full text + `is_expected_unique` for every chunk, every sample; drop the
  existence queries.
- Load: change `process_sample_vcf_headers.py:21` from a blind INSERT to an **anti-join
  INSERT** (`INSERT ... SELECT ... LEFT JOIN <target> USING(key) WHERE <target>.key IS NULL`,
  with `GROUP BY`/`DISTINCT` on the source), so it dedups *and* is idempotent across batches
  (satisfies constraint #4). Anti-join INSERT rather than `MERGE` to match GVS convention
  (GVS uses `MERGE` nowhere) and avoid mutating-DML concerns.
- **Cost:** max (transient) scratch bloat; one large dedup query scanning the text
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
- `false` (shared blob): content-address it (Option 2) → stored once, then **loaded
  (deduplicated) into `vcf_header_lines`** — the same final table as Option 1, so the header
  text remains fully inspectable by the sanity check.
- `true` (command lines): write naively inline with associations → no small-object
  explosion, no dedup benefit lost (they were never going to dedup).
- **Pros:** cost-optimal; uses the schema flag as intended.
- **Cons:** two code paths in the writer.
- **Scaffold caveat:** the current `HYBRID` code (`VcfHeaderLineScratchCreator`, dormant behind
  the hardcoded strategy switch) is a *stub* — it writes the blob's association only and drops
  the text with a TODO. A complete Option 3 must actually store the blob and load its text into
  `vcf_header_lines`; until then `HYBRID` is not functional and would leave the shared header
  lines uninspectable. Only Option 1 is fully implemented.

### Option 4 — Pre-load local consolidation
- Dedup header Parquet by hash in a consolidation step before load (git precedent:
  "Experimenting with local consolidation and deduping of parquet data").
- **Cons:** already paid N× GCS writes; only saves BQ-side cost. Weaker than 2/3.

**Recommendation:** implement **Option 1 first** (low risk, unblocks TSPS/AoU, benchmarks
the "naive" cost for AC1), then evaluate **Option 3** if the naive cost is unacceptable.
Options 1 and 3 share the same writer-schema and WDL wiring work; only the dedup locus differs.

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
- **Phasing (per §1.5):** the AoU-style gate is achieved by running the workflow **twice**
  (headers-only, sanity-check, then data-only). The single-run intra-workflow gate
  (approach A with a split scatter) and phase-aware parquet deletion are **not** implemented —
  in a combined run headers and data load together.
- Validated with `womtool validate`. End-to-end coverage is the Terra `quickstart*` tests.

### 4.3 Python — scratch → final load
- `process_sample_vcf_headers.py:21` — change the `vcf_header_lines` INSERT to an **anti-join
  INSERT**: `INSERT ... SELECT s.hash, ANY_VALUE(s.text), ANY_VALUE(s.is_expected_unique)
  FROM scratch s LEFT JOIN vcf_header_lines t USING(vcf_header_lines_hash) WHERE s.text IS
  NOT NULL AND t.vcf_header_lines_hash IS NULL GROUP BY s.hash`. This dedups the naive-loaded
  scratch rows *and* makes the step idempotent across retries and incremental batches
  (satisfies constraints #4 and #6).
- `process_sample_vcf_headers.py:24` — the `sample_vcf_header` INSERT **also** needs the same
  anti-join treatment, keyed on (`sample_id`, `vcf_header_lines_hash`). *(Correction: this is
  not idempotent as a blind INSERT — a task retry or re-ingest would duplicate the
  (sample_id, hash) associations. See §7.)*
- Anti-join INSERT is used rather than `MERGE` to match GVS convention (GVS uses `MERGE`
  nowhere) and sidestep mutating-DML concerns; `INSERT` also never touches a streaming buffer.
- `load_parquet_to_bq.py` — no change needed for Option 1 (raw load into scratch), but the
  scratch loads rely on deterministic BQ load job IDs for retry-idempotency
  (`load_parquet_to_bq.py:166`, `_make_job_id`) — see §7.
- Option 2/3 would need a content-addressed loader.

### 4.4 Table schema
- No schema change to the three tables. `vcf_header_lines_scratch` stays as the load
  target for Option 1. *(Option 2/3 could drop scratch for the Parquet path.)*

---

## 5. AC1 — cost of a naive header load into BQ

"Naive" = Option 1. Measure and compare against the current BQ Write API path.

**Option 1 (Parquet) cost components:**
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
3. **GCS transient storage** — held until `DeleteParquetFiles` removes it after loading (the
   14-day GCS lifecycle rule is only an off-by-default backstop).
4. **GCS object count** — drives discovery `ls` cost and BQ small-file load overhead.
5. **BQ load jobs** — free in dollars, but subject to quota + small-file overhead.
6. **BQ query bytes scanned for dedup** — the on-demand dollar driver.
7. **Final BQ storage** — deduped; ~identical and negligible across all options.

Representative unit prices (US, *verify current*): on-demand query **$6.25 / TiB scanned**;
GCS Standard **~$0.02 / GB-month**; BQ Storage Write API **~$0.025 / GB** (streaming); BQ
batch **load jobs are free**.

### 6.2 Comparison table (dominant term per cell)

| Dimension | Current BQ Write API | Opt 1 Naive+anti-join INSERT | Opt 2 Content-addr (all) | Opt 3 Hybrid | Opt 4 Pre-load consolidate |
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
  (every sample's copy of the blob) for up to 14 days, plus one dedup (anti-join INSERT) scan
  of **N·B** once. Both are ~N·B — cheap in absolute dollars (§6.4) but the largest
  storage/scan footprint of the Parquet options. Zero ingest queries; lowest complexity.
- **Option 2 (content-addressed, all chunks)** — eliminates N·B: the blob is stored once,
  dedup is free via the object namespace, no dedup scan. But `is_expected_unique = true`
  lines can explode to **~N tiny objects**, which hits (a) `ls` discovery cost/time, (b) BQ
  small-file load overhead, and (c) GCS Class-A op counts from N precondition writes.
  Lowest bytes, worst object-count profile, highest complexity (bespoke loader).
- **Option 3 (hybrid)** — keeps the blob content-addressed (stored once, no dedup scan for
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
financial: Option 1's downside is a ~20 GB transient footprint (auto-cleaned by
`DeleteParquetFiles` after loading, with the 14-day lifecycle rule as an off-by-default
backstop) and one ~20 GB scan; Options 2/3 shrink that to ~80 MB at the price of complexity
(Opt 2 also trades it for an ~N-tiny-object problem). This argues for **Option 1 first**
(cheapest to build, only real cost is transient storage that cleanup already handles) and
**Option 3** only if the transient footprint or object count becomes a problem at full AoU
scale.

### 6.5 How to measure each (extends the AC1 harness)

Reuse the `GvsQuickstartVcfIntegration.wdl` harness at a fixed sample count and extrapolate
by N:

- **Ingest queries:** count BQ jobs per ingest step (`cost_observability` / INFORMATION_SCHEMA.JOBS);
  confirm 0 for all Parquet options vs ~2(k+1)·N for the BQ path.
- **GCS bytes/objects:** `gcloud storage du` and object counts on the parquet output dir
  *before* `DeleteParquetFiles` runs.
- **Dedup scan:** `totalBytesProcessed` / `totalSlotMs` from the anti-join INSERT job.
- **Load:** job count + per-job duration from `load_parquet_to_bq.py` `stats.json`.
- **Empirically measure the drivers** from one real gVCF header: **B** and **k** (they
  dominate every formula), and **d_B** via
  `SELECT COUNT(DISTINCT vcf_header_lines_hash) ... WHERE is_expected_unique = false`.

---

## 7. Idempotency

**Goal:** re-running header loading converges to the same final tables — no duplicate rows,
no partial corruption — under three triggers: (a) a WDL **task retry** (preemption /
`maxRetries`), (b) a Cromwell **workflow resume / re-run**, and (c) a deliberate
**re-ingest** of the same batch.

### 7.1 Existing idempotency, and the gap the Parquet path leaves

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
  see §9.)*
- **The promotion step is not self-contained.** The original scratch→final blind INSERT relied
  on the scratch DELETE having run to avoid re-inserting; and because the Parquet path does no
  write-time hash dedup, its scratch holds full text for *every* (sample, hash), which a blind
  INSERT would duplicate in `vcf_header_lines`.

The anti-join INSERT (§7.2) closes both: it dedups by key and refuses to re-insert rows already
in the target, so the promotion is deduplicating and idempotent regardless of the scratch
DELETE or the (absent) Parquet per-sample skip.

### 7.2 Where idempotency must be enforced per layer

The Parquet path deliberately moves dedup/idempotency *off* the ingest write (no BQ
round-trips there) and onto the scratch → final load. Each layer needs a keyed, re-runnable
operation:

| Layer | Non-idempotent form | Idempotent form |
|---|---|---|
| Parquet write (offline, per sample) | (already safe — writes local files only; re-run overwrites) | deterministic file contents |
| GCS upload of header parquet | `cp` (overwrite — safe, identical bytes) | idempotent by nature; Opt 2/3 blob uses `ifGenerationMatch=0` |
| BQ load scratch ← parquet | `WRITE_APPEND` retried → double-load | **deterministic load job ID** (`load_parquet_to_bq.py:166`) so BQ dedups the retried job |
| Load scratch → `vcf_header_lines` | blind INSERT | **anti-join INSERT** on `vcf_header_lines_hash` |
| Load scratch → `sample_vcf_header` | blind INSERT | **anti-join INSERT** on (`sample_id`, `vcf_header_lines_hash`) |
| Scratch cleanup | — | DELETE is naturally idempotent (re-delete = no-op) |

Key point: with both final-table populations as **anti-join INSERTs keyed on their natural
keys**, the whole chain becomes idempotent *independently* of the scratch-delete step —
removing the fragile dependency in §7.1.

### 7.3 Idempotency by option

- **Option 1 (naive + anti-join INSERT):** idempotent once both loads are anti-joined (§7.2)
  and scratch loads use deterministic job IDs. Note scratch can accumulate rows across
  retries; the anti-join skips keys already in the target (and `GROUP BY`/`DISTINCT` collapses
  within-batch dups), which is exactly why an anti-join INSERT — not a blind
  `INSERT ... SELECT DISTINCT`, which would still re-insert rows already loaded — is the right
  primitive here.
- **Option 2 / 3 (content-addressed):** strongest form — the expensive header **text** is
  idempotent at the *storage* layer (path = hash, identical bytes, `ifGenerationMatch=0`
  skips re-writes), so re-runs are no-ops for the blob. The BQ association / `vcf_header_lines`
  loads still need the keyed anti-join INSERT / skip-existing, same as Option 1.
- **Option 4 (consolidation):** the consolidation output must be deterministic (stable
  file naming by hash) for the downstream load to be idempotent; otherwise re-runs produce
  new files that re-load.

### 7.4 Recommendation

Make idempotency structural, not incidental:
1. Both final-table populations become **anti-join INSERT** statements keyed on natural keys.
2. Rely on **deterministic BQ load job IDs** for scratch-load retry safety.
3. Keep the scratch DELETE as cleanup, but the design must be correct **without** depending
   on it having run.

This is small extra work on top of Option 1 (one of the two merges was already going to
change for dedup) and it satisfies constraint #6 for every option.

### 7.5 Table persistence & incremental ingest

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

---

## 8. AC2 — header-analysis queries WIP based on VS-1215

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

## 9. Open questions / risks

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
  - **Fix applied:** point `samples_with_header_data` at the durable
    `sample_vcf_header` table instead of transient scratch (`GvsImportGenomes.wdl`). This is
    path-agnostic and removes the reliance on `HEADERS_LOADED` (aligns with the code's note
    that `sample_load_status` "may no longer be needed").
  - **Safe for the BQ path:** the view is a shared BQ/Parquet artifact, but the change doesn't
    regress BQ — the `HEADERS_LOADED` status clause (still UNION'd in) covers the during-ingest
    window before the scratch→final merge, and post-merge `sample_vcf_header` is populated for
    BQ too, so no BQ samples drop out of the view.
  - *Alternative (not chosen):* have the Parquet path write `HEADERS_LOADED` — but only at a
    post-merge step (writing it during generation would mark "loaded" before the data is in
    BigQuery). Left as an option if per-sample header status is wanted for parity.
  - **Verify in the Terra `quickstart*` run** (traced statically; not executed here).
- **Small-object cost (Option 2/3):** how many distinct `is_expected_unique = true` chunks
  arise at AoU scale? Determines whether content-addressing the unique lines is viable.
- **Dedup-scan cost at scale (Option 1):** the anti-join INSERT scans the text column;
  confirm this is cheaper than the query storm it replaces.
- **`billing_project_id` not threaded to the loader** (`GvsImportGenomes.wdl:1248`,
  VS-1955) — verify header loads don't need it.
- **Scratch retention:** Option 1 keeps scratch; confirm `clean_up_scratch_table`
  (`process_sample_vcf_headers.py:28-31`) still applies cleanly after the anti-join INSERT.

---

## 10. File / line reference index

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
