-- =====================================================================================
-- Step 6 -- patch the delivered participant mapping table
-- =====================================================================================
--
-- Removes from `foxtrot.vid_to_participant_mapping_2026_08_28` every `(vid, person_id)`
-- entry named in `foxtrot.exclusions` (Step 5, 2,344,756,787 rows), writing the result as
-- the next dated version of the table. The correction is purely subtractive: the new
-- predicate is the old one plus `GQ != 0` plus `FT`, so a corrected `person_ids` array is
-- always a subset of the delivered one.
--
--
-- =====================================================================================
-- RUNBOOK -- the order to actually run things in. Everything else in this file is either
-- a diagnostic that has already been run or the reasoning behind one of these lines.
-- =====================================================================================
--
--   #   what                                       cost     status
--   --  -----------------------------------------  -------  --------------------------
--   P1  delivered table's physical shape           free     DONE 2026-09-05
--   P2  profile the exclusion sets                 ~$0.30   DONE 2026-09-05
--   P3  reconcile the VID sets                     ~$0.40   DONE 2026-09-05
--   P3b fixup VIDs' excluded pair count            ~$0.30   DONE 2026-09-05  37,702
--   P3c confirm the 993 overlap                    free     DONE 2026-09-05
--   P4  inspect the extreme VIDs (A, C, D1)        ~$1      DONE 2026-09-05
--       -- P4 QUERY B and QUERY D2 were not needed and can be skipped
--
--   6.1 build the 1% sample                        ~$1.20   DONE 2026-09-05
--   6.2 run 6.3 against the sample, then:          ~$3.60   GREEN 2026-09-05
--       -- (a)..(d) below are LITERAL MARKERS in this file: grep -n '^-- (' to find them.
--         (a) job stats, then per-stage            ~$0.01   run 1 FAILED, run 2 PASSES
--         (b) V1a/V1b/V1c                          ~$0.05   V1a FAILED -- the CHECK was
--                                                           wrong; fixed, and now passes by
--                                                           deduction from (d)
--         (c) V2a/V2b/V2c                          ~$0.60   V2b FAILED -- same cause, fixed
--         (d) diagnostic for a V2b failure         ~$0.60   DONE -- BENIGN, one fixup VID
--   6.2b materialize the exclusion side as tables  ~$1.30   DONE 2026-09-05  (5 tables)
--   6.2' re-run 6.2, then (a) again                ~$3.60   RUN 2 DONE -- PASSES
--       -- do NOT re-run 6.1; see below
--   6.3 BUILD the corrected table                  ~$119    BUILT 2026-09-05, job
--                                                           job_hv-YiMnw24rdq1X5lqkXVVyYP3qf
--                                                           -- unvalidated until 6.4
--   6.5 build `emptied_vids_t` FIRST               ~$0.90   DONE 2026-09-05. 172 VIDs
--       -- ORDER SWAPPED: 6.4 now depends on it              emptied where ZERO was expected.
--                                                           ADJUDICATED BENIGN: 0 SNPs, 172
--                                                           indels, all case-4 misses by the
--                                                           fixup workflows. 3,112,714 pairs.
--                                                           See the RESULT block before 6.6.
--   6.4 reconcile from the audit columns           ~$0.30   GREEN 2026-09-05. All 3 checks
--                                                           pass. Removed 2,341,606,369,
--                                                           inside the predicted band;
--                                                           excess 930 against VS-2010's
--                                                           ceiling of 932. See the RESULT
--                                                           block at the end of 6.4.
--   6.6 full validation                            ~$119    GREEN 2026-09-05. All 3 pass.
--                                                           THE TABLE IS FULLY VALIDATED.
--   6.6R reducer rewrite of 6.6, same price        ~$119    NOT NEEDED -- 6.6 took a
--                                                           32-stage plan and co-located on
--                                                           vid. Keep as a fallback.
--   6.7 extract the audit columns to a side table  ~$0.30   <- NEXT. MUST precede 6.8, which
--                                                           makes the columns unreadable.
--   6.8 DROP the audit columns -- promotes it      free     HELD INDEFINITELY 2026-09-05 at
--       -- see the HOLD note in the 6.8 section                the user's request. NOT before
--                                                             STEP 8 either -- STEP 8 uses
--                                                             EXPORT DATA ... AS SELECT so
--                                                             the schema need not change.
--                                                             Run only once collaborator
--                                                             questions have died down.
--   6.9 drop the rehearsal and 6.2b scratch tables free     KEEP `emptied_vids_t`.
--
-- 6.2 IS THREE SEPARATE THINGS AND (a) IS THE ONE WITH TEETH. (a), (b) and (c) are not
-- section headings anywhere below -- they are the three queries STEP 6.2 asks you to run in
-- order, and the line numbers above locate them. (b) and (c) are named for the `check_name`
-- strings their rows carry: V1a/V1b/V1c and V2a/V2b/V2c. Building the sample output
-- proves nothing by itself -- correctness was never in doubt, the membership filter being a
-- total function of its two inputs. The rehearsal exists to measure EXECUTION SHAPE, because
-- that is what failed for the scrub: its analogous statement held all 2,000 on-demand slots
-- for 82 minutes and never finished. So run the job-stats query BEFORE the V-checks; if the
-- shuffle bytes are wrong, the V-checks will pass and the full-scale run will still fail.
--
-- THAT IS EXACTLY WHAT HAPPENED ON RUN 1, WHICH IS WHY 6.2b EXISTS. The sample output built
-- cleanly in 44 seconds with zero spill, and the per-stage breakdown still says do not spend
-- the $119: four stages moved the ENTIRE sampled table's `person_ids` (215.7 GiB each at
-- ~14,462 bytes per output record), where the three-branch design was supposed to move only
-- the 0.34% affected rows. The full result and the reasoning are in STEP 6.2 under "RESULT,
-- RUN 1". The V-checks would have passed. Do not run them until (a) does.
--
-- 6.2b IS THE REMEDY AND IT WAS PRE-WRITTEN. STEP 6.3's shape notes already said that if the
-- probe side got repartitioned anyway, the fix was to materialize the exclusion side so the
-- planner has real statistics. 6.2b does that; 6.3 now reads the resulting `_t` tables
-- instead of CTEs. Nothing about the predicate or the branch logic changed.
--
-- AFTER 6.2b, RE-RUN 6.2 BUT NOT 6.1. The OUTPUT table must be rebuilt by the new 6.3 --
-- re-reading run 1's job id shows run 1's stages forever, so the job-stats query is only
-- measuring the new plan if a new job produced it. The INPUT sample must NOT be rebuilt:
-- `TABLESAMPLE SYSTEM` would draw different blocks, and then run 1 and run 2 would differ in
-- both the plan and the data, which is precisely the comparison being made. Keeping
-- `vid_to_participant_mapping_2026_08_28_sample` fixed makes run 2 a controlled experiment
-- and saves $1.20.
--
--
-- TABLE NAMES
--
--   READ   `foxtrot.vid_to_participant_mapping_2026_08_28`         the delivered table, untouched
--   WRITE  `foxtrot.vid_to_participant_mapping_2026_09_05`         the corrected table
--   TEMP   `foxtrot.vid_to_participant_mapping_2026_08_28_sample`  1% rehearsal input
--   TEMP   `foxtrot.vid_to_participant_mapping_2026_09_05_sample`  1% rehearsal output
--   WRITE  `foxtrot.mapping_correction_audit`                      per-VID impact, kept
--
-- The rehearsal tables are the production names plus a `_sample` suffix, so a stray one is
-- self-describing in a table listing and sorts next to the table it was drawn from. Note the
-- consequence for STEP 6.2's sed, which appends rather than substitutes and is therefore not
-- idempotent -- see the warning there. STEP 6.9 drops both.
--
-- The output follows the dated-version convention the dataset already uses, so there is no
-- rename and no in-place swap: the delivered table stays exactly as it is and becomes its
-- own rollback target. `2026_09_05` should be the date the build actually runs. Every
-- occurrence of that string names the output table and nothing else, so if the build slips
-- a day, `sed -i '' 's/2026_09_05/<new date>/g'` on this file is the entire change.
--
-- Restructured from the plan's original sketch per QUERY C of `step5_d1_reconcile.sql`.
-- The sketch ran D2 first (per-VID `delivered / excluded / corrected`, ~20 TB) and Step 6
-- second (~20 TB), paying twice for the same `person_ids` column. This version computes
-- `delivered` and `n_removed` as columns of the output, so the impact report and the
-- `corrected < 0` gate are both read off the artifact that would actually ship, for
-- pennies, and the expensive column is read once.
--
--
-- COST MODEL
--
-- Everything expensive here is expensive for one reason: touching `person_ids`. The scrub
-- measured one full pass at ~19 TiB; call that T, about $119 at $6.25/TiB on demand.
-- Bytes are billed per column per STATEMENT, so a statement that references `person_ids`
-- twice -- as STEP 6.3's UNION ALL does -- still costs 1 x T.
--
--   P1  layout preflight        INFORMATION_SCHEMA                    -- free
--   P2  exclusion size profile  `exclusions` only                     -- ~$0.30
--   P3  vid reconciliation      vid columns only                      -- ~$0.40
--   P4  extreme-VID inspection  clustered backup, pruned              -- ~$3
--   6.1 build the 1% sample     0.01 x T                              -- ~$1.20
--   6.2 rehearse on the sample  ~0.03 x T (build + both validations)  -- ~$3.60
--   6.3 BUILD                   1 x T                                 -- ~$119
--   6.4 free validation         scalar columns only                   -- ~$0.30
--   6.5 cheap validation        vid columns only                      -- ~$0.40
--   6.6 full validation, opt.   1 x T (reads the corrected copy only) -- ~$119
--   6.7 audit side-table        scalar columns only                   -- ~$0.30
--   6.8 drop audit columns      DDL metadata                          -- free
--   6.9 housekeeping            DROP TABLE                            -- free
--
-- So the floor is ~$128 and the recommended path with 6.6 is ~$247. Step 5 cost $509, so
-- 6.6 is 23% on top of what has already been spent to establish the predicate, on a table
-- that is about to be delivered to an external collaborator for the third time. Run it.
--
--
-- TWO DESIGN DECISIONS, BOTH LOAD-BEARING
--
-- 1. FILTER ARRAY MEMBERSHIP; NEVER SUBTRACT COUNTS.
--
-- QUERY B of `step5_d1_reconcile.sql` found 932 duplicate `(vid, person_id)` pairs in
-- `exclusions` across 109 people -- the genome-wide extent of VS-2010's duplicated ingest
-- batches. Membership filtering is idempotent under those duplicates; subtracting
-- `COUNT(*)` per VID is not, and would over-remove 932 entries with no way to notice.
--
-- The same property makes the build robust to two other things the delivered table carries
-- no constraint against: a repeated `vid` row (each row is filtered independently, where a
-- `GROUP BY vid` would silently merge them) and a repeated `person_id` within an array
-- (every copy of an excluded person is removed, and copies of a retained person are left
-- alone -- see the note under STEP 6.3).
--
-- 2. THE 5,776 FIXUP VIDS ARE EXCLUDED FROM THIS STEP AND LEFT TO STEP 6a.
--
-- Rows written by `GvsMapUnmappedVIDs` and `GvsMapDroppedDuplicateVIDs` name carriers found
-- at a DIFFERENT allele representation than the VID's own coordinates. `exclusions` was
-- computed at the VID's own coordinates, via `allele_status`, so for these VIDs the two are
-- not talking about the same genotypes.
--
-- THE TWO SOURCE TABLES ARE NOT DISJOINT, which matters for how the next two paragraphs
-- are read. P3 measured 5,776 distinct VIDs across them on 2026-09-05, against 5,227
-- unmapped and 1,542 dropped-duplicate individually -- so 993 VIDs are in BOTH, two thirds
-- of the dropped-duplicate table. For those 993 the unsafe reading below governs. Nothing
-- about the carve-out changes, since it is a `UNION DISTINCT` and always absorbed the
-- overlap; what changes is that "the unmapped VIDs" and "the dropped-duplicate VIDs" are
-- not two populations with opposite properties, they are two labels a VID can both carry.
--
-- For a VID that is ONLY unmapped this is harmless either way: the VID's own coordinates
-- match nothing in `alt_allele` -- that is the definition of unmapped -- so Step 5 emitted
-- no exclusion rows for it and there is nothing to apply. P3 confirms this exactly: of the
-- 121 fixup VIDs present in `exclusions`, ZERO are unmapped-sourced.
--
-- For a VID carrying the dropped-duplicate label it is NOT harmless, and this is the reason
-- for the boundary. Those VIDs did match at their own coordinates -- that is why the base join
-- wrote a row before `GvsMapDroppedDuplicateVIDs` deleted and re-inserted it -- so Step 5
-- has exclusion rows for them, computed against the left-aligned representation, while the
-- delivered array holds carriers drawn from the `input_*` representations. Applying one to
-- the other can remove a person whose delivered membership was justified by a genotype at a
-- representation that passes `FT`. That is over-removal of real carriers, the one failure
-- mode this whole exercise exists to avoid, and it fails silently.
--
-- 5,776 VIDs out of 1,601,242,198. They pass through this step untouched and are corrected
-- exactly in Step 6a, which evaluates `FT` per `(vid, input_location, input_ref,
-- input_alt)` and unions the survivors. **Step 6a must run before Step 8 re-exports.**
--
--
-- WHAT IS DELIBERATELY NOT DONE HERE
--
-- Duplicate `person_id` entries within a delivered array are left in place unless the
-- person is excluded. Roughly 1.1 million such entries exist genome-wide, out of ~2.56
-- trillion -- the VS-2010 duplicated vet rows whose alleles reached the VAT and passed
-- `FT`. Deduplicating them is a one-word change (`ARRAY(SELECT DISTINCT p ...)`) and is
-- tempting while the table is open, but it would touch rows this correction otherwise
-- leaves alone, and it would falsify the statement that the output is exactly the delivered
-- table minus the excluded members -- which is what STEP 6.2's multiset comparison and
-- STEP 6.4's reconciliation both test. Left as a separate decision; VS-2011 fixes it at
-- source for the next callset.
--
-- =====================================================================================


-- =====================================================================================
-- P1 -- what the delivered table physically looks like right now. Free.
-- =====================================================================================
--
-- Worth reading before writing a replacement, because the table in the dataset today is
-- NOT what `GvsCreateParticipantMappingTable.wdl` produced. That workflow's CTAS emits a
-- table ordered by `(location, ref, alt)` and then applies clustering with a separate
-- `bq update --clustering_fields=vid` (line 90), which sets the clustering spec as metadata
-- and does not recluster the data already written. The scrub then replaced the whole thing
-- with a plain CTAS that has no `ORDER BY`, no `PARTITION BY` and no `CLUSTER BY`.
--
-- So the expectation is: unpartitioned, and either unclustered or carrying a clustering
-- spec that the data does not honor. Two consequences if that is what comes back.
--
-- The open question in the plan -- "does the delivered Parquet need to be ordered?" -- has
-- already been answered in practice rather than in principle. The scrub's output was
-- unordered and it has been delivered. Unless the collaborator has reported a problem with
-- that second delivery, ordering is not required and STEP 6.3 should not pay for it. A
-- global `ORDER BY` over 1.6 billion rows is a single-reducer sort of the entire
-- `person_ids` column and is the most likely thing in this file to fail outright.
--
-- And `CLUSTER BY vid` on the new table is a free improvement rather than a change: it
-- gives the delivered table the physical property the WDL only ever asserted.

SELECT
  table_name,
  ddl
FROM `foxtrot.INFORMATION_SCHEMA.TABLES`
WHERE table_name LIKE "vid_to_participant_mapping%"
   OR table_name IN ("exclusions", "allele_status");

-- `TABLE_STORAGE` is a region-scoped view, not a dataset-scoped one, so it is qualified
-- differently from every other INFORMATION_SCHEMA reference in this file and filtered on
-- `table_schema`. Adjust `region-us` if the dataset lives elsewhere.
SELECT
  table_name,
  total_rows,
  ROUND(total_logical_bytes / POW(1024, 4), 2)        AS logical_tib,
  ROUND(total_logical_bytes / POW(1024, 3) * 0.02, 2) AS usd_per_month_active_storage
FROM `region-us`.INFORMATION_SCHEMA.TABLE_STORAGE
WHERE table_schema = "foxtrot"
  AND (table_name LIKE "vid_to_participant_mapping%" OR table_name IN ("exclusions", "allele_status"))
ORDER BY total_logical_bytes DESC;

-- The `vid_to_participant_mapping%` wildcard is there to ENUMERATE THE VERSIONS, so that
-- the one being read is known to be the current one.
--
--
-- RESULT, RUN 2026-09-05
--
--   table                                              rows           TiB    $/month
--   vid_to_participant_mapping                         1,601,270,624  18.82  385.36   clustered
--   vid_to_participant_mapping_2026_06_01_prescrub...  1,601,242,198  18.82  385.36   clustered
--   vid_to_participant_mapping_2026_08_28              1,601,242,198  18.65  381.90   neither
--   vid_to_participant_mapping_2025_10_28                 34,158,793   0.40    8.18   clustered
--   vid_to_participant_mapping_2026_06_01_sample          15,655,638   0.19    3.85
--   vid_to_participant_mapping_2026_06_01_sample_s...     15,655,638   0.19    3.81
--   four zero-row `_2026_06_01*` tables                            0   0.00    0.00
--
--   57.07 TiB, ~$1,168/month, before `2026_09_05` adds another ~$382.
--
-- THE P1 PREDICTION HELD. `2026_08_28` is unpartitioned AND unclustered, which is the shape
-- of the scrub's plain CTAS and not of anything `GvsCreateParticipantMappingTable.wdl`
-- produces. So `CLUSTER BY vid` in STEP 6.3 restores a property the delivered table lost.
--
-- THE LINEAGE IS CONFIRMED BY THE DELTA, not just by the names. `_prescrub_backup_vs_2000`
-- and `2026_08_28` differ by 0.17 TiB, which is 0.90% of 18.82 -- against the 0.903% of
-- entries the withdrawn/control scrub removed.
--
-- AND THE COST MODEL IS CONFIRMED. `person_ids` is essentially all of the 18.65 TiB (1.6
-- billion `vid` strings are ~0.03 TiB), and 2.56 trillion entries x 8 bytes is 18.6 TiB. So
-- T is 18.65 TiB and $116.56, marginally under the $119 used above; no line item moves.
--
-- TWO OTHER TABLES LOOK RELEVANT AND ARE NOT.
--
-- `vid_to_participant_mapping`, unsuffixed and clustered, has 28,426 MORE rows than either
-- delivery. It is the ORIGINAL mapping table, built when fewer samples had been withdrawn,
-- and its canonical-looking name is an accident of being first. Not a version of anything
-- current; do not read from it and do not reconcile against it.
--
-- `vid_to_participant_mapping_2025_10_28`, 34 M rows, is likewise old and superseded.
--
-- STORAGE IS NOT A CONCERN HERE. The family totals 57.07 TiB and ~$1,168/month and this step
-- adds ~$382, which the collaborator has never been sensitive to. Nothing in this file drops
-- a production table; STEP 6.9 drops only the two rehearsal tables it created.
--
--
-- THE SCRUB'S ROW COUNTS ARE EQUAL FOR A REASON THAT MATTERS TO STEP 6.
--
-- `_prescrub_backup_vs_2000` and `2026_08_28` both hold 1,601,242,198 rows, which looks like
-- the scrub emptied no VIDs. It did empty two -- the case-4 representation mismatches
-- documented at the bottom of `scrub_participant_mapping_table.sql` -- and those two rows
-- were then ADDED BACK MANUALLY. The counts match because 2 out and 2 in.
--
-- So `2026_08_28` contains two rows whose `person_ids` were not produced by the base join.
-- That is the same hazard as the 5,776 fixup VIDs and it needs the same treatment: whatever
-- those arrays hold, `exclusions` was computed at the VID's own coordinates and does not
-- describe them, so applying it can empty the rows a second time. Two rows out of 1.6
-- billion, and precisely the two that someone already had to intervene on by hand.
--
-- The two rows are:
--
--   17-72146919-GAGAAGAAGAAGAAGAAGAAGAAGAAGAAGA-G  ->  [<PERSON_A>]   30bp deletion in a GAA repeat
--   1-143186828-T-TG                               ->  [<PERSON_B>]   1bp insertion
--
-- BUT THE CARVE-OUT TURNS OUT TO BE UNNECESSARY, and the reason is worth stating because it
-- is an argument rather than a measurement. Step 5 built `exclusions` under
-- `si.withdrawn IS NULL AND si.is_control = false` (`step5_build_exclusions.sql:265-268`),
-- so it only ever names active non-control people. The scrub emptied these two VIDs because
-- EVERY person in their delivered arrays was withdrawn or control -- that is what emptying
-- means here. So the active non-control base join yields nothing at these coordinates, and
-- `exclusions` therefore contains no rows for either VID. They land in STEP 6.3's
-- `unaffected` branch and are copied through untouched, with or without a carve-out.
--
-- VERIFIED, 2026-09-05. All three checks below were run and all three passed: (a) returned
-- zero rows, (b) returned both people as active and non-control, (c) returned n = 1 for each
-- VID. So these two rows are inert to STEP 6.3, the manual fix used valid participants rather
-- than restoring what the scrub had removed, and the INSERTs did not duplicate a vid. The
-- queries are retained because the argument above is only as good as its premises.
--
--   -- (a) neither VID is in `exclusions`. Expect ZERO ROWS.
--   SELECT e.vid, COUNT(*) AS exclusion_rows,
--          COUNTIF(e.person_id IN (<PERSON_A>, <PERSON_B>)) AS names_the_manual_person
--   FROM `foxtrot.exclusions` AS e
--   WHERE e.vid IN ("17-72146919-GAGAAGAAGAAGAAGAAGAAGAAGAAGAAGA-G", "1-143186828-T-TG")
--   GROUP BY e.vid;
--
--   -- (b) the two manually inserted people are themselves active and non-control. Two rows,
--   -- `withdrawn` NULL and `is_control` false. Anything else means the manual fix put back
--   -- what the scrub had just taken out, which is a finding in its own right -- and one
--   -- Step 6 will NOT correct, since (a) says these VIDs are outside its reach.
--   SELECT si.sample_name, si.withdrawn, si.is_control
--   FROM `foxtrot.sample_info` AS si
--   WHERE SAFE_CAST(si.sample_name AS INT64) IN (<PERSON_A>, <PERSON_B>);
--
--   -- (c) the INSERTs added one row each rather than duplicating an existing vid. ~$0.19,
--   -- because the delivered table is unclustered so a vid filter still scans the vid column.
--   SELECT vid, COUNT(*) AS n
--   FROM `foxtrot.vid_to_participant_mapping_2026_08_28`
--   WHERE vid IN ("17-72146919-GAGAAGAAGAAGAAGAAGAAGAAGAAGAAGA-G", "1-143186828-T-TG")
--   GROUP BY vid;
--
-- The `-- MANUAL:` marker in each `fixup_vids` CTE (STEP 6.3, STEP 6.4, STEP 6.6) holds the
-- two VIDs ready to uncomment. Leave them commented -- check (a) says they are unreachable.
-- They are kept only so that a future re-run against a rebuilt `exclusions` has the carve-out
-- to hand rather than having to rediscover the two VIDs. If they are not carved out, STEP 6.5's emptied-VID list is the backstop
--     and will name them -- but that is after paying $117, and the whole point of STEP 6.5
--     is that a nonzero count there should be news.


-- =====================================================================================
-- P2 -- profile the exclusion sets. Reads `exclusions` only, ~$0.30.
-- =====================================================================================
--
-- This sizes the one part of STEP 6.3 whose execution shape is genuinely uncertain: the
-- per-element membership test `p NOT IN UNNEST(x.ids)`. Its cost per affected row is
-- (delivered array length) x (exclusion set size) if BigQuery scans the array linearly, and
-- (delivered array length) if it hashes it. The scrub established that the linear case is
-- real -- `p IN UNNEST(g.ids)` against a 450,000-element array took 12+ minutes on 1% of
-- the table -- so the number that matters is `max_excluded_per_vid`.
--
-- The mean is already known to be 431.2, which is nothing. A p999 in the low thousands is
-- also nothing. A max in the hundreds of thousands, on a VID that also has a very long
-- delivered array, is the case that stalls, and STEP 6.2's rehearsal is where it surfaces.
--
-- `distinct_person_ids` versus `rows` restates QUERY B's 932 surplus rows at whole-table
-- scale; it is here so the STEP 6.4 reconciliation has its expected excess in hand.

SELECT
  COUNT(*)                                                       AS exclusion_rows,
  COUNT(DISTINCT vid)                                            AS affected_vids,
  COUNT(DISTINCT FORMAT("%s|%d", vid, person_id))                AS distinct_pairs,
  COUNT(*) - COUNT(DISTINCT FORMAT("%s|%d", vid, person_id))     AS surplus_rows
FROM `foxtrot.exclusions`;

WITH per_vid AS (
  SELECT vid, COUNT(DISTINCT person_id) AS excluded
  FROM `foxtrot.exclusions`
  GROUP BY vid
)
SELECT
  COUNT(*)                                              AS affected_vids,
  ROUND(AVG(excluded), 1)                               AS mean_excluded,
  APPROX_QUANTILES(excluded, 1000)[OFFSET(500)]         AS p50,
  APPROX_QUANTILES(excluded, 1000)[OFFSET(990)]         AS p99,
  APPROX_QUANTILES(excluded, 1000)[OFFSET(999)]         AS p999,
  MAX(excluded)                                         AS max_excluded_per_vid,
  SUM(excluded)                                         AS total_distinct_pairs,
  -- Size of the array side of STEP 6.3's join, which decides broadcast versus shuffle.
  ROUND(SUM(excluded) * 8 / POW(1024, 3), 1)            AS ids_payload_gib,
  -- A LOWER bound on the membership test's total element comparisons, since the true cost
  -- is SUM(delivered x excluded) and delivered >= excluded for every VID. This is the
  -- number that condemned the single-branch array test.
  ROUND(SUM(POW(excluded, 2)) / 1e12, 2)                AS sum_excluded_sq_e12
FROM per_vid;

-- Picks the heavy/light boundary for STEP 6.3, by showing where the excluded entries
-- actually sit. Run it with the two above; it is the same scan.
WITH per_vid AS (
  SELECT vid, COUNT(DISTINCT person_id) AS excluded
  FROM `foxtrot.exclusions`
  GROUP BY vid
)
SELECT t, COUNTIF(excluded > t) AS vids_above, SUM(IF(excluded > t, excluded, 0)) AS entries_above
FROM per_vid, UNNEST([1000, 5000, 10000, 50000, 100000, 250000]) AS t
GROUP BY t ORDER BY t;

-- -------------------------------------------------------------------------------------
-- RESULT, RUN 2026-09-05. This is the section that changed STEP 6.3's design.
--
--   exclusion_rows 2,344,756,787   affected_vids 5,437,616
--   distinct_pairs 2,344,755,855   surplus_rows 932
--   mean 431.2   p50 3   p99 5,614   p999 87,664   max 535,659
--   ids_payload_gib 17.5   sum_excluded_sq_e12 247.27
--
--   t          vids_above   entries_above
--   1,000        152,924    2,190,144,230
--   5,000         58,786    1,978,552,482
--   10,000        36,422    1,820,821,825
--   50,000        10,053    1,250,705,584
--   100,000        4,702      872,812,375
--   250,000          686      255,834,268
--
-- TWO THINGS CAME OUT OF THIS.
--
-- The 932 surplus rows are the VS-2010 duplicates at whole-table scale; they are inert
-- here (`DISTINCT` in STEP 6.3's `ex_pairs`) and are the basis of the genome-wide
-- ~109-113 sample estimate now on that ticket.
--
-- The tail is the real finding, in two senses at once. As a performance matter it made the
-- single-branch array test untenable, which is why STEP 6.3 now splits at t = 1000 -- see
-- the reasoning there. As a DATA matter it is the headline of the whole investigation:
-- `max_excluded_per_vid` = 535,659 against 535,662 active non-control participants
-- (`SELECT COUNT(*) FROM foxtrot.sample_info WHERE is_control IS FALSE AND withdrawn IS
-- NULL`, 2026-09-05), and the top 20 all sit between 532,570 and 535,659. Those VIDs
-- currently name essentially every participant in the callset as a carrier and should name
-- a handful. 686 VIDs exceed 250,000.
--
-- TWO SEPARATE THINGS NEED EXPLAINING HERE, and they are easy to run together.
--
-- WHY IS NEARLY EVERYONE A CARRIER? Because `delivered >= excluded`, every one of the top
-- 20 has at least 99.4% of the cohort in its array. That is the class of site the team has
-- run into before -- hg38 carries what is actually the minor allele, so essentially every
-- sample is called non-reference against it. Nothing about the mapping table causes this
-- and nothing in Step 6 changes it; the array is a faithful record of who was called.
--
-- WHY DOES NEARLY EVERY ONE OF THOSE GENOTYPES GET EXCLUDED? This is the open part, and
-- the reference-error reading actually argues AGAINST it: a common reference-error site
-- should be in the truth resources, earn `POSITIVE_TRAIN_SITE` and hence `yng_status =
-- "Y"`, and pass `FT` outright through `LOGICAL_OR(is_yes)` whatever its calibration
-- sensitivity -- the term that rescues 64% of over-threshold Foxtrot alleles. So these
-- genotypes ought to pass, and they do not.
--
-- DO NOT REACH FOR THE SITE-LEVEL FILTER RULE HERE. CLAUDE.md's "a site is judged by its
-- best allele of a class" is real but belongs to `ExtractCohortEngine`, the VCF/PGEN
-- extract. The VDS takes its site filters from `filter_set_sites` instead
-- (`import_gvs.py:339`), and that table never holds calibration-sensitivity failures -- only
-- EXCESS_ALLELES, ExcessHet, LowQual and NO_HQ_GENOTYPES. So THERE IS NO VETS SITE-LEVEL
-- RULE IN THE VDS AT ALL, and a VID exists in the VAT for exactly one reason: `gvs_all_ac >
-- 0`, meaning at least one FT-passing carrier.
--
-- WHICH MAKES THE ARGUMENT TIGHTER, NOT WEAKER. `FT` is `~any_no & (any_yes | all_ok)`
-- folded over the genotype's called non-ref allele indices (`import_gvs.py:370-381`), so it
-- is a deterministic function of the ALLELE SET. Two carriers of A1 whose called alleles
-- match get identical `FT`. Had all 535,662 been `0/1` for A1 they would pass or fail
-- together -- 3 versus 535,659 is impossible. So the survivors and the excluded must differ
-- in one of exactly two ways, and there is no third:
--
--   (a) GQ 0, which is per-genotype and independent of the alleles; or
--   (b) a different allele set -- the excluded are `1/2`, carrying a partner allele that
--       fails while the ~3 survivors are plain `0/1`.
--
-- That is why QUERY A's breakdown is DECISIVE rather than suggestive: `is_gq0` and
-- `is_ft_fail` exhaust the possibilities. QUERY D -- the site's allele structure from
-- `filter_set_info` -- is worth running only if `ft_fail` dominates.
--
-- MECHANISM (b) REQUIRES A1 TO BE `G`, NOT `Y`, and that is where the reference-error
-- reading comes back in as a PRECONDITION rather than an objection. `any_yes` is a
-- whole-genotype short-circuit, not a per-allele exclusion -- `ExtractCohortEngine.java:939`
-- passes the genotype outright on a single `Y` among its own alleles, matching Hail's
-- `any_yes | all_ok`. So if A1 carried `Y`, every `1/2` carrier would pass however bad the
-- partner allele was, and (b) could not operate at all. A common allele that the truth
-- resource maintainers have not catalogued -- which is what an unrecognized hg38
-- reference-error site is -- is precisely a common allele scored `G`.
-- -------------------------------------------------------------------------------------


-- =====================================================================================
-- P3 -- reconcile the VID sets before touching any arrays. vid columns only, ~$0.40.
-- =====================================================================================
--
-- Three questions, one statement so the `vid` column of the mapping table is billed once.
-- Do NOT add `person_ids` to any of these -- referencing the column bills all 19 TiB no
-- matter how few rows come back.
--
-- `orphan_exclusion_vids` are exclusion rows whose VID has no row in the delivered table at
-- all. A handful is expected and benign: the scrub dropped VIDs whose every carrier was
-- withdrawn. A large number means `allele_status` and the delivered table disagree about
-- which VIDs exist, which must be understood before patching anything, because it would
-- mean the two were built from different VAT versions.
--
-- `fixup_vids_in_exclusions` is the number that justifies decision 2 above. The prediction
-- is that essentially all of them are dropped-duplicate VIDs and essentially none are
-- unmapped VIDs, for the reason given there. If unmapped VIDs show up in `exclusions` in
-- quantity, the premise that they match nothing at their own coordinates is wrong and
-- Step 6a needs rethinking before Step 6 runs.
--
-- `duplicate_vid_rows` re-runs the scrub's STEP 1 precondition, and as of the three-branch
-- rewrite it is a HARD GATE rather than a diagnostic. The membership filter is row-local
-- and handles duplicate vid rows correctly, but STEP 6.3's heavy branch reaggregates with
-- `GROUP BY vid` and would silently merge two rows sharing a `vid` into one. If this comes
-- back nonzero, do not run STEP 6.3 -- find out what produced the duplicates first.

WITH
fixup_vids AS (
  SELECT vid, "unmapped" AS source
  FROM `foxtrot.unmapped_vid_mapping_2026_06_01`
  GROUP BY vid
  UNION ALL
  SELECT vid, "dropped_duplicate"
  FROM `foxtrot.dropped_duplicate_mapping_2026_06_01`
  GROUP BY vid
),
exclusion_vids AS (
  SELECT DISTINCT vid FROM `foxtrot.exclusions`
),
mapping_vids AS (
  SELECT vid FROM `foxtrot.vid_to_participant_mapping_2026_08_28`
)
SELECT "mapping table rows"                       AS measure, CAST(COUNT(*) AS STRING) AS value FROM mapping_vids
UNION ALL
SELECT "distinct mapping vids",                   CAST(COUNT(DISTINCT vid) AS STRING) FROM mapping_vids
UNION ALL
SELECT "duplicate_vid_rows",
       CAST((SELECT COUNT(*) FROM (SELECT vid FROM mapping_vids GROUP BY vid HAVING COUNT(*) > 1)) AS STRING)
UNION ALL
SELECT "affected vids in exclusions",             CAST(COUNT(*) AS STRING) FROM exclusion_vids
UNION ALL
SELECT "orphan_exclusion_vids (not in mapping)",
       CAST((SELECT COUNT(*) FROM exclusion_vids e
             LEFT JOIN mapping_vids m USING (vid) WHERE m.vid IS NULL) AS STRING)
UNION ALL
SELECT "fixup vids, distinct",                    CAST((SELECT COUNT(DISTINCT vid) FROM fixup_vids) AS STRING)
UNION ALL
SELECT "fixup vids present in exclusions",
       CAST((SELECT COUNT(*) FROM (SELECT DISTINCT vid FROM fixup_vids) f
             JOIN exclusion_vids e USING (vid)) AS STRING)
UNION ALL
SELECT "  of those, unmapped-sourced",
       CAST((SELECT COUNT(DISTINCT f.vid) FROM fixup_vids f
             JOIN exclusion_vids e USING (vid) WHERE f.source = "unmapped") AS STRING)
UNION ALL
SELECT "  of those, dropped-duplicate-sourced",
       CAST((SELECT COUNT(DISTINCT f.vid) FROM fixup_vids f
             JOIN exclusion_vids e USING (vid) WHERE f.source = "dropped_duplicate") AS STRING);

-- -------------------------------------------------------------------------------------
-- RESULT, RUN 2026-09-05. Every gate passes; one recorded figure was wrong.
--
--   mapping table rows                    1,601,242,198
--   distinct mapping vids                 1,601,242,198
--   duplicate_vid_rows                                0   <- HARD GATE for STEP 6.3, passes
--   affected vids in exclusions               5,437,616   <- agrees with P2
--   orphan_exclusion_vids                             0
--   fixup vids, distinct                          5,776
--   fixup vids present in exclusions                121
--     of those, unmapped-sourced                      0
--     of those, dropped-duplicate-sourced           121
--
-- `duplicate_vid_rows` = 0 clears STEP 6.3's heavy branch to reaggregate with `GROUP BY
-- vid`. Nothing else in the file depends on it.
--
-- `orphan_exclusion_vids` = 0 is better than the "handful" this section predicted, and it
-- corroborates the hand-edit finding from P1: the scrub emptied exactly two VIDs and both
-- were re-inserted by hand, so no VID was left dropped for `exclusions` to point at. It
-- also simplifies STEP 6.4 -- the reconciliation's "present in the delivered table"
-- qualifier is now vacuous and can be read as covering every excluded pair.
--
-- 121 AND 0 IS DECISION 2's PREDICTION LANDING EXACTLY. The carve-out was justified by an
-- asymmetry -- unmapped VIDs match nothing at their own coordinates so Step 5 emits nothing
-- for them, while dropped-duplicate VIDs do match and so carry exclusions computed against
-- the wrong representation. Zero unmapped and all 121 dropped-duplicate is that asymmetry
-- measured rather than argued. It also quantifies the V2b rehearsal caveat: 1% of 121 is
-- ~1, so a V2b mismatch attributable to fixup VIDs should be zero or one VID, and anything
-- more is a different problem.
--
-- 5,776 IS NOT 6,769, and the earlier figure was never measured -- it was 5,227 + 1,542,
-- assuming the two source tables are disjoint. They are not: 993 VIDs are in both, two
-- thirds of the dropped-duplicate table. See the note on the `fixup_vids` CTE in STEP 6.3.
-- Nothing in Step 6 changes, but Step 6a's scope description does, and the reasoning that
-- treats the two tables as populations with opposite safety properties has to be read as
-- being about LABELS a VID may carry both of.
--
-- NEEDED FOR STEP 6.4's RECONCILIATION: how many excluded PAIRS belong to those 121 VIDs.
-- `SUM(n_removed)` is expected to equal distinct excluded pairs MINUS the fixup VIDs' pairs
-- PLUS the delivered arrays' duplicate multiplicity. RUN 2026-09-05, closing the last of the
-- three terms -- see the RESULT block below the query.
-- -------------------------------------------------------------------------------------

WITH fixup_vids AS (
  SELECT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`       GROUP BY vid
  UNION DISTINCT
  SELECT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01`  GROUP BY vid
)
SELECT
  COUNT(DISTINCT FORMAT("%s|%d", e.vid, e.person_id)) AS fixup_excluded_pairs,
  COUNT(DISTINCT e.vid)                               AS fixup_vids_with_pairs
FROM `foxtrot.exclusions` AS e
JOIN fixup_vids AS f ON f.vid = e.vid;

-- -------------------------------------------------------------------------------------
-- RESULT, 2026-09-05
-- -------------------------------------------------------------------------------------
--   fixup_excluded_pairs     37,702
--   fixup_vids_with_pairs       121
--
-- INTERNALLY CONSISTENT: `fixup_vids_with_pairs` = 121 matches P3's `fixup vids present in
-- exclusions` exactly, so every fixup VID that appears in `exclusions` at all contributes at
-- least one pair -- no VID is counted by one query and not the other.
--
-- THE HAZARD DID NOT MATERIALIZE. 37,702 over 121 VIDs is a mean of 311.6, BELOW the global
-- 431.2 (P2). The reason for measuring rather than assuming was that a dropped-duplicate VID
-- is exactly the kind that can hold a very long array, and all 121 are dropped-duplicate-
-- sourced -- but they turn out to be ordinary-weight. At 0.0016% of the 2,344,755,855
-- distinct pairs, the carve-out barely moves the reconciliation.
--
-- STEP 6.4's IDENTITY IS NOW CLOSED. Substituting all three measured terms:
--
--   SUM(n_removed)  ==  2,344,755,855 - 37,702 + [0 .. 932]
--                   ==  2,344,718,153 to 2,344,719,085
--
-- A 932-wide window on 2.34 billion. Anything outside it is a defect, not noise.


-- -------------------------------------------------------------------------------------
-- P3c -- CONFIRM THE 993 OVERLAP DIRECTLY. Two tiny tables, effectively free.
-- -------------------------------------------------------------------------------------
-- 993 is currently ARITHMETIC, not a measurement: `5,227 + 1,542 - 5,776`. That subtraction
-- is only valid if 5,227 and 1,542 were themselves DISTINCT-VID counts rather than row
-- counts, and the whole point of the 6,769 correction was that an unverified assumption
-- about these two tables had gone unchallenged for a long time. So measure the three
-- quantities in one statement instead of deriving the third from the other two. Costs
-- nothing; run it with P4's QUERY A.
--
-- EXPECT 5,227 / 1,542 / 993, and `unmapped + dropped_dup - in_both` = 5,776 to agree with
-- P3. If `unmapped_vids` or `dropped_dup_vids` comes back BELOW its published figure, that
-- figure was a row count and Step 6a's scope description needs revisiting again -- Step 6's
-- carve-out is still unaffected either way, since `UNION DISTINCT` never depended on any of
-- this.

WITH
u AS (SELECT DISTINCT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`),
d AS (SELECT DISTINCT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01`)
SELECT
  (SELECT COUNT(*) FROM u)                                       AS unmapped_vids,
  (SELECT COUNT(*) FROM d)                                       AS dropped_dup_vids,
  (SELECT COUNT(*) FROM u JOIN d USING (vid))                    AS in_both,
  (SELECT COUNT(*) FROM (SELECT vid FROM u UNION DISTINCT SELECT vid FROM d)) AS union_vids;

-- -------------------------------------------------------------------------------------
-- RESULT, 2026-09-05 -- CONFIRMED, all four figures as predicted.
-- -------------------------------------------------------------------------------------
--   unmapped_vids      5,227
--   dropped_dup_vids   1,542
--   in_both              993
--   union_vids         5,776
--
-- So 5,227 and 1,542 WERE distinct-VID counts, the `5,227 + 1,542 - 5,776` derivation was
-- sound, and 993 is now measured rather than inferred. `union_vids` agreeing with P3's
-- independently computed 5,776 closes the loop. Nothing downstream changes; this retires the
-- question rather than answering it differently.


-- =====================================================================================
-- P4 -- inspect the extreme VIDs. Three small statements, a few dollars all in.
-- =====================================================================================
--
-- P2 found VIDs excluding up to 535,659 of 535,662 possible participants. This section
-- looks at a handful of them directly. It is NOT a performance rehearsal -- STEP 6.3's
-- heavy branch is linear and these VIDs cost it nothing in particular.
--
-- WHAT IS ALREADY SETTLED, AND SO NOT WORTH A QUERY. `delivered >= excluded` holds for
-- every VID by construction, and `delivered <= 535,662`. For the worst VID that pins its
-- delivered array to [535,659, 535,662] -- a three-element range -- and for the 20th to
-- within 0.6%. So the inflation factors at the top of the distribution do not need
-- measuring; they follow from P2. They are also the signature the team has seen before, of
-- hg38 carrying a minor allele at the site: at least 99.4% of the cohort is named as a
-- carrier at every one of the top 20.
--
-- WHAT IS NOT SETTLED, AND IS WHY QUERY A NOW BREAKS DOWN BY CAUSE. A common
-- reference-error site should sit in the truth resources, earn `POSITIVE_TRAIN_SITE` and
-- therefore `yng_status = "Y"`, and pass `FT` through the `LOGICAL_OR(is_yes)` term
-- regardless of its calibration sensitivity -- that term rescues 64% of over-threshold
-- Foxtrot alleles. "Nearly everyone is non-reference here" thus predicts these genotypes
-- PASS. They do not. `exclusions` carries `is_gq0` and `is_ft_fail`, so the breakdown says
-- which, for the price of one more aggregate on a scan already being paid for.
--
-- THE ONE THING WORTH A FULL-SCALE CHECK is the argument STEP 6.3 relies on to claim
-- nothing gets dropped -- that every VAT VID has `gvs_all_ac > 0` and so retains at least
-- one carrier. At `delivered - excluded` of ~3 that argument has almost no margin left, and
-- if `gvs_all_ac` for these VIDs is single digits while the mapping table names half a
-- million people, the Step 5 predicate is confirmed at ~1e5x inflation -- a far stronger
-- validation than the median case can give. That is QUERY C, and it is the cheap one.
--
-- WHY THE PRESCRUB BACKUP AND NOT THE DELIVERED TABLE, for QUERY B. P1 found that
-- `vid_to_participant_mapping_prescrub_backup_vs_2000` is CLUSTERED BY vid while
-- `2026_08_28` is not, so a literal `WHERE vid IN (...)` against the backup gets real block
-- pruning and costs a few dollars where the same query against the delivered table costs a
-- full 1 x T. Its arrays are pre-scrub and so up to ~0.9% longer, which for this purpose is
-- fine -- it slightly overstates the delivered length, in the conservative direction.
--
-- A LITERAL IN-LIST IS REQUIRED for that pruning. A join to a subquery does not prune, so
-- run QUERY A, paste its VIDs into QUERY B, and do not be tempted to combine them.

-- QUERY A -- the 20 worst VIDs, BROKEN DOWN BY CAUSE. Reads `exclusions` only, ~$0.30.
--
-- The breakdown is the point, not the ranking. See the note under P2's RESULT block: these
-- VIDs are near-universal non-reference sites, and a common reference-error site ought to
-- sit in the truth resources, earn `POSITIVE_TRAIN_SITE`, and so pass `FT` through the
-- `LOGICAL_OR(is_yes)` term whatever its calibration sensitivity. That predicts these
-- genotypes PASS. They do not, so which column is doing the excluding is the question.
--
--   `gq0` dominant        -> a caller-confidence story. The site is non-reference for
--                            everyone but DRAGEN will not commit to the call, which points
--                            at repetitive or segdup context rather than at the filter
--                            model. QUERY D is then not worth running.
--   `ft_fail` dominant    -> the filter model is failing genotypes at a site it let
--                            through, which is the multi-allelic partner mechanism. Run
--                            QUERY D to see the site's allele structure.
--   both, per VID         -> expected to be rare; `is_gq0` and `is_ft_fail` are not
--                            mutually exclusive per row, so count them separately rather
--                            than assuming they partition.
SELECT
  vid,
  COUNT(DISTINCT person_id)                                    AS excluded,
  COUNT(DISTINCT IF(is_gq0,     person_id, NULL))              AS gq0,
  COUNT(DISTINCT IF(is_ft_fail, person_id, NULL))              AS ft_fail,
  COUNT(DISTINCT IF(is_gq0 AND is_ft_fail, person_id, NULL))   AS both
FROM `foxtrot.exclusions`
GROUP BY vid
ORDER BY excluded DESC
LIMIT 20;

-- -------------------------------------------------------------------------------------
-- RESULT, 2026-09-05 -- `ft_fail` IS THE WHOLE STORY. GQ 0 IS NOISE.
-- -------------------------------------------------------------------------------------
--   vid                       excluded    gq0   ft_fail   both
--   7-152451810-T-G            535,659      0   535,659      0
--   10-70272137-T-G            535,627      0   535,627      0
--   4-40960447-G-T             535,521      0   535,521      0
--   16-9293170-T-TCAGGTTTTTCA  535,508      0   535,508      0
--   2-61866200-C-T             535,108      1   535,108      1
--   14-89409352-C-G            534,966      2   534,966      2
--   3-65490858-G-A             534,941      0   534,941      0
--   12-4195211-C-G             534,882      0   534,882      0
--   3-94641999-A-G             534,877    163   534,877    163
--   12-28707196-C-A            534,758      2   534,758      2
--   10-38959806-A-AAAAG        534,745      0   534,745      0
--   13-37650636-T-C            534,703     20   534,703     20
--   16-29292938-T-C            534,591      0   534,591      0
--   2-55184106-A-T             534,512      0   534,512      0
--   19-6994372-A-T             534,292      1   534,292      1
--   7-65019556-G-A             534,253      0   534,253      0
--   1-27036512-G-A             534,245      0   534,245      0
--   16-74491696-C-G            533,451      0   533,451      0
--   2-208094726-G-C            532,726      0   532,726      0
--   13-80792215-G-T            532,570     12   532,570     12
--
-- TWO EXACT IDENTITIES HOLD IN EVERY ROW: `ft_fail = excluded`, and `both = gq0`. The first
-- says the filter model accounts for 100% of these exclusions. The second says GQ 0 is a
-- strict SUBSET of the FT failures -- there is not one person anywhere in the top 20 excluded
-- for being a no-call who would otherwise have passed. So GQ 0 contributes nothing
-- independent here and the caller-confidence reading is dead: it is the filter model, alone.
--
-- ============================ THE MECHANISM, CORRECTED ============================
--
-- This INVERTS the orientation written down before the query ran. The old reading had the
-- excluded people as `1/2` carriers dragged down by a failing partner, with the VID's own
-- allele passing. RULE THAT OUT ON POPULATION GENETICS, not on the filter model: it needs
-- 535,659 of 535,662 people to be het-var at one site, which requires two alt alleles both
-- near fixation AND near-universal heterozygosity. Hardy-Weinberg caps heterozygosity at 50%
-- for a 50/50 site. 99.999% is impossible.
--
-- SCOPE THAT ARGUMENT NARROWLY -- it holds only where the site's ARITY IS LOW, and QUERY D1
-- later found arity running from 2 to 39 across these twenty. With k alleles at equal
-- frequency heterozygosity reaches 1 - 1/k, so at a 39-allele site compound heterozygosity
-- can be the norm and HW forbids nothing. It is decisive at 7-152451810-T-G, where D1
-- measured `alleles_at_site` = 2. The conclusion below survives anyway, because D1 confirmed
-- it DIRECTLY at nineteen of twenty -- but one VID does turn out to run on the original
-- mechanism, and HW is not what separates them. See D1's RESULT block.
--
-- SO THE VID'S OWN ALLELE IS WHAT FAILS, for essentially everyone. And that closes the
-- question of who the survivors are, because `FT` is `~any_no & (any_yes | all_ok)`: if the
-- VID's allele is over threshold then `all_ok` is FALSE for every genotype containing it,
-- whatever else that genotype carries. The only surviving term is `any_yes`. The VID's own
-- allele is not `Y` -- a `Y` there would pass all 535,662 -- so the `Y` must sit on a PARTNER
-- allele, and only a multi-allelic carrier has one.
--
--   THE TOP-20 VIDs ARE IN THE VAT ONLY BECAUSE A HANDFUL OF `1/2` CARRIERS WERE RESCUED
--   BY THE `any_yes` WHOLE-GENOTYPE SHORT-CIRCUIT ON A PARTNER ALLELE.
--
-- Which makes the reference-error reading a matched pair rather than a single claim: at the
-- same site sit an UNCATALOGUED common allele, scored `G` and failing calibration sensitivity
-- because the model sees a site where nearly every sample is non-reference, and a CATALOGUED
-- one scored `Y`. That is exactly the shape of an hg38 position where the reference carries
-- the minor allele and the truth resources have only partly caught up. `any_yes` being a
-- whole-genotype short-circuit rather than a per-allele rescue -- the thing this file had
-- backwards until `ExtractCohortEngine.java:939` was read properly -- is load-bearing.
--
-- QUERY C IS NOW A SAFETY GATE, NOT A CONFIRMATION. The reading above predicts `gvs_all_ac`
-- in the single digits to low tens for these VIDs. If it comes back in the hundreds of
-- thousands, the Step 5 `FT` reconstruction is WRONG at these sites and the patch would gut
-- 535,000 legitimate entries per VID. The chr21:20-22Mb VDS validation reached 100.0000%
-- agreement but had no reason to contain a site of this kind, so the top of the distribution
-- is genuinely unvalidated. RUN QUERY C BEFORE SPENDING THE $119 ON STEP 6.3.

-- QUERY B -- OPTIONAL, and lower value than it looks. Their delivered array lengths, which
-- P2 already pins to within three for the worst VID and 0.6% for the 20th. Run it only if
-- an exact denominator is wanted for the write-up, or if QUERY C comes back surprising and
-- the bound needs confirming rather than inferring.
--
-- PASTE QUERY A's VIDs INTO THE IN-LIST. Check the bytes-billed estimate before running: it
-- should be single-digit GB. If the editor says terabytes, the pruning did not happen --
-- confirm the table name and that the list is literal rather than a subquery.
SELECT vid, ARRAY_LENGTH(person_ids) AS delivered_prescrub
FROM `foxtrot.vid_to_participant_mapping_prescrub_backup_vs_2000`
WHERE vid IN (
  -- PASTE HERE
  "7-152451810-T-G"
)
ORDER BY delivered_prescrub DESC;

-- QUERY C -- what the VAT says the allele count actually is. RUN THIS BEFORE STEP 6.3.
--
-- Promoted from confirmation to SAFETY GATE by QUERY A's result -- see the note there.
-- `gvs_all_ac` is the target the correction aims at, so `delivered - excluded` should land
-- near it. Order-of-magnitude agreement is the pass condition; an exact match is not
-- expected, since AC counts alleles and the mapping table counts people, and a homozygous
-- carrier contributes 2 to one and 1 to the other.
--
--   single digits to low tens  -> PASS, and a ~1e5x inflation confirmation far sharper than
--                                 the median case can give. Proceed to STEP 6.3.
--   hundreds of thousands      -> STOP. The Step 5 `FT` reconstruction is wrong at these
--                                 sites and the patch would remove ~535,000 legitimate
--                                 entries per VID. Do not spend the $119 until this is
--                                 understood.
--   VID absent from the VAT    -> it should not be in the mapping table either. A finding in
--                                 its own right, and one STEP 6.5 would otherwise surface
--                                 only after the money was spent.
-- COLLAPSED TO ONE ROW PER VID. The VAT is keyed `(vid, transcript)`, not `vid`, so a VID
-- overlapping a dozen transcripts returns a dozen identical copies of these three columns --
-- `gvs_all_ac`, `gvs_all_an` and `gvs_all_af` are VID-level facts stored redundantly on every
-- transcript row. `GROUP BY vid` is the whole fix.
--
-- `distinct_ac` EARNS ITS PLACE: collapsing with `ANY_VALUE` silently picks a winner if the
-- copies ever disagree, so count them instead of trusting them. Expect 1 in every row.
-- 0 means `gvs_all_ac` is NULL throughout for that VID; anything above 1 means the VAT
-- carries conflicting allele counts for one VID, which would be a finding of its own and
-- would make `ANY_VALUE` the wrong reduction. `vat_rows` is free alongside it and gives the
-- transcript multiplicity, which is worth an eyeball -- these are near-universal
-- non-reference sites, so a 0 there (no transcript overlap at all) is unremarkable.
SELECT
  vid,
  ANY_VALUE(gvs_all_ac)      AS gvs_all_ac,
  ANY_VALUE(gvs_all_an)      AS gvs_all_an,
  ANY_VALUE(gvs_all_af)      AS gvs_all_af,
  COUNT(*)                   AS vat_rows,
  COUNT(DISTINCT gvs_all_ac) AS distinct_ac
FROM `foxtrot.foxtrot_v4_2025_07_29_vat_v9_r2_p2`
WHERE vid IN (
  "7-152451810-T-G",
  "10-70272137-T-G",
  "4-40960447-G-T",
  "16-9293170-T-TCAGGTTTTTCA",
  "2-61866200-C-T",
  "14-89409352-C-G",
  "3-65490858-G-A",
  "12-4195211-C-G",
  "3-94641999-A-G",
  "12-28707196-C-A",
  "10-38959806-A-AAAAG",
  "13-37650636-T-C",
  "16-29292938-T-C",
  "2-55184106-A-T",
  "19-6994372-A-T",
  "7-65019556-G-A",
  "1-27036512-G-A",
  "16-74491696-C-G",
  "2-208094726-G-C",
  "13-80792215-G-T"
)
GROUP BY vid
ORDER BY gvs_all_ac;

-- -------------------------------------------------------------------------------------
-- RESULT, 2026-09-05 -- GATE PASSES. Inflation 444x to 534,942x.
-- -------------------------------------------------------------------------------------
-- `excluded` is carried over from QUERY A; `AN/2` and the last two columns are derived.
--
--   vid                          excluded     AC    AN/2   exc+AN/2  short   inflation
--   7-152451810-T-G               535,659      2       3    535,662      0    267,830x
--   10-70272137-T-G               535,627      4      21    535,648     14    133,908x
--   4-40960447-G-T                535,521     16      56    535,577     85     33,471x
--   16-9293170-T-TCAGGTTTTTCA     535,508      8     140    535,648     14     66,940x
--   2-61866200-C-T                535,108      9     507    535,615     47     59,457x
--   14-89409352-C-G               534,966    496     533    535,499    163      1,080x
--   3-65490858-G-A                534,941      1     432    535,373    289    534,942x
--   12-4195211-C-G                534,882     26     779    535,661      1     20,573x
--   3-94641999-A-G                534,877     63      85    534,962    700      8,491x
--   12-28707196-C-A               534,758     61      67    534,825    837      8,768x
--   10-38959806-A-AAAAG           534,745    835     829    535,574     88        641x
--   13-37650636-T-C               534,703      2     669    535,372    290    267,352x
--   16-29292938-T-C               534,591    205   1,043    535,634     28      2,609x
--   2-55184106-A-T                534,512    228   1,011    535,523    139      2,345x
--   19-6994372-A-T                534,292  1,205   1,285    535,577     85        444x
--   7-65019556-G-A                534,253    976   1,185    535,438    224        548x
--   1-27036512-G-A                534,245    299     848    535,093    569      1,788x
--   16-74491696-C-G               533,451      2   2,208    535,659      3    266,726x
--   2-208094726-G-C               532,726      5      49    532,775  2,887    106,546x
--   13-80792215-G-T               532,570  1,040   1,283    533,853  1,809        513x
--
-- MECHANICS FIRST. `distinct_ac` = 1 in all 20 rows, so `ANY_VALUE` was a safe reduction and
-- the VAT does not carry conflicting allele counts for a VID. `vat_rows` ran 1 to 14, so the
-- collapse was doing real work -- 14-89409352-C-G was returning fourteen identical copies.
--
-- THE GATE: PASSES, WITH ENORMOUS MARGIN. The failure condition was `gvs_all_ac` in the
-- hundreds of thousands, which would have meant Step 5's `FT` reconstruction was wrong here
-- and the patch was about to delete half a million legitimate entries per VID. The largest
-- AC in the set is 1,205. STEP 6.3 IS CLEARED TO RUN on this evidence.
--
-- ======================= `gvs_all_an` CLOSES THE COHORT, AND THAT IS =======================
-- ======================= A STRONGER RESULT THAN THE GATE ITSELF     =======================
--
-- `AN/2` is the number of participants holding a passing called genotype at the site. Add it
-- to `excluded` and EVERY ROW LANDS AT OR JUST BELOW THE 535,662 ACTIVE NON-CONTROL COHORT:
-- exact for 7-152451810-T-G, one short for 12-4195211-C-G, three short for 16-74491696-C-G,
-- median shortfall 139, worst 2,887 (0.54%). The shortfalls are participants with no call at
-- all at that site, which is ordinary.
--
-- THE DIRECTION OF THE ERROR IS THE POINT. If Step 5 were OVER-excluding -- if any person it
-- put in `exclusions` actually holds a passing genotype -- that person would be counted on
-- both sides and the sum would EXCEED 535,662. It never does. Not once, in twenty sites
-- chosen precisely because they are the most extreme in the callset. Over-exclusion is ruled
-- out ARITHMETICALLY here, not by sampling, and this is the first check in the whole
-- investigation that constrains the tail rather than the median. The chr21 VDS validation
-- established fidelity where the data is ordinary; this establishes it where it is not.
--
-- 7-152451810-T-G IS THE PUREST CASE IN THE CALLSET. `AN` = 6: exactly THREE participants of
-- 535,662 hold a passing call, and `excluded + 3` is the cohort to the person. Nobody is
-- homozygous reference. The entire cohort is non-reference at that position -- the hg38
-- reference-error signature with nothing left over, and 267,830x inflation in the delivered
-- mapping table.
--
-- 3-65490858-G-A IS THE WORST INFLATION AND THE SHARPEST REMAINING TEST. `AC` = 1: ONE alt
-- allele in the entire callset, so exactly one heterozygous participant, against 534,941
-- named as carriers. 534,942x. It is also the tightest case for STEP 6.5's no-empty argument
-- -- if that single person is in the delivered array the VID keeps them and survives with a
-- margin of one, and if the array empties instead, that is a genuine case-4 representation
-- finding rather than a rounding artifact. It is a SNP, so left-alignment cannot be the
-- cause and the finding would be real. Watch this VID specifically in STEP 6.5.
--
-- WHAT THIS DOES NOT SETTLE. Whether the `AN/2` participants are homozygous reference or are
-- carriers of a DIFFERENT alt allele at the site is not distinguishable from these three
-- columns, and the two readings differ: 16-74491696-C-G has 2,208 passing participants but
-- `AC` = 2, so either ~2,206 people are genuinely reference there (0.4% of the cohort) or
-- they carry something else. QUERY D1's `partners` column answers it -- do not guess from
-- `AN` alone.

-- QUERY D1 -- RUN THIS ONE FIRST. One row per VID, verdict column included.
--
-- QUERY D2 below returns one row per ALLELE at the site, which is its payload and cannot be
-- collapsed without losing the thing it is for. But the PREDICTION is a per-VID yes/no, so
-- evaluating it per VID is strictly easier to read across 20 of them, and the detail is then
-- only needed for whichever VIDs come back `NO`. Same scan, so running D1 first costs
-- nothing.
--
-- `own_allele_rows` SHOULD BE 1. A 0 means `filter_set_info` has no row at the VID's own ref
-- and alt -- for the two indels read the left-alignment caveat below before concluding
-- anything, and for a SNP it is a Step 5 defect. Above 1 would mean a duplicate key in
-- `filter_set_info`, which `build_allele_status.sql:90` already warns fans out joins.
WITH v AS (
  SELECT vid FROM UNNEST([
    "7-152451810-T-G",
    "10-70272137-T-G",
    "4-40960447-G-T",
    "16-9293170-T-TCAGGTTTTTCA",
    "2-61866200-C-T",
    "14-89409352-C-G",
    "3-65490858-G-A",
    "12-4195211-C-G",
    "3-94641999-A-G",
    "12-28707196-C-A",
    "10-38959806-A-AAAAG",
    "13-37650636-T-C",
    "16-29292938-T-C",
    "2-55184106-A-T",
    "19-6994372-A-T",
    "7-65019556-G-A",
    "1-27036512-G-A",
    "16-74491696-C-G",
    "2-208094726-G-C",
    "13-80792215-G-T"
  ]) AS vid
),
loc AS (
  SELECT
    vid,
    (CASE SPLIT(vid, "-")[OFFSET(0)]
       WHEN "X" THEN 23
       WHEN "Y" THEN 24
       ELSE CAST(SPLIT(vid, "-")[OFFSET(0)] AS INT64)
     END) * 1000000000000 + CAST(SPLIT(vid, "-")[OFFSET(1)] AS INT64) AS location,
    SPLIT(vid, "-")[OFFSET(2)] AS own_ref,
    SPLIT(vid, "-")[OFFSET(3)] AS own_alt
  FROM v
),
joined AS (
  SELECT
    l.vid,
    (f.ref = l.own_ref AND f.alt = l.own_alt)                                   AS is_own,
    f.yng_status,
    f.calibration_sensitivity,
    f.calibration_sensitivity > IF(LENGTH(f.alt) = LENGTH(f.ref), 0.997, 0.990) AS over_threshold
  FROM loc AS l
  JOIN `foxtrot.filter_set_info` AS f USING (location)
  WHERE f.filter_set_name = "foxtrot_v4_2025_07_29"
)
SELECT
  vid,
  COUNT(*)                                                        AS alleles_at_site,
  COUNTIF(is_own)                                                 AS own_allele_rows,
  ANY_VALUE(IF(is_own, yng_status, NULL))                         AS own_yng,
  ANY_VALUE(IF(is_own, calibration_sensitivity, NULL))            AS own_cal_sens,
  LOGICAL_OR(is_own AND over_threshold)                           AS own_over_threshold,
  COUNTIF(NOT is_own)                                             AS partners,
  COUNTIF(NOT is_own AND yng_status = "Y")                        AS partners_yes,
  IF(COUNTIF(is_own) = 1
       AND LOGICAL_OR(is_own AND over_threshold)
       AND LOGICAL_OR(is_own AND yng_status = "G")
       AND COUNTIF(NOT is_own AND yng_status = "Y") >= 1,
     "predicted", "NO -- see D2")                                 AS verdict
FROM joined
GROUP BY vid
ORDER BY verdict, vid;

-- -------------------------------------------------------------------------------------
-- RESULT, 2026-09-05 -- 19 of 20 `predicted`. The mechanism is confirmed.
-- -------------------------------------------------------------------------------------
-- THREE THINGS HOLD IN ALL TWENTY ROWS, before the verdict column is even read:
--
--   `own_allele_rows` = 1     Every VID found exactly one `filter_set_info` row at its own
--                             ref and alt. No missing rows, no duplicate keys -- and BOTH
--                             INDELS FOUND THEIRS, so the left-alignment caveat below did
--                             not bite and the 200bp widening is not needed.
--   `own_yng` = "G"           Not one of these near-fixed alleles is catalogued as a
--                             positive training site. 20 for 20.
--   `alleles_at_site` >= 2    Every one is multi-allelic; the range is 2 to 39, median 6.
--
-- AND `partners_yes` >= 1 IN 19 OF 20. So at nineteen of the twenty most extreme sites in
-- the callset, the site holds an UNCATALOGUED near-fixed allele scored `G` alongside a
-- CATALOGUED one scored `Y`. That is the matched pair the corrected mechanism predicted,
-- confirmed nineteen times independently. Your reading of the truth resources is the whole
-- explanation: the maintainers have catalogued something at these positions, just not the
-- allele that nearly everyone carries.
--
-- THE ONE EXCEPTION, AND IT IS NOT A DEFECT: `10-38959806-A-AAAAG`. Its own allele scores
-- 0.9897 against the 0.990 indel threshold -- UNDER by three ten-thousandths, so it passes
-- -- and it has no `Y` partner. The corrected mechanism cannot be what excludes 534,745
-- people there. The ORIGINAL mechanism can: with the VID's own allele ok, a `0/1` carrier
-- passes on `all_ok`, but a `1/2` carrier whose partner is over threshold fails, and both
-- partners here are `G`. Consistent with QUERY C, where this is the only VID whose `AC` (835)
-- EXCEEDS its `AN/2` (829) -- at least six homozygotes, and essentially every passing
-- participant carries the allele, exactly as an allele that passes on its own merits should
-- look. STEP 5 IS NOT WRONG HERE; the site simply fails by the other route.
--
-- WHICH MEANS THE HARDY-WEINBERG ARGUMENT UNDER QUERY A WAS SCOPED TOO BROADLY. It rules out
-- the partner mechanism at LOW-ARITY sites only. At 7-152451810-T-G, `alleles_at_site` = 2,
-- so heterozygosity caps at 50% and the argument is airtight. It says nothing at a 39-allele
-- site, where compound heterozygosity can be the norm. The nineteen do not depend on it --
-- they are confirmed directly by `own_over_threshold` and `partners_yes` -- but do not carry
-- the HW claim forward as though it applied to every extreme VID. It does not.
--
-- ============ EVERY ONE OF THESE FAILS BY A HAIRSBREADTH, WHICH IS WORTH ============
-- ============ RAISING WITH THE TEAM SEPARATELY FROM THIS CORRECTION      ============
--
-- The margins over threshold, ascending: +0.0001 four times (12-4195211-C-G, 13-80792215-G-T,
-- 14-89409352-C-G, 3-94641999-A-G), +0.0002, +0.0003 twice, +0.0004, +0.0006 twice, +0.0009,
-- +0.0010, +0.0011, +0.0013 three times, +0.0016, +0.0022, and +0.0023 at the very top
-- (10-70272137-T-G at 0.9993). THE LARGEST MARGIN IN THE SET IS 0.0023. Four VIDs fail by one
-- ten-thousandth. A threshold moved from 0.997 to 0.998 would flip most of the top 20 from
-- excluded to included, changing `gvs_all_ac` at each from single digits to ~535,000.
--
-- THIS DOES NOT CHANGE THE CORRECTION, AND MUST NOT BE ALLOWED TO. The mapping table's job is
-- to agree with the VAT that ships beside it, not to be independently right about biology.
-- The VAT says `gvs_all_ac` = 2 at 7-152451810-T-G; a mapping table naming 535,659 carriers
-- contradicts the artifact it indexes, and that is true whether or not the filter model made
-- the right call. Proceed with STEP 6.3 unchanged.
--
-- What it IS is a filter-model finding to hand to the team on its own: VETS is marginally
-- excluding near-fixed alleles at apparent hg38 reference-error sites, on knife-edge
-- calibration sensitivities, and because it is marginal the effect is unstable across model
-- rebuilds. Raise it separately; do not fold it into the mapping-table work.


-- QUERY D2 -- the per-allele detail. Run for whichever VIDs D1 returns `NO` for; narrowing
-- the list below to those is the intended use, the full 20 being left in only so it is
-- runnable as-is. Its trigger condition was met absolutely: QUERY A returned
-- `ft_fail = excluded` in all 20 rows. The site's whole allele structure.
--
-- THE PREDICTION IS THE INVERSE OF WHAT THIS COMMENT USED TO SAY -- see the corrected
-- mechanism under QUERY A's RESULT block. Each of these locations should hold AT LEAST TWO
-- alleles:
--
--   the VID's own allele    scored `G`, and OVER its threshold. This is the one failing
--                           535,659 people, so anything else here is a Step 5 defect.
--   a partner allele        scored `Y`. This is the rescuer; without it `gvs_all_ac` is 0
--                           and the VID could not be in the VAT at all.
--
-- Four ways it can come back negative, each meaning something different:
--
--   VID's allele UNDER threshold      -> it should be passing everyone. Flatly contradicts
--                                        QUERY A. Step 5 defect; go to the coalesce defaults
--                                        at `import_gvs.py:346-383`.
--   no row for the VID's allele       -> missing `filter_set_info` coalesces to PASSING in
--                                        `import_gvs.py`, so it should be passing everyone
--                                        too. Same conclusion, same place to look. For the
--                                        two INDEL VIDs see the left-alignment caveat below
--                                        before concluding anything.
--   only one allele, over threshold   -> then NOBODY can pass, `gvs_all_ac` should be 0, and
--                                        the VID should not exist in the VAT. Either QUERY C
--                                        disagrees -- in which case Step 5 is wrong -- or
--                                        this VID's array EMPTIES under the patch and it is
--                                        a real case-4 finding of the kind STEP 6.5 hunts.
--   a partner exists but is `G`       -> it cannot rescue anyone. Same dead end as above:
--                                        no survivors, so reconcile against QUERY C.
--
-- LEFT-ALIGNMENT CAVEAT, and it bites exactly two of the 20. The join is on `location` only,
-- decoded from the VID string -- but a VID is left-aligned while `filter_set_info` is keyed
-- on the non-left-aligned `(location, ref, alt)` from `alt_allele`. For a SNP the two
-- coincide, since a SNP is already minimal and left-aligned, so 18 of the top 20 join
-- exactly. `16-9293170-T-TCAGGTTTTTCA` and `10-38959806-A-AAAAG` are INDELS and their
-- `filter_set_info` rows may sit at a different position entirely. If either returns no row
-- for its own allele, widen to `BETWEEN location AND location + 200` before calling it a
-- defect -- 200bp rightward is the window the production synonym search uses.
--
-- Thresholds are 0.997 SNP / 0.990 indel; class is `LENGTH(alt) = LENGTH(ref)`. Note this
-- is the per-ALLELE rule, which is the one that governs the VDS. Do not apply the
-- min-over-class SITE rule from `ExtractCohortEngine` -- see the note under P2's RESULT
-- block; it does not reach the VDS or the VAT.
--
-- The filter set name is `foxtrot_v4_2025_07_29`, the same one Steps 3, 4 and 5 all used.
-- It is not optional: `filter_set_info` holds every model ever built for this dataset and
-- omitting the predicate silently mixes them. If the bytes estimate is large, the join is
-- not pruning partitions; compute the locations from QUERY A by hand and paste them as an
-- integer literal list instead.
WITH v AS (
  SELECT vid FROM UNNEST([
    "7-152451810-T-G",
    "10-70272137-T-G",
    "4-40960447-G-T",
    "16-9293170-T-TCAGGTTTTTCA",
    "2-61866200-C-T",
    "14-89409352-C-G",
    "3-65490858-G-A",
    "12-4195211-C-G",
    "3-94641999-A-G",
    "12-28707196-C-A",
    "10-38959806-A-AAAAG",
    "13-37650636-T-C",
    "16-29292938-T-C",
    "2-55184106-A-T",
    "19-6994372-A-T",
    "7-65019556-G-A",
    "1-27036512-G-A",
    "16-74491696-C-G",
    "2-208094726-G-C",
    "13-80792215-G-T"
  ]) AS vid
),
loc AS (
  SELECT
    vid,
    (CASE SPLIT(vid, "-")[OFFSET(0)]
       WHEN "X" THEN 23
       WHEN "Y" THEN 24
       ELSE CAST(SPLIT(vid, "-")[OFFSET(0)] AS INT64)
     END) * 1000000000000 + CAST(SPLIT(vid, "-")[OFFSET(1)] AS INT64) AS location
  FROM v
)
SELECT
  l.vid                                        AS looking_at,
  f.location,
  f.ref,
  f.alt,
  IF(LENGTH(f.alt) = LENGTH(f.ref), "SNP", "INDEL")                       AS class,
  f.calibration_sensitivity,
  f.yng_status,
  f.calibration_sensitivity > IF(LENGTH(f.alt) = LENGTH(f.ref), 0.997, 0.990) AS over_threshold
FROM loc AS l
JOIN `foxtrot.filter_set_info` AS f USING (location)
WHERE f.filter_set_name = "foxtrot_v4_2025_07_29"
ORDER BY f.location, over_threshold DESC, f.calibration_sensitivity;


-- =====================================================================================
-- STEP 6.1 -- build a 1% sample of the mapping table. 0.01 x T, ~$1.20.
-- =====================================================================================
--
-- TABLESAMPLE really does read fewer blocks, unlike LIMIT. At 1% this is ~54,000 affected
-- VIDs, which is far more than enough to establish the execution shape and to run both
-- validations at full fidelity for the price of a coffee.
--
-- IT EXERCISES ALL THREE OF STEP 6.3's BRANCHES, which was not true of the two-branch
-- design this replaced. 1% of the 152,924 heavy VIDs is ~1,529, enough that the explode ->
-- anti-join -> reaggregate path shows up in the plan with real volume behind it and its
-- exploded row count extrapolates x100. That number is the one thing about STEP 6.3 still
-- unmeasured: `sum(delivered)` over the heavy VIDs is not knowable from `exclusions`, only
-- bounded, so read it off the rehearsal's plan before committing to the full build.
--
-- Two TABLESAMPLE restrictions, both of which reject the query outright rather than
-- degrading: the sampled table may not be referenced elsewhere in the same query, and
-- sampling may not appear inside an IN/EXISTS subquery. The alias must come BEFORE
-- TABLESAMPLE; the other order is a syntax error.

CREATE OR REPLACE TABLE `foxtrot.vid_to_participant_mapping_2026_08_28_sample` AS
SELECT * FROM `foxtrot.vid_to_participant_mapping_2026_08_28` AS m TABLESAMPLE SYSTEM (1 PERCENT);


-- =====================================================================================
-- STEP 6.2 -- rehearse on the sample. ~$3.60 all in.
-- =====================================================================================
--
-- Run STEP 6.3 verbatim against `vid_to_participant_mapping_2026_08_28_sample`, then the two validations below
-- against its output. This is the step that decides whether the full-scale build is safe,
-- and the reason it exists is the scrub's history: its first attempt at the analogous
-- statement held all 2,000 on-demand slots for 82 minutes and did not finish. Correctness
-- was never the problem there; execution shape was.
--
-- Extract STEP 6.3's statement -- the DDL and the query, none of the comment block above
-- it -- into `step6_section_63.sql`, then:
--
--   sed -e 's/vid_to_participant_mapping_2026_08_28/vid_to_participant_mapping_2026_08_28_sample/g' \
--       -e 's/vid_to_participant_mapping_2026_09_05/vid_to_participant_mapping_2026_09_05_sample/g' \
--       step6_section_63.sql > /tmp/step6_sample.sql
--   bq query --nouse_legacy_sql --project_id="${PROJECT}" --label=step6:rehearsal < /tmp/step6_sample.sql
--
-- BOTH EXPRESSIONS ARE REQUIRED. With only the first, the rehearsal reads the 1% sample and
-- then CREATE OR REPLACEs the REAL `2026_09_05` output table with 1% of the data -- silently,
-- and it looks like it worked.
--
-- RUN THE SED EXACTLY ONCE, ON PRISTINE TEXT. The rehearsal names are now SUFFIXES of the
-- production names, which is what makes them readable but also makes this substitution
-- non-idempotent: a second pass yields `..._sample_sample` and the query fails on a missing
-- table. Re-render from `step6_section_63.sql` rather than re-running sed on its own output.
--
-- Two checks on the rendered file, both of which must come back EMPTY:
--
--   grep -nE '`foxtrot\.vid_to_participant_mapping_2026_0[0-9]_[0-9]{2}`' /tmp/step6_sample.sql
--   grep -n  '_sample_sample' /tmp/step6_sample.sql
--
-- The first catches a bare name that survived the rewrite; the second catches sed applied
-- twice. Do not use the old `grep 2026_0` check -- the rehearsal names now contain the date
-- too, so it matches everything and tells you nothing.
--
-- THE FIRST CHECK MUST BE ANCHORED ON BACKTICKS, and the looser form that was here before --
-- `mapping_2026_0[0-9]_[0-9]{2}[^_]` -- is a FALSE-POSITIVE GENERATOR. It fires on the
-- `OPTIONS(description = ...)` string, which says "Derived from
-- vid_to_participant_mapping_2026_08_28 ..." in prose. That is documentation, not a table
-- reference, and rewriting it would be wrong. Only a backtick-quoted name is a reference, so
-- match the backticks.
--
-- WHAT TO READ OFF THE JOB, not just the result:
--
--   * Shuffle bytes. The design intent is that 99.66% of rows pass through the UNION ALL's
--     first branch with a ~110 MB broadcast anti-join and no repartitioning, and only the
--     0.34% affected rows are shuffled by `vid` to meet their exclusion arrays. If the
--     execution graph shows the mapping table's `person_ids` being repartitioned wholesale,
--     the full-scale run will shuffle ~21 TB -- a third of what killed the scrub's STEP 3,
--     and not worth finding out about at $119. The fallback is under STEP 6.3.
--   * Slot-milliseconds, scaled by 100. Compare against the scrub's successful 18m33s run,
--     which did strictly more per row than this does.
--   * `n_removed` distribution, which previews STEP 6.4.
--
-- (a) HOW TO ACTUALLY READ THEM. The console's execution graph shows this, but it is easier
-- more quotable to pull it from INFORMATION_SCHEMA. Run it straight after the rehearsal
-- build. Not literally free -- the JOBS views bill a 10 MB minimum and this dry-runs at a
-- 666 MB upper bound -- but that is fractions of a cent. `region-us`; adjust if the foxtrot
-- dataset lives elsewhere.

SELECT
  job_id,
  TIMESTAMP_DIFF(end_time, start_time, SECOND)                             AS wall_secs,
  ROUND(total_bytes_billed / POW(2, 40), 3)                                AS tib_billed,
  ROUND(total_slot_ms / 1000 / 60, 1)                                      AS slot_minutes,
  ROUND(total_slot_ms / 1000 / 60 * 100, 0)                                AS projected_full_slot_min,
  ARRAY_LENGTH(job_stages)                                                 AS stages,
  ROUND((SELECT SUM(s.shuffle_output_bytes)         FROM UNNEST(job_stages) AS s) / POW(2, 30), 1) AS shuffle_gib,
  ROUND((SELECT SUM(s.shuffle_output_bytes_spilled) FROM UNNEST(job_stages) AS s) / POW(2, 30), 1) AS spilled_gib
FROM `region-us`.INFORMATION_SCHEMA.JOBS_BY_PROJECT
WHERE creation_time > TIMESTAMP_SUB(CURRENT_TIMESTAMP(), INTERVAL 12 HOUR)
  AND statement_type = "CREATE_TABLE_AS_SELECT"
  AND query LIKE "%vid_to_participant_mapping_2026_09_05_sample%"
ORDER BY creation_time DESC
LIMIT 5;

-- HOW TO JUDGE WHAT COMES BACK. At 1% the sample holds ~26 billion entries and ~210 GB of
-- `person_ids`, so:
--
--   `shuffle_gib` in the low tens        -> AS DESIGNED. Only the 0.34% affected rows are
--                                           being repartitioned by vid; the unaffected 99.66%
--                                           are taking the broadcast anti-join. PROCEED.
--   `shuffle_gib` approaching ~210       -> THE WHOLE TABLE IS BEING SHUFFLED. The branch
--                                           split is not doing its job and the full-scale run
--                                           will move ~21 TB. Do NOT spend the $119; go to
--                                           STEP 6.3's fallback note.
--   `spilled_gib` > 0                    -> shuffle spilled to disk, which at 1% means the
--                                           full run certainly will. Same conclusion.
--   `projected_full_slot_min`            -> slot-minutes x 100. Sanity-check against the
--                                           scrub's successful 18m33s run, which did strictly
--                                           more work per row than this does. Same order of
--                                           magnitude is a pass; 10x is a warning.
--
-- Slot time is the WEAKER signal of the two -- per the Step 5 audit, on-demand throughput
-- varied 4x within a single overnight run on identical SQL, so do not treat a slow rehearsal
-- as diagnostic on its own. Shuffle bytes are a property of the plan and do not drift.
--
-- -------------------------------------------------------------------------------------
-- THE TOTAL IS NOT SUFFICIENT TO JUDGE BY. Read the PER-STAGE breakdown.
-- -------------------------------------------------------------------------------------
-- `shuffle_output_bytes` summed over stages counts a row once per stage that moves it, so a
-- 69-stage plan legitimately re-shuffles and the total can exceed the table's size several
-- times over without anything being wrong. The threshold stated above -- "approaching 210 GB
-- means the whole table is moving" -- assumed a SINGLE shuffle and is not usable as written.
--
-- WHAT ACTUALLY MATTERS IS *WHICH* DATA MOVES, and `bytes_per_out_record` separates the two
-- cases cleanly because they differ by ~450x:
--
--   ~10,000-13,000 B/record  -> the WIDE payload. Whole `person_ids` arrays (mean 1,600
--                               entries x 8 B) are being repartitioned. THIS is what scales
--                               100x to ~20 TB at full scale. Minimize it.
--   ~20-30 B/record          -> NARROW (vid, person_id) pairs from `exclusions`. Harmless,
--                               and critically it DOES NOT SCALE: `exclusions` is read at
--                               full size in the rehearsal already, because only the mapping
--                               table is sampled. Shuffling those pairs four times costs the
--                               same at 100% as it does at 1%.
--
-- WHICH ALSO MEANS `projected_full_slot_min` (slot-minutes x 100) IS PESSIMISTIC BY
-- CONSTRUCTION. Roughly a quarter of the rehearsal's bytes are `exclusions` and that quarter
-- is already full size, so it does not multiply. Treat the projection as an upper bound.

SELECT
  s.id                                                        AS stage_id,
  s.name,
  s.records_read,
  s.records_written,
  ROUND(s.shuffle_output_bytes / POW(2, 30), 1)               AS shuffle_gib,
  SAFE_DIVIDE(s.shuffle_output_bytes, s.records_written)      AS bytes_per_out_record,
  ROUND(s.slot_ms / 1000 / 60, 1)                             AS slot_minutes
FROM `region-us`.INFORMATION_SCHEMA.JOBS_BY_PROJECT AS j, UNNEST(j.job_stages) AS s
WHERE j.job_id = "<JOB_ID FROM THE QUERY ABOVE>"
  AND j.creation_time > TIMESTAMP_SUB(CURRENT_TIMESTAMP(), INTERVAL 12 HOUR)
  AND s.shuffle_output_bytes > 0
ORDER BY s.shuffle_output_bytes DESC
LIMIT 15;

-- -------------------------------------------------------------------------------------
-- RESULT, RUN 1, 2026-09-05. job_aMsSvfEPmh7ZUUTtk3Jlst8ZStMe. VERDICT: DO NOT RUN 6.3
-- AS WRITTEN. Rebuild the exclusion side as real tables (STEP 6.2b) and rehearse again.
-- -------------------------------------------------------------------------------------
--   wall_secs 44 | tib_billed 0.251 | slot_minutes 1049.9 | stages 69
--   shuffle_gib 1319.3 | spilled_gib 0.0
--
-- The three clean signals all read WELL. Billed bytes landed exactly where predicted
-- (~210 GB for the 1% sample + the ~66 GB `exclusions` scan). Zero spill. And it FINISHED,
-- in 44 seconds -- the scrub's failed attempt held 2,000 slots for 82 minutes and did not.
-- Whatever is wrong here is not the scrub's failure mode recurring.
--
-- THE PER-STAGE BREAKDOWN SPLITS EXACTLY ON THE 450x BOUNDARY PREDICTED ABOVE, with nothing
-- ambiguous in between -- every stage is either ~14,000-19,000 B/record or ~34 B/record:
--
--   stage             records_written  shuffle_gib  B/record   what it is
--   S44 Compute            16,015,400        215.7    14,462   WIDE, whole sample
--   S40 Join+              16,015,400        215.7    14,462   WIDE, whole sample
--   S37 Input              16,015,401        215.7    14,459   WIDE, whole sample
--   S2D Input              16,015,401        215.7    14,459   WIDE, whole sample
--   S06 Join+           2,344,719,084         74.2        34   narrow, `ex_pairs`
--   S18 Aggregate       2,344,718,153         74.2        34   narrow, `ex_pairs` DISTINCT
--   S41 Repartition         2,590,115         46.0    19,087   WIDE, unexplained -- see below
--   S0E Aggregate         518,209,967         41.1        85   narrow
--   S0D Aggregate         518,209,967         17.1        36   narrow
--   ...nine more Repartition stages, all ~34 B/record, 11-17 GiB each
--
-- SO 862.8 GiB OF THE 1,319.3 IS THE ENTIRE SAMPLED TABLE'S `person_ids`, MOVED FOUR TIMES.
-- `records_written` on those four stages is 16,015,401 -- the whole 1% sample, not the
-- ~54,000 affected rows the three-branch split was built to isolate. Only that 862.8 scales:
-- x100 is ~86 TB at full size, against ~21 TB for the scrub's one known-good full-table
-- shuffle. The narrow stages are `exclusions` at full size already and do not multiply.
--
-- ONE OF THE FOUR IS INHERENT AND MUST BE BUDGETED FOR, WHICH CORRECTS SOMETHING SAID
-- EARLIER IN THIS FILE. `CLUSTER BY vid` cannot be satisfied without repartitioning the
-- output by `vid`, so it costs one full-table WIDE shuffle -- 215.7 GiB here, ~21.6 TB at
-- full scale. It is free in DOLLARS, which is what P1 was talking about, but it is not free
-- in execution shape and should never have been described as simply "free". It buys the
-- delivered table a physical property it lost, at the cost of exactly the shuffle the scrub
-- already demonstrated is survivable. Keep it; just do not count it as zero.
--
-- THE OTHER THREE ARE THE FAILURE, AND THIS FILE ALREADY PREDICTED IT. See the paragraph in
-- STEP 6.3 beginning "IF THE SAMPLE REHEARSAL SHOWS THE PROBE SIDE BEING REPARTITIONED
-- ANYWAY". Two bare `Input` stages reading 16,015,401 rows and shuffling all 215.7 GiB with
-- no reduction (records_read == records_written) are branches 1 and 2 repartitioning the
-- probe side by `vid` instead of streaming past a broadcast build side. `ex_vids` really is
-- only ~110 MB, but it is a CTE reached through `SELECT DISTINCT` over 2.34 BILLION rows,
-- and BigQuery will not stake a broadcast on a cardinality it can only guess at. It is not
-- an optimizer defect; it is missing statistics, and the remedy is to supply them.
--
-- `ex_pairs` IS ALSO BEING RE-EVALUATED, exactly as the note in STEP 6.3 feared. S0D and
-- S0E are two SEPARATE aggregate stages each reading all 2,344,719,084 rows, on top of S06
-- and S18. Cheap in dollars, but it is real repeated work and materializing removes it too.
--
-- FREE RESULT WORTH BANKING: `records_written` ON S18 IS 2,344,718,153 EXACTLY -- the low
-- end of STEP 6.4's reconciliation window, confirmed directly off the query plan rather
-- than inferred. So the carve-out leaves 2,344,719,084 pairs, of which DISTINCT keeps
-- 2,344,718,153, meaning exactly 931 duplicate pairs collapse. VS-2010 counted 932
-- duplicated pairs genome-wide, so precisely one of them sits on a fixup VID and was carved
-- out before the DISTINCT saw it. Both numbers tie out, from opposite directions.
--
-- S41 IS NOT YET EXPLAINED AND SHOULD BE RE-READ AFTER THE FIX. 2,590,115 records at 19,087
-- B/record is WIDE and skews wider than average, which fits affected VIDs (they are affected
-- precisely because many people carry them). But 2.59M is 16% of the sample, where the
-- affected fraction is 0.34%, so the row count does not fit any branch. It is 3.5% of the
-- total and not worth chasing ahead of the four 215.7 GiB stages; if it survives STEP 6.2b
-- unchanged, chase it then.
--
-- WHAT "GOOD" LOOKS LIKE ON THE RE-RUN. One wide stage at ~215.7 GiB (the clustered write),
-- the affected branches down at single-GiB, narrow stages roughly unchanged: total in the
-- 650-700 GiB region rather than 1,319.3. Full-scale wide shuffle ~21.6 TB, matching the
-- scrub's known-good shape instead of quadrupling it.
--
-- -------------------------------------------------------------------------------------
-- RESULT, RUN 2, 2026-09-05. bqjob_r55e684abd0f0541c_000001a072a51f37_1, after STEP 6.2b.
-- VERDICT: PROCEED TO STEP 6.3. Not because the target above was met -- it was not -- but
-- because the target was measuring the wrong thing. Read the reasoning before accepting
-- the verdict: the raw numbers look like a failure and are not one.
-- -------------------------------------------------------------------------------------
--   wall_secs 34 | tib_billed 0.253 | slot_minutes 651.8 | stages 40
--   shuffle_gib 1122.2 | spilled_gib 0.0        (run 1: 44 / 0.251 / 1049.9 / 69 / 1319.3)
--
-- WHAT 6.2b BOUGHT, EXACTLY: the entire 197.1 GiB drop is NARROW traffic, and it is the
-- repeated `ex_pairs` aggregation. Run 1 evaluated that CTE across four stages (S06, S18,
-- S0D, S0E -- 206.6 GiB between them); run 2 reads the materialized table once, in S00, for
-- 74.2 GiB. 29 stages and 38% of slot time went with it. That part worked as intended.
--
-- WHAT IT DID NOT BUY: THE WIDE TRAFFIC IS UNCHANGED. Four stages still move the whole
-- sampled table at ~14,460 B/record -- S15 and S19 (bare `Input`, 16,015,401 in and out),
-- S23 (`Join+`) and S27 (`Compute`, the clustered write). So `ex_vids_t` at 137.0 MiB and
-- `light_vids_t` at 133 MiB are BOTH above whatever BigQuery's broadcast threshold actually
-- is -- commonly cited as ~100 MB, and this is the first hard evidence in this file of where
-- it sits. Materializing gave the planner exact statistics and it still declined. The
-- statistics were necessary and not sufficient.
--
-- THE SEMI-JOIN REDUCER DID WORK, and S1A proves it: 52,898 records written at 283,835
-- B/record. 52,898 is the affected VID count at 1% (5,437,495 x 0.01 = 54,375), and 283,835
-- B/record is what an affected VID's array weighs -- these hold the biggest arrays in the
-- table, which is what made them affected in the first place. The affected branches are
-- correctly reduced. What is not reduced is the full-table pass each branch makes BEFORE it.
--
-- ===== WHY THIS IS A GO ANYWAY, AND WHY THE EARLIER FRAMING WAS WRONG =====
--
-- Run 1's analysis called 4 x 215.7 GiB "86 TB at full scale, four times the scrub's
-- known-good 21 TB" and treated that as the risk. IT IS NOT THE RISK. Shuffle exhaustion is
-- a PER-STAGE property: a stage spills or fails when ITS OWN footprint outgrows what the
-- slots holding it can carry. Four sequential stages of 215.7 GiB do not co-exist -- each
-- drains before the next fills. The number to compare against the scrub is therefore the
-- PEAK SINGLE STAGE:
--
--   peak single stage, run 2        215.7 GiB at 1%        ->  21.1 TiB at full scale
--   scrub's successful 18m33s run   one full-table shuffle ->  ~21 TiB
--
-- Identical. The extra passes cost TIME, not headroom, and `spilled_gib` = 0.0 confirms the
-- headroom is there. What quadrupling buys is a longer job, and that is quantifiable:
--
--   naive projection (scale everything x100)   65,180 slot-min  = 33 min at 2,000 slots
--   holding the narrow 10.4% of slot time fixed 34,203 slot-min = 17 min at 2,000 slots
--
-- The truth is between the two, because `exclusions` and the `_t` tables are already full
-- size in the rehearsal and do not scale. Against the scrub's successful 18m33s and its
-- failed 82-minute non-finish, 17-33 minutes sits comfortably inside known-good territory.
--
-- SO THE GATE PASSES ON ITS OWN TERMS. It was never "is this plan optimal", it was "will
-- this fail the way the scrub failed". A peak stage equal to a known-good run, zero spill,
-- and a wall-clock projection bracketing that run's actual time is a no.
--
-- IF SOMEONE WANTS TO OPTIMISE FURTHER ANYWAY, the lever is getting the two key sets under
-- the broadcast threshold, and `FARM_FINGERPRINT(vid)` does it: 5,437,495 INT64s is 43.5 MB
-- against 137.0 MiB. It is NOT free -- a 64-bit collision between an affected and an
-- unaffected VID would drop the row from every branch and lose it silently. Expected
-- collisions are 1.6e9 x 5.44e6 / 2^64 ~ 5e-4, so it is very unlikely, and more usefully it
-- is CHEAPLY CHECKABLE: a query joining hash to hash and re-testing the string reads the
-- `vid` column alone, ~26 bytes a row against `person_ids`' ~12.8 KB, so it costs cents
-- rather than 1 x T. Verify zero collisions first, then hash. Roughly $4 all in, to maybe
-- halve the wall time of a job that already fits.
--
-- (b) V1 -- no excluded pair survives, and no VID was emptied without being counted.
-- Full fidelity on the sample, unlike the full-scale version which is optional.

-- V1a JOINS `ex_pairs_t`, NOT `exclusions`, AND THE DIFFERENCE IS THE WHOLE CHECK.
-- `ex_pairs_t` is by construction the exact set STEP 6.3 is contracted to remove: distinct
-- exclusion pairs with the 5,776 fixup VIDs carved out. `exclusions` is a strict superset of
-- it, so testing against `exclusions` asks 6.3 to have done something it deliberately does
-- not do, and a fixup VID landing in the sample then reports FAIL against a correct build.
-- That is not hypothetical -- it is what happened on 2026-09-05, failing V1a and V2b
-- together off a single cause. Do not "tighten" this back to `exclusions`.
SELECT "V1a: no excluded pair survives" AS check_name,
       IF(NOT EXISTS(
            SELECT 1
            FROM `foxtrot.vid_to_participant_mapping_2026_09_05_sample` AS c
            CROSS JOIN UNNEST(c.person_ids) AS p
            JOIN `foxtrot.ex_pairs_t` AS e
              ON e.vid = c.vid AND e.person_id = p
          ), "pass", "FAIL") AS result
UNION ALL
SELECT "V1b: no empty person_ids array",
       IF(NOT EXISTS(
            SELECT 1 FROM `foxtrot.vid_to_participant_mapping_2026_09_05_sample`
            WHERE ARRAY_LENGTH(person_ids) = 0
          ), "pass", "FAIL")
UNION ALL
SELECT "V1c: n_removed agrees with the arrays",
       IF(NOT EXISTS(
            SELECT 1 FROM `foxtrot.vid_to_participant_mapping_2026_09_05_sample`
            WHERE n_removed != delivered - ARRAY_LENGTH(person_ids)
          ), "pass", "FAIL");

-- (c) V2 -- the output is exactly the delivered table minus its excluded members, and nothing
-- else moved. This is the check that cannot be afforded at full scale (2 x T), which is
-- precisely why it is worth running here.
--
-- Compared as multisets -- both sides sorted -- rather than as sets, so that a build which
-- dropped or added a repeat of a retained person id still fails. `TO_JSON_STRING` over a
-- sorted array is the same device the scrub's STEP 4b used.
--
-- `expected` is built by EXPLODE then ANTI-JOIN then REAGGREGATE, which is a genuinely
-- different formulation from STEP 6.3's per-VID membership test rather than a copy of it --
-- that is most of the value of running this at all. It is also the formulation that cannot
-- be afforded at full scale: at 1% it explodes ~26 billion entries, which is fine, and at
-- 100% it is the 2.58-trillion-row shuffle that killed the scrub's STEP 3.
--
-- Do NOT write it the compact way, as `ARRAY(SELECT p FROM UNNEST(o.person_ids) AS p WHERE
-- NOT EXISTS (SELECT 1 FROM exclusions ...))`. That parses, and then fails at execution
-- against real tables with "Correlated subqueries that reference other tables are not
-- supported unless they can be de-correlated" -- a restriction a CTE-literal fixture does
-- not reproduce, which is exactly how it gets shipped untested.

WITH fixup_vids AS (
  -- MANUAL: two case-4 VIDs were fixed by hand and are NOT in either mapping table. They are
  -- inert here -- see the STEP 6.3 note -- so they are deliberately not listed.
  SELECT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`      GROUP BY vid
  UNION DISTINCT
  SELECT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01` GROUP BY vid
),
ex_pairs_check AS (
  -- THE CARVE-OUT MUST BE HERE. Without it `expected` demands removals STEP 6.3 deliberately
  -- does not perform on the 5,776 fixup VIDs, and V2b reports FAIL against a correct build --
  -- which is exactly what happened on 2026-09-05.
  --
  -- RE-DERIVED FROM `exclusions`, NOT READ FROM `ex_pairs_t`, ON PURPOSE. `ex_pairs_t` holds
  -- precisely this and joining it here would be cheaper and shorter. It would also mean no
  -- check anywhere re-computes 6.2b's contents, so an error in 6.2b would propagate into both
  -- the build and its validation and pass. Independence is the entire value of V2; do not
  -- trade it for the saved aggregation.
  SELECT DISTINCT e.vid, e.person_id
  FROM `foxtrot.exclusions` AS e
  LEFT JOIN fixup_vids AS f ON f.vid = e.vid
  WHERE f.vid IS NULL
),
kept AS (
  SELECT x.vid, x.person_id
  FROM (
    SELECT o.vid, p AS person_id
    FROM `foxtrot.vid_to_participant_mapping_2026_08_28_sample` AS o
    CROSS JOIN UNNEST(o.person_ids) AS p
  ) AS x
  LEFT JOIN ex_pairs_check AS e
    ON e.vid = x.vid AND e.person_id = x.person_id
  WHERE e.vid IS NULL
),
expected AS (
  -- No DISTINCT: multiplicity of a RETAINED person id must be preserved, so that a build
  -- which quietly deduplicated arrays fails V2b rather than passing it.
  SELECT vid, ARRAY_AGG(person_id ORDER BY person_id) AS person_ids
  FROM kept
  GROUP BY vid
)
SELECT "V2a: no vid gained a row" AS check_name,
       IF(NOT EXISTS(
            SELECT 1 FROM `foxtrot.vid_to_participant_mapping_2026_09_05_sample` AS c
            LEFT JOIN `foxtrot.vid_to_participant_mapping_2026_08_28_sample` AS o USING (vid)
            WHERE o.vid IS NULL
          ), "pass", "FAIL") AS result
UNION ALL
SELECT "V2b: every surviving array is the original minus its excluded members",
       IF(NOT EXISTS(
            SELECT 1
            FROM `foxtrot.vid_to_participant_mapping_2026_09_05_sample` AS c
            JOIN expected AS x USING (vid)
            WHERE TO_JSON_STRING(ARRAY(SELECT p FROM UNNEST(c.person_ids) AS p ORDER BY p))
               != TO_JSON_STRING(x.person_ids)
          ), "pass", "FAIL")
UNION ALL
SELECT "V2c: every vid dropped from the output was emptied, not lost",
       IF(NOT EXISTS(
            SELECT 1
            FROM expected AS x
            LEFT JOIN `foxtrot.vid_to_participant_mapping_2026_09_05_sample` AS c USING (vid)
            WHERE c.vid IS NULL AND ARRAY_LENGTH(x.person_ids) > 0
          ), "pass", "FAIL");

-- RESULT, RUN 2026-09-05. V1a and V2b both FAILED; V1b, V1c, V2a and V2c passed. The cause
-- was the missing fixup-VID carve-out in BOTH checks, now fixed above in V1a's join and in
-- `ex_pairs_check`. The two failures were never independent: a surviving excluded pair makes
-- the array differ from `expected`, so one cause fails both. V2a and V2c passing bounded the
-- blast radius immediately -- no VID gained a row, none was lost -- which confined the
-- problem to array CONTENTS before any diagnostic ran.
--
-- V1a FAILING DOES NOT ACCOUNT FOR V2b FAILING, and assuming it does is the trap here. V2b
-- compares multisets, so it also fires if the build dropped a RETAINED person, deduplicated
-- an array, or changed multiplicity -- none of which V1a can see. Every V1a-failing VID is a
-- V2b-failing VID; the reverse does not hold. So the question to answer is not "is V1a
-- benign" but "is the V2b-failing SET exactly the fixup VIDs and nothing else".
--
-- (d) DIAGNOSTIC. Run this when V2b fails, before re-running anything. Same cost as V2
-- itself, because it rebuilds the same `expected`. It reports every mismatching VID with the
-- two facts that classify it:

WITH fixup_vids AS (
  SELECT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`      GROUP BY vid
  UNION DISTINCT
  SELECT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01` GROUP BY vid
),
kept AS (
  SELECT x.vid, x.person_id
  FROM (
    SELECT o.vid, p AS person_id
    FROM `foxtrot.vid_to_participant_mapping_2026_08_28_sample` AS o
    CROSS JOIN UNNEST(o.person_ids) AS p
  ) AS x
  LEFT JOIN (SELECT DISTINCT vid, person_id FROM `foxtrot.exclusions`) AS e
    ON e.vid = x.vid AND e.person_id = x.person_id
  WHERE e.vid IS NULL
),
expected AS (
  SELECT vid, ARRAY_AGG(person_id ORDER BY person_id) AS person_ids
  FROM kept
  GROUP BY vid
),
mismatches AS (
  SELECT c.vid,
         ARRAY_LENGTH(c.person_ids) AS got_len,
         ARRAY_LENGTH(x.person_ids) AS want_len,
         c.delivered,
         c.n_removed
  FROM `foxtrot.vid_to_participant_mapping_2026_09_05_sample` AS c
  JOIN expected AS x USING (vid)
  WHERE TO_JSON_STRING(ARRAY(SELECT p FROM UNNEST(c.person_ids) AS p ORDER BY p))
     != TO_JSON_STRING(x.person_ids)
)
SELECT m.vid,
       f.vid IS NOT NULL       AS is_fixup_vid,
       m.got_len = m.delivered
         AND m.n_removed = 0   AS untouched_by_6_3,
       m.delivered, m.n_removed, m.got_len, m.want_len
FROM mismatches AS m
LEFT JOIN fixup_vids AS f ON f.vid = m.vid
ORDER BY is_fixup_vid, untouched_by_6_3, m.vid
LIMIT 200;

-- NOTE the diagnostic's `kept` deliberately OMITS the carve-out -- it must reproduce the
-- failing comparison in order to explain it. Do not copy the corrected `ex_pairs_check` into
-- it; that would make it report nothing and look like a pass.
--
-- HOW TO READ IT. The ORDER BY puts every bad outcome above every benign one.
--
--   BENIGN, and the only acceptable result: every row has is_fixup_vid = true AND
--     untouched_by_6_3 = true, with got_len > want_len. That is a fixup VID taking branch 1
--     as designed -- n_removed = 0, output array identical to the delivered one, `expected`
--     wanting removals 6.3 was never asked to make.
--   REAL BUG: any row with is_fixup_vid = false. 6.3 mangled a VID it does own.
--   REAL BUG, subtler: is_fixup_vid = true but untouched_by_6_3 = false. A fixup VID went
--     down a branch it should not have, so the carve-out leaked in STEP 6.3 itself.
--   WRONG DIRECTION: got_len < want_len on any row. The build removed someone `expected`
--     kept. This is the serious direction and is never benign.
--   TOO MANY: P3 puts 121 fixup VIDs in `exclusions`, so 1% is ~1.2 rows. A dozen means the
--     carve-out model is wrong, even if every row looks individually benign.
--
-- RESULT, 2026-09-05: ONE ROW, BENIGN, and the rehearsal is therefore green.
--
--   vid                is_fixup  untouched  delivered  n_removed  got_len  want_len
--   22-50028710-A-AC   true      true       60495      0          60495    60480
--
-- Every classifier reads the benign way: the sole mismatch is a fixup VID, it took branch 1
-- untouched, and got_len > want_len is the safe direction. One row against ~1.2 expected.
--
-- THE CORRECTED V1a PASSES BY DEDUCTION, so re-running (b) is optional. A surviving
-- `ex_pairs_t` pair would have to produce a V2b mismatch on a NON-fixup VID -- `ex_pairs_t`
-- has the fixup VIDs carved out -- and (d) found no such row. Re-running (c) is worth the
-- $0.60 only because `ex_pairs_check` is new SQL that has never executed.
--
-- WHAT THIS PROVES AND WHAT IT DOES NOT. It proves STEP 6.3 does exactly what it is
-- contracted to do. It also makes explicit that the contract stops short of a correct table:
-- the 121 fixup VIDs KEEP their exclusions after 6.3. At the 15 pairs seen on
-- 22-50028710-A-AC that is order 1,800 residual over-reported (vid, person) pairs genome-wide
-- -- 7.7e-07 of the correction, but not zero, and concentrated on VIDs that already needed
-- hand-holding. STEP 6a is what closes them, and STEP 6a must therefore run before STEP 8
-- re-exports Parquet. Do not read a green STEP 6.4/6.5/6.6 as "the table is fixed".
--
-- STEP 6.4's BOUNDS ALREADY ACCOUNT FOR THIS, which is a useful independent confirmation:
-- its lower bound 2,344,718,153 is exactly `ex_pairs_t`'s row count -- the carve-out applied
-- -- and its upper bound 2,344,719,085 is that plus VS-2010's 932 duplicate pairs. Neither
-- bound needs adjusting for the fixup VIDs; they were never in the expectation.


-- =====================================================================================
-- STEP 6.2b -- materialize the exclusion side as real clustered tables. ~$1.00 total.
-- REQUIRED. Run 1 of the rehearsal showed the CTE form does not give the optimizer
-- enough to plan a broadcast, and the whole three-branch design depends on one.
-- =====================================================================================
--
-- WHAT THIS BUYS, IN ONE SENTENCE: exact cardinality and byte-size statistics on the build
-- side of every join in STEP 6.3, so BigQuery can see that `ex_vids` is 5.4 million short
-- strings and broadcast it, instead of guessing at the output of a `SELECT DISTINCT` over
-- 2.34 billion rows and defensively repartitioning 21 TB of `person_ids` to be safe.
--
-- IT ALSO REMOVES THE REPEATED AGGREGATION. In run 1, `ex_pairs` was aggregated over all
-- 2,344,719,084 rows in at least four separate stages (S06, S18, S0D, S0E). Written once to
-- a table, it is read once per consumer and never re-derived.
--
-- THESE ALL READ `exclusions` AND NOTHING ELSE, so none of them touches `person_ids` and
-- none costs anything near 1 x T. That is the whole reason this is affordable: the note in
-- STEP 6.3 warning against a pre-materialized "affected subset" table still stands, and for
-- exactly the same reason -- such a table would select `person_ids` and cost a second T.
-- Materializing the EXCLUSION side is the cheap half of that trade, not the expensive half.
--
-- CLUSTER BY vid on all three: it is the join key everywhere downstream, and clustering the
-- probe-side inputs lets BigQuery prune rather than scan when it does choose a hash join.
--
-- RUN ALL FIVE STATEMENTS, IN ORDER. `ex_pairs_t` must exist before the other four.

CREATE OR REPLACE TABLE `foxtrot.ex_pairs_t`
CLUSTER BY vid
OPTIONS (description = "Step 6.2b scratch: distinct (vid, person_id) exclusion pairs with the 5,776 Step 6a fixup VIDs carved out. Materialized so Step 6.3's joins get real statistics. Drop after Step 6.8.")
AS
WITH fixup_vids AS (
  SELECT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`          GROUP BY vid
  UNION DISTINCT
  SELECT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01`     GROUP BY vid
)
SELECT DISTINCT e.vid, e.person_id
FROM `foxtrot.exclusions` AS e
LEFT JOIN fixup_vids AS f ON f.vid = e.vid
WHERE f.vid IS NULL;

-- EXPECT EXACTLY 2,344,718,153 ROWS. That is not an estimate: run 1's S18 stage wrote
-- precisely that many. `SELECT COUNT(*)` on a table is free from metadata, so check it, and
-- if it disagrees do not proceed -- something about the carve-out or `exclusions` changed.

CREATE OR REPLACE TABLE `foxtrot.light_ex_t`
CLUSTER BY vid
OPTIONS (description = "Step 6.2b scratch: exclusion arrays for VIDs with <= 1000 excluded people. Step 6.3 branch 2 build side. Drop after Step 6.8.")
AS
SELECT vid, ARRAY_AGG(person_id) AS ids
FROM `foxtrot.ex_pairs_t`
GROUP BY vid
HAVING COUNT(*) <= 1000;

-- ADDED AFTER MEASURING `light_ex_t` AT 1,313.1 MiB ON 2026-09-05. That is nine times too
-- large to broadcast, so a direct `mapping JOIN light_ex_t` leaves BigQuery no choice but to
-- repartition the PROBE side -- the whole mapping table -- because it cannot know which rows
-- will match until it has joined them. That is one of run 1's four 215.7 GiB stages and
-- materializing the arrays does not fix it; the arrays are the problem.
--
-- SO GIVE BRANCH 2 A VID-ONLY KEY SET TO MEET THE MAPPING TABLE WITH. At ~133 MiB this is
-- the same order as `ex_vids_t` and broadcastable on the same terms, and it reduces the
-- mapping table to the ~54,000 light-affected rows BEFORE the arrays are ever mentioned.
-- The second join, from those 54,000 rows to `light_ex_t`, has a probe side small enough
-- that shuffling both sides costs nothing worth counting.
--
-- The 1,313.1 MiB is not waste, incidentally -- it is 154.6M excluded entries x 8 bytes plus
-- the vid strings, which is exactly what the light branch has to consume. The mistake was
-- making the table carrying it the thing that first touches 21 TB of `person_ids`.
CREATE OR REPLACE TABLE `foxtrot.light_vids_t`
CLUSTER BY vid
OPTIONS (description = "Step 6.2b scratch: keys only of VIDs with <= 1000 excluded people. Step 6.3 branch 2's SEMI-JOIN REDUCER -- exists so the mapping table meets a broadcastable key set before it meets light_ex_t's 1.3 GiB of arrays. Drop after Step 6.8.")
AS
SELECT vid
FROM `foxtrot.ex_pairs_t`
GROUP BY vid
HAVING COUNT(*) <= 1000;

CREATE OR REPLACE TABLE `foxtrot.heavy_ex_t`
CLUSTER BY vid
OPTIONS (description = "Step 6.2b scratch: keys of VIDs with > 1000 excluded people. Step 6.3 branch 3 build side. Drop after Step 6.8.")
AS
SELECT vid
FROM `foxtrot.ex_pairs_t`
GROUP BY vid
HAVING COUNT(*) > 1000;

CREATE OR REPLACE TABLE `foxtrot.ex_vids_t`
CLUSTER BY vid
OPTIONS (description = "Step 6.2b scratch: the ~5.4M distinct affected VIDs. Step 6.3 branch 1 anti-join build side; must be small enough to broadcast. Drop after Step 6.8.")
AS
SELECT vid
FROM `foxtrot.ex_pairs_t`
GROUP BY vid;

-- SANITY, FREE FROM METADATA. `light_ex_t` + `heavy_ex_t` must equal `ex_vids_t`, and P2 put
-- those at roughly 5,275,000 / 152,924 / 5,428,000. A `heavy_ex_t` far off 152,924 means the
-- HEAVY_VID_THRESHOLD reasoning needs revisiting before branch 2's quadratic test runs.
-- `ex_pairs_t` must read EXACTLY 2,344,718,153 -- run 1's plan reported that figure, so this
-- is an equality test, not a range.
--
-- SIZE MATTERS MORE THAN ROW COUNT HERE. Broadcast eligibility is about bytes, and the
-- ~110 MB estimate for `ex_vids_t` is the load-bearing assumption of this whole step. VIDs
-- are variable-length strings and the estimate came from an average; if it comes back at
-- several hundred MB, expect the anti-join to stay a shuffle and see the fallback at the end
-- of STEP 6.3's shape notes.

SELECT
  table_id,
  row_count,
  ROUND(size_bytes / POW(2, 20), 1) AS mib
FROM `foxtrot.__TABLES__`
WHERE table_id IN ("ex_pairs_t", "light_ex_t", "heavy_ex_t", "ex_vids_t", "light_vids_t")
ORDER BY row_count DESC;

-- USE `<dataset>.__TABLES__`, NOT `region-us`.INFORMATION_SCHEMA.TABLE_STORAGE. The
-- region-qualified view is PROJECT-scoped: it lists only datasets in the project the query
-- runs in, and only if the region argument matches the dataset's location. Asking it for
-- `foxtrot` returned zero rows on 2026-09-05 -- no error, just nothing, which is the worst
-- possible failure for a sanity check since an empty result reads like "the tables are
-- missing". `__TABLES__` is dataset-qualified and resolves wherever the dataset actually
-- lives. Note this does NOT invalidate the `region-us`.INFORMATION_SCHEMA.JOBS_BY_PROJECT
-- queries elsewhere in this file: JOBS is scoped by where the JOB ran, which is the current
-- project, and those work. It is the TABLE-scoped views that need the dataset qualifier.


-- =====================================================================================
-- STEP 6.3 -- BUILD the corrected copy alongside the original. 1 x T, ~$119.
-- Does not touch the original.
-- =====================================================================================
--
-- THE SHAPE, and why it is a three-way UNION ALL rather than one LEFT JOIN.
--
-- The natural formulation is `mapping LEFT JOIN per_vid_exclusions USING (vid)` with an
-- `IF(ids IS NULL, person_ids, ARRAY(...))`. It is correct, and it is a trap: the exclusion
-- side carries 2.34 billion INT64s, roughly 18 GB of array payload (P2 measures it), which
-- is too large for BigQuery to broadcast. A shuffle join then repartitions BOTH sides by
-- `vid` -- moving all ~21 TB of `person_ids` across the network to discover that 99.66% of
-- it had no match. That is the same failure mode as the scrub's STEP 3, at a third of the
-- volume.
--
-- Splitting on the join outcome fixes it, because the branches need different things. The
-- unaffected branch needs only to know that a `vid` is absent from a 5.4-million-entry key
-- set -- ~110 MB, comfortably broadcastable, so those rows stream through with no
-- repartitioning at all. The affected branches are 0.34% of VIDs, and P2's numbers put
-- their entire array payload on the order of 200 GB, which shuffles without incident.
--
-- All three branches reference `person_ids`, and bytes are billed per column per statement,
-- so the UNION ALL costs 1 x T and not 3 x T. This is the same accounting the scrub relied
-- on to run its validations as a single statement.
--
-- WHY THE AFFECTED SIDE IS ITSELF SPLIT, at HEAVY_VID_THRESHOLD = 1000 excluded people.
--
-- `p NOT IN UNNEST(x.ids)` is a LINEAR SCAN of `ids` per probe, not a hash lookup, so a
-- single affected row costs `delivered x excluded` element comparisons and that whole cost
-- lands on ONE slot. It cannot be subdivided, so no amount of slot capacity helps.
--
-- P2 measured the tail, and it is not survivable. `max_excluded_per_vid` is 535,659 against
-- an active non-control cohort of 535,662 (`sample_info`, 2026-09-05) -- three short of
-- every participant alive in the callset. The top 20 all sit between 532,570 and 535,659,
-- and 686 VIDs exceed 250,000. The worst single row is therefore ~2.9e11 comparisons; the
-- scrub's one calibration point (a ~450,000-element array, 12+ minutes) puts BigQuery near
-- 1e8 comparisons/sec/slot, which makes that row alone roughly an hour of unsplittable
-- work. Summed, `sum_excluded_sq_e12` = 247.27 is a LOWER bound of 2.47e14 comparisons,
-- since `delivered >= excluded` for every VID.
--
-- So the heavy VIDs go through explode -> hash anti-join -> reaggregate instead, which is
-- linear rather than quadratic and parallelizes across slots. That is the formulation
-- already written and exercised in STEP 6.2's V2 `expected` CTE, so it is known-good SQL.
--
-- 1000 is where P2's `by_threshold` puts the knee. Above it sit 152,924 VIDs -- 2.8% of
-- affected VIDs, 0.0096% of the table -- carrying 93.4% of all excluded entries
-- (2,190,144,230 of 2,344,755,855). Routing only those caps the array test's worst possible
-- row at 535,662 x 1000 ~ 5.4e8, about five seconds, and those rows are spread across 5.28
-- million VIDs averaging 29 excluded each, so they parallelize freely. Moving the threshold
-- lower buys little and grows the explode side; moving it higher reintroduces the quadratic
-- tail the split exists to remove.
--
-- THE HEAVY BRANCH IS NOT ROW-LOCAL, unlike the array test: it groups by `vid`, so two rows
-- sharing a `vid` would be merged into one. P3's `duplicate_vid_rows` must be zero before
-- this runs. It is expected to be -- `vid` is the delivered table's primary key -- but the
-- array-test formulation tolerated a violation silently and this one does not, so P3 is now
-- a hard gate rather than a diagnostic.
--
-- `ex_pairs` IS REFERENCED FOUR TIMES, and BigQuery does not guarantee a CTE is evaluated
-- once. Re-scanning `exclusions` costs ~$0.30 a time, so the dollars do not matter, but the
-- repeated aggregation might. If the plan shows it, materialize `ex_pairs`, `light_ex` and
-- `heavy_ex` as real tables clustered by `vid` first -- all three come from `exclusions`
-- alone and so are cheap -- and reference those instead.
--
-- IF THE SAMPLE REHEARSAL SHOWS THE PROBE SIDE BEING REPARTITIONED ANYWAY, the optimizer
-- has collapsed the branches back into one join. Force the issue by materializing
-- `ex_vids` and `light_ex` as real clustered tables first (cheap, from `exclusions` alone)
-- so the optimizer has exact statistics for both, and re-read the plan. Do not reach for a
-- pre-materialized "affected subset" table: selecting it reads `person_ids` and so costs a
-- second 1 x T, which is the whole thing this restructuring was for.
--
-- THAT IS EXACTLY WHAT RUN 1 SHOWED, on 2026-09-05: four stages at 215.7 GiB and ~14,462
-- bytes per output record, each writing all 16,015,401 sample rows. The CTEs below have
-- therefore been REPLACED BY THE `_t` TABLES BUILT IN STEP 6.2b, which must be run first.
-- The SQL is otherwise unchanged -- same predicate, same branches, same semantics -- so
-- nothing about correctness moved; only the statistics available to the planner did.
--
-- IF THE `_t` TABLES ARE NOT ENOUGH EITHER and branch 1 still repartitions the probe side,
-- the remaining lever is to drop `CLUSTER BY vid` from branch 1's perspective by writing the
-- three branches to three separate tables and using `bq cp --append`, which costs nothing to
-- read and lets the unaffected 99.66% be written without ever entering a join. That is more
-- moving parts and more ways to deliver a partial table, so it is a last resort, not a
-- default -- but it is the shape that cannot be talked out of streaming.
--
-- `CLUSTER BY vid` and no `ORDER BY` -- see P1. No `PARTITION BY`: `vid` is a STRING and
-- there is no integer or date column to range- or time-partition on, which is why the
-- delivered table has never been partitioned either.

CREATE OR REPLACE TABLE `foxtrot.vid_to_participant_mapping_2026_09_05`
CLUSTER BY vid
OPTIONS (description = "Step 6: Foxtrot participant mapping with GQ 0 no-calls and FT-failing genotypes removed. Derived from vid_to_participant_mapping_2026_08_28 by removing the (vid, person_id) pairs in the exclusions table. Predicate validated against the delivered VDS in Step 4; exclusions built in Step 5. The 5,776 GvsMapUnmappedVIDs/GvsMapDroppedDuplicateVIDs fixup VIDs are NOT corrected here -- see Step 6a. `delivered` and `n_removed` are audit columns, dropped in Step 6.8; while they are present this table has not been validated.")
AS

WITH
-- THE FOUR CTEs THAT USED TO STAND HERE -- `fixup_vids`, `ex_pairs`, `light_ex`, `heavy_ex`
-- and `ex_vids` -- ARE NOW THE `_t` TABLES BUILT BY STEP 6.2b. Run 1 of the rehearsal showed
-- BigQuery would not broadcast a build side it could only reach through a `SELECT DISTINCT`
-- over 2.34 billion rows, and repartitioned all of `person_ids` instead. The definitions
-- moved verbatim; the reasoning behind each one moved with them and is worth reading there,
-- particularly why the fixup carve-out is a LEFT JOIN and not `NOT IN` (a single NULL from
-- the subquery would silently empty it and produce an output identical to the input, which
-- passes every validation below except STEP 6.4's reconciliation), and why `DISTINCT` is
-- what makes the 932 VS-2010 duplicate pairs inert rather than merely harmless.
--
-- ONE THING DID NOT MOVE. STEP 6.2b's `fixup_vids` omits the commented-out MANUAL block for
-- the two case-4 VIDs restored by hand after the scrub
-- ("17-72146919-GAGAAGAAGAAGAAGAAGAAGAAGAAGAAGA-G" and "1-143186828-T-TG"). They are
-- VERIFIED INERT as of 2026-09-05 -- P1 check (a) returned zero rows, so `exclusions` cannot
-- reach them. If that check is ever re-run and returns rows, they must be added to STEP
-- 6.2b's `fixup_vids` and `ex_pairs_t` rebuilt, not patched here.

-- BRANCH 1 -- 99.66% of VIDs. Nothing to remove, so nothing is computed: the arrays are
-- copied through untouched and `n_removed` is a literal.
unaffected AS (
  SELECT
    m.vid                          AS vid,
    m.person_ids                   AS person_ids,
    ARRAY_LENGTH(m.person_ids)     AS delivered,
    0                              AS n_removed
  FROM `foxtrot.vid_to_participant_mapping_2026_08_28` AS m
  LEFT JOIN `foxtrot.ex_vids_t` AS v ON v.vid = m.vid
  WHERE v.vid IS NULL
),

-- BRANCH 2 -- the light affected VIDs, IN TWO JOINS RATHER THAN ONE. `p NOT IN UNNEST(x.ids)`
-- is the membership test, row-local and quadratic, which the <= 1000 bound is what makes
-- acceptable.
--
-- THE SPLIT IS ABOUT WHICH TABLE MEETS `person_ids` FIRST, and it is the whole reason
-- `light_vids_t` exists. `light_ex_t` measured 1,313.1 MiB on 2026-09-05 -- 154.6M excluded
-- entries plus vid strings -- so BigQuery cannot broadcast it, and a direct
-- `mapping JOIN light_ex_t` forces the probe side to be repartitioned: all 21 TB of
-- `person_ids`, to find the 0.34% that match. `light_vids_t` is the same key set at ~133 MiB
-- with the arrays projected away, so `light_rows` reduces the mapping table to the affected
-- rows against a broadcastable build side, and only then does `light` fetch the arrays for
-- them. The second join's probe side is ~54,000 rows in the rehearsal, where a shuffle is
-- free in practice.
--
-- KEEP THE TWO TABLES CONSISTENT. `light_vids_t` and `light_ex_t` are built by the same
-- `HAVING COUNT(*) <= 1000` over the same `ex_pairs_t`, so their key sets are identical by
-- construction and the inner join below cannot drop a row. If one is ever rebuilt without
-- the other, this branch silently loses rows and only STEP 6.4's reconciliation notices.
light_rows AS (
  SELECT
    m.vid                          AS vid,
    m.person_ids                   AS person_ids,
    ARRAY_LENGTH(m.person_ids)     AS delivered
  FROM `foxtrot.vid_to_participant_mapping_2026_08_28` AS m
  JOIN `foxtrot.light_vids_t` AS k ON k.vid = m.vid
),

light AS (
  SELECT
    r.vid                                                AS vid,
    ARRAY(SELECT p FROM UNNEST(r.person_ids) AS p
          WHERE p NOT IN UNNEST(x.ids))                  AS person_ids,
    r.delivered                                          AS delivered
  FROM light_rows AS r
  JOIN `foxtrot.light_ex_t` AS x ON x.vid = r.vid
),

-- BRANCH 3 -- the heavy affected VIDs, one row per delivered person rather than per VID.
-- `delivered` is carried down the explosion so the reaggregation can report it without a
-- second read of `person_ids`.
heavy_exploded AS (
  SELECT
    m.vid                          AS vid,
    ARRAY_LENGTH(m.person_ids)     AS delivered,
    p                              AS person_id
  FROM `foxtrot.vid_to_participant_mapping_2026_08_28` AS m
  JOIN `foxtrot.heavy_ex_t` AS h ON h.vid = m.vid
  CROSS JOIN UNNEST(m.person_ids) AS p
),

-- The anti-join is a hash probe against `ex_pairs`, so it is linear in the exploded row
-- count and spread across slots, where the array test would have been quadratic on one.
--
-- A heavy VID whose every member is excluded contributes no group and so no output row,
-- which is exactly branch 2's `ARRAY_LENGTH > 0` filter arriving by a different route --
-- see the note on drop semantics below. `ANY_VALUE(delivered)` is safe because `delivered`
-- is constant within a `vid`, given P3 confirms `vid` is unique.
heavy AS (
  SELECT
    x.vid                          AS vid,
    ARRAY_AGG(x.person_id)         AS person_ids,
    ANY_VALUE(x.delivered)         AS delivered
  FROM heavy_exploded AS x
  LEFT JOIN `foxtrot.ex_pairs_t` AS e
    ON e.vid = x.vid AND e.person_id = x.person_id
  WHERE e.vid IS NULL
  GROUP BY x.vid
)

SELECT vid, person_ids, delivered, n_removed FROM unaffected
UNION ALL
SELECT
  vid,
  person_ids,
  delivered,
  delivered - ARRAY_LENGTH(person_ids)   AS n_removed
FROM light
-- Reproduces the pipeline's drop semantics: a VID with no surviving carrier gets no row at
-- all, not a row with an empty array. That is what the fixed workflow's join produces for
-- such a VID, and what the scrub did. Branch 3 gets the same semantics for free, since a
-- fully excluded VID leaves `heavy` with nothing to group.
--
-- This SHOULD remove nothing. Every VAT VID has `gvs_all_ac > 0` -- confirmed, zero of
-- 1,601,242,198 have zero -- so every VID has at least one FT-passing non-GQ-0 carrier in
-- the VDS, and if the mapping table names that carrier the array cannot empty. A nonzero
-- count here is therefore a FINDING, not an outcome: it means the mapping table named only
-- people who are not the VDS carrier, which is the case-4 representation mismatch
-- documented at the bottom of `scrub_participant_mapping_table.sql`. The scrub hit exactly
-- two such VIDs. STEP 6.5 counts them for free and lists them for cents.
--
-- P2 SHOWED HOW LITTLE MARGIN THAT ARGUMENT HAS AT THE TOP OF THE DISTRIBUTION. The worst
-- VID excludes 535,659 people out of a cohort of 535,662, so `delivered - n_removed` for it
-- is somewhere between 0 and 3, and the top 20 are all within ~3,100 of the same ceiling.
-- Those VIDs are the sharpest available test of the `gvs_all_ac > 0` argument: if the
-- predicate is right they keep a handful of carriers each, and if any of them empties, that
-- is a real case-4 finding rather than a rounding artifact. STEP 6.5 will name them.
WHERE ARRAY_LENGTH(person_ids) > 0
UNION ALL
SELECT
  vid,
  person_ids,
  delivered,
  delivered - ARRAY_LENGTH(person_ids)   AS n_removed
FROM heavy;

-- ONE THING THE ARRAY FILTER DOES THAT SUBTRACTION WOULD NOT, stated so it is not "fixed"
-- later. If a delivered array contains the same excluded person twice -- which happens for
-- the ~932 VS-2010 pairs, since the duplicated `alt_allele` rows both fed the original
-- `ARRAY_AGG` -- every copy is removed and `n_removed` counts 2. That is correct, and it is
-- why STEP 6.4's reconciliation expects a small positive excess rather than an exact match.


-- =====================================================================================
-- STEP 6.4 -- FREE validation, from the audit columns. ~$0.30. All rows must read 'pass'.
-- =====================================================================================
--
-- These do not reference `person_ids`, so they cost the two INT64 columns and nothing else.
-- Between them they are the strongest evidence available short of STEP 6.6, because
-- `SUM(n_removed)` has an exact expected value.
--
-- THE RECONCILIATION IS THE ONE TO READ FIRST. `exclusions` is a subset of the delivered
-- pairs by construction -- both come from `alt_allele` joined to the same active
-- non-control sample population and the same VAT alleles, with `exclusions` additionally
-- requiring GQ 0 or an `FT` failure -- so every excluded pair must have been present to
-- remove, and:
--
--   SUM(n_removed)  ==  distinct excluded pairs, for non-fixup VIDs present in the
--                       delivered table, PLUS the delivered arrays' duplicate multiplicity
--
-- ALL THREE TERMS ARE NOW MEASURED as of 2026-09-05, so this is no longer a loose
-- expectation but a numeric prediction with a 932-wide window on 2.34 billion:
--
--   distinct excluded pairs                     2,344,755,855  (P2)
--   less pairs on fixup VIDs                       -   37,702  (P3 follow-up)
--   plus delivered duplicate multiplicity          0 to   932  (P2 `surplus_rows`)
--                                               -------------
--   SUM(n_removed) expected             2,344,718,153 to 2,344,719,085
--
-- The floor holds if no delivered array carries a duplicate copy of an excluded person; the
-- ceiling if all 932 VS-2010 surplus pairs do. That is a tolerance of 0.00004%, which is why
-- the `BETWEEN 0 AND 2000` gate below is generous rather than lax -- anything outside this
-- window by more than a rounding artifact is a real defect.
--
-- "present in the delivered table" is VACUOUS: P3 measured `orphan_exclusion_vids` = 0, so
-- every excluded pair has a row to be removed from. The fixup subtraction is 37,702 pairs
-- over 121 VIDs, a mean of 311.6 -- BELOW the global 431.2. The concern that motivated
-- measuring rather than assuming (a dropped-duplicate VID is exactly the kind that can hold
-- a very long array) did not materialize: these are ordinary-weight VIDs, and the term is
-- 0.0016% of the total.
--
-- A SHORTFALL means excluded pairs were not found in the arrays: the join is wrong, the
-- carve-out removed too much, or the two tables were built from different populations.
-- Investigate before promoting the table in STEP 6.8; this is the sharpest signal in the file.
--
-- An EXCESS is expected and should be small -- bounded by P2's `surplus_rows`, 932. An
-- excess materially larger than that means the delivered arrays carry duplicate person ids
-- from some source other than VS-2010, which is worth understanding but does not block:
-- removing every copy of an excluded person is right either way.
--
-- THIS DELIBERATELY DOES NOT READ `foxtrot.ex_pairs_t`, AND MUST NOT BE "SIMPLIFIED" TO.
-- STEP 6.2b's table is exactly the input STEP 6.3 consumed, so validating one against the
-- other would be an identity, not a check: if `ex_pairs_t` is wrong, both sides are wrong
-- together and this passes. Re-deriving from `exclusions` costs ~$0.30 and is the only
-- reason the reconciliation can catch a bad carve-out at all. The same applies to STEP 6.6.
-- What `ex_pairs_t` IS good for here is a free prior: `SELECT COUNT(*)` on it is metadata
-- only and must read 2,344,718,153, so check that before trusting anything downstream.

WITH
fixup_vids AS (
  SELECT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`      GROUP BY vid
  UNION DISTINCT
  SELECT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01` GROUP BY vid
  -- MANUAL: the two case-4 VIDs restored by hand after the scrub. VERIFIED INERT 2026-09-05
  -- -- P1 check (a) returned zero rows, so `exclusions` cannot reach them. Leave commented;
  -- uncomment BOTH lines only if that check is ever re-run and returns rows.
  -- UNION DISTINCT
  -- SELECT vid FROM UNNEST(["17-72146919-GAGAAGAAGAAGAAGAAGAAGAAGAAGAAGA-G",
  --                         "1-143186828-T-TG"]) AS vid
),
expected AS (
  -- `SELECT DISTINCT vid` on the mapping side, not `SELECT vid`: a duplicate vid row would
  -- otherwise fan this out and inflate the expectation. P3 says there are none; do not rely
  -- on a check running before this one.
  SELECT COUNT(*) AS n
  FROM (SELECT DISTINCT vid, person_id FROM `foxtrot.exclusions`)  AS e
  JOIN (SELECT DISTINCT vid FROM `foxtrot.vid_to_participant_mapping_2026_08_28`)       AS m USING (vid)
  LEFT JOIN fixup_vids AS f USING (vid)
  WHERE f.vid IS NULL
),
dropped AS (
  -- THE 172 EMPTIED VIDs ARE NOT IN THE OUTPUT, SO THEIR REMOVALS ARE NOT IN SUM(n_removed).
  -- `expected` below counts pairs for every non-fixup vid in the DELIVERED table, which
  -- includes them, so without this term the two populations differ and the reconciliation
  -- goes NEGATIVE on a correct build -- which is exactly what happened on 2026-09-05 before
  -- this was added. Built by STEP 6.5; see the note there.
  SELECT IFNULL(SUM(excluded_pairs), 0) AS n FROM `foxtrot.emptied_vids_t`
),
actual AS (
  SELECT
    SUM(delivered)                             AS entries_delivered,
    SUM(n_removed)                             AS entries_removed,
    SUM(delivered) - SUM(n_removed)            AS entries_corrected,
    COUNTIF(n_removed > 0)                     AS vids_changed,
    COUNTIF(n_removed > delivered)             AS vids_negative,
    COUNTIF(delivered = 0)                     AS vids_delivered_empty,
    MAX(n_removed)                             AS worst_vid_removals,
    APPROX_QUANTILES(IF(n_removed > 0, SAFE_DIVIDE(delivered, delivered - n_removed), NULL), 100)[OFFSET(50)] AS median_inflation_affected,
    APPROX_QUANTILES(IF(n_removed > 0, SAFE_DIVIDE(delivered, delivered - n_removed), NULL), 100)[OFFSET(99)] AS p99_inflation_affected,
    MAX(IF(n_removed > 0, SAFE_DIVIDE(delivered, delivered - n_removed), NULL))                               AS max_inflation
  FROM `foxtrot.vid_to_participant_mapping_2026_09_05`
)
SELECT "corrected count is never negative"        AS check_name,
       IF(vids_negative = 0, "pass", "FAIL")      AS result,
       CAST(vids_negative AS STRING)              AS value
FROM actual
UNION ALL
SELECT "no delivered array was already empty",
       IF(vids_delivered_empty = 0, "pass", "FAIL"),
       CAST(vids_delivered_empty AS STRING) FROM actual
UNION ALL
SELECT "removals reconcile against exclusions (0 <= excess <= ~932)",
       IF(((SELECT entries_removed FROM actual) + (SELECT n FROM dropped)
           - (SELECT n FROM expected)) BETWEEN 0 AND 2000, "pass", "FAIL"),
       CONCAT(CAST((SELECT entries_removed FROM actual) AS STRING), " removed + ",
              CAST((SELECT n FROM dropped) AS STRING), " on emptied vids vs ",
              CAST((SELECT n FROM expected) AS STRING), " expected, excess ",
              CAST((SELECT entries_removed FROM actual) + (SELECT n FROM dropped)
                   - (SELECT n FROM expected) AS STRING))
UNION ALL
SELECT "-- report --", "", ""
UNION ALL
SELECT "entries delivered",  "", CAST((SELECT entries_delivered FROM actual) AS STRING)
UNION ALL
SELECT "entries removed",    "", CAST((SELECT entries_removed FROM actual) AS STRING)
UNION ALL
SELECT "removals on emptied vids", "", CAST((SELECT n FROM dropped) AS STRING)
UNION ALL
SELECT "vids emptied and dropped", "", CAST((SELECT COUNT(*) FROM `foxtrot.emptied_vids_t`) AS STRING)
UNION ALL
SELECT "entries corrected",  "", CAST((SELECT entries_corrected FROM actual) AS STRING)
UNION ALL
SELECT "entry rate",         "", CAST(ROUND(100 * SAFE_DIVIDE((SELECT entries_removed FROM actual),
                                                              (SELECT entries_delivered FROM actual)), 4) AS STRING)
UNION ALL
SELECT "vids changed",       "", CAST((SELECT vids_changed FROM actual) AS STRING)
UNION ALL
SELECT "worst vid removals", "", CAST((SELECT worst_vid_removals FROM actual) AS STRING)
UNION ALL
SELECT "median inflation, affected vids", "", CAST((SELECT median_inflation_affected FROM actual) AS STRING)
UNION ALL
SELECT "p99 inflation, affected vids",    "", CAST((SELECT p99_inflation_affected FROM actual) AS STRING)
UNION ALL
SELECT "max inflation",                   "", CAST((SELECT max_inflation FROM actual) AS STRING);

-- `entries_delivered` is the exact denominator D2 was going to produce, so the provisional
-- genome-wide rate of ~0.092% becomes a measured one here, for free, on the artifact that
-- ships. That was the one number QUERY C's restructuring appeared to lose.

-- RESULT 2026-09-05 -- ALL THREE CHECKS PASS. THE GATE IS GREEN.
--
--     corrected count is never negative                        pass (0)
--     no delivered array was already empty                     pass (0)
--     removals reconcile against exclusions                    pass, EXCESS 930
--
--     entries delivered           2,559,382,818,198
--     entries removed             2,341,606,369     (predicted 2,341,605,439-2,341,606,371)
--     removals on emptied vids        3,112,714
--     entries corrected           2,557,041,211,829
--     entry rate                          0.0915 %
--     vids changed                    5,437,323     (0.34 % of 1,601,242,198)
--     vids emptied and dropped              172
--     worst vid removals                535,659
--     inflation on affected vids    median 1.0244, p99 7.5968, max 534,942
--
-- EXCESS 930 AGAINST A CEILING OF 932 IS THE STRONGEST CORROBORATION IN THIS WHOLE EFFORT,
-- because the two numbers were derived with nothing in common. 932 came from counting
-- duplicate `(vid, person_id)` pairs in `exclusions` for the VS-2010 scope update. 930 falls
-- out of an arithmetic identity over a 2.56-trillion-entry table that never looks at
-- duplicates at all. Do not chase the 2-pair gap: a duplicated pair raises the excess only if
-- BOTH copies actually sat in the delivered array, and several benign situations prevent that
-- -- one copy a GQ 0 no-call and the other passing, or the pair landing on a carved-out fixup
-- VID. Landing at 930 with no term for any of it is the expected shape, not a shortfall.
--
-- The 0.0915 % entry rate sits just above the 0.0825 % of `alt_allele` rows that `exclusions`
-- covers, which is the right direction: the mapping table is restricted to VAT alleles, where
-- exclusions concentrate.
--
-- WHAT 6.4 DOES NOT ESTABLISH. It proves the right NUMBER of entries went away, not that they
-- were the right ONES -- removing an equal count of wrong pairs would pass here. That
-- identity is checked at 1 % by STEP 6.2's V1a and at full scale only by STEP 6.6.


-- =====================================================================================
-- STEP 6.5 -- CHEAP validation and the emptied-VID list. vid columns only, ~$0.40.
-- =====================================================================================
--
-- `COUNT(*)` is answered from metadata, so the row-count delta is free. Per STEP 6.3's
-- closing note the expectation is zero, and anything else is a finding rather than an
-- outcome -- so identify them individually rather than reporting a number.
--
-- ONE VID TO CHECK BY NAME: `3-65490858-G-A`. P4 QUERY C measured `gvs_all_ac` = 1 there --
-- a SINGLE alt allele in the entire callset, so exactly one heterozygous participant, against
-- 534,941 named as carriers in the delivered table. That is the narrowest margin anywhere in
-- the top 20 and therefore the sharpest test of the no-empty argument: the VID survives if
-- and only if that one person is in its delivered array. It is a SNP, so a left-alignment
-- mismatch cannot explain a failure and an empty here would be a genuine case-4 finding.
-- 7-152451810-T-G and 13-37650636-T-C are next tightest at `AC` = 2.

SELECT
  (SELECT COUNT(*) FROM `foxtrot.vid_to_participant_mapping_2026_08_28`)           AS vid_rows_before,
  (SELECT COUNT(*) FROM `foxtrot.vid_to_participant_mapping_2026_09_05`) AS vid_rows_after,
  (SELECT COUNT(*) FROM `foxtrot.vid_to_participant_mapping_2026_08_28`)
    - (SELECT COUNT(*) FROM `foxtrot.vid_to_participant_mapping_2026_09_05`) AS vids_emptied;

-- The emptied VIDs themselves. Anti-join on `vid` alone; do NOT add `person_ids`.
SELECT o.vid
FROM `foxtrot.vid_to_participant_mapping_2026_08_28` AS o
LEFT JOIN `foxtrot.vid_to_participant_mapping_2026_09_05` AS c ON o.vid = c.vid
WHERE c.vid IS NULL;

-- For each one that comes back, the question to ask is the scrub's case 4: does the VDS
-- hold this variant under a different representation, carried by someone the mapping table
-- never named? `check_vds_for_alleles.py` answers it directly. Such a VID is not made worse
-- by being dropped here -- it was already wrong, naming a non-carrier -- but it is a
-- pointer at the unfixed half of that finding, and it belongs in the delivery note rather
-- than being silently absorbed.
--
-- COVERAGE, the other direction. Every active non-control sample should still have mapping
-- entries after the correction. Run this against the SCRUBBED SAMPLE output from STEP 6.2,
-- not the full table: the correction is row-local, so correcting a subset equals a subset
-- of correcting, which makes the sampled output a true subset of the full output. Zero
-- missing on the sample is therefore conclusive for the whole table, and only a NONZERO
-- result is ambiguous. Same asymmetry the scrub's STEP 4c relied on.

WITH present AS (
  SELECT DISTINCT p AS person_id
  FROM `foxtrot.vid_to_participant_mapping_2026_09_05_sample` AS m, UNNEST(m.person_ids) AS p
)
SELECT COUNT(*) AS active_samples_with_no_mapping_entry
FROM `foxtrot.sample_info` AS si
LEFT JOIN present ON present.person_id = SAFE_CAST(si.sample_name AS INT64)
WHERE si.withdrawn IS NULL
  AND si.is_control = false
  AND present.person_id IS NULL;


-- RESULT, 2026-09-05: 172 VIDs EMPTIED against an expectation of ZERO, and the finding is
-- BENIGN FOR THIS BUILD but is a real gap in the fixup workflows' coverage.
--
--   0 SNPs, 172 indels, 0 malformed vids, 0 absent from `ex_pairs_t`, 3,112,714 pairs total.
--
-- The SNP count is the whole discriminator. A SNP is already minimal and left-aligned, so an
-- empty there could not be blamed on representation and would be a genuine case-4 finding.
-- There are none. The length spectrum splits into two mechanisms:
--
--   single-base insertions      99 VIDs (57.6%)   3,092,468 pairs (99.35%)
--   insertions 2-3 bp            8 VIDs            9,810 pairs
--   deletions 197-314 bp        50 VIDs (29.1%)    2,028 pairs (0.07%)
--   deletions under 190 bp      15 VIDs            8,408 pairs
--
-- The 99 single-base insertions are homopolymer slippage -- an inserted base in a run can be
-- written anywhere in the run, so the left-aligned vid and the non-left-aligned `alt_allele`
-- row sit at different coordinates. They carry essentially all the affected entries. The 50
-- deletions at 197-314 bp across 35 distinct lengths are Alu-scale: a polymorphic Alu
-- insertion reads as a deletion against a reference carrying the element, in repetitive
-- context where normalization is hardest. Negligible by entries, a third of the VIDs.
--
-- NONE OF THE 172 CAN BE A FIXUP VID, and that follows without a query: fixup VIDs are carved
-- out of STEP 6.3 and take branch 1 untouched, so they cannot be emptied. So these are 172
-- case-4 VIDs that `GvsMapUnmappedVIDs` and `GvsMapDroppedDuplicateVIDs` did NOT catch. That
-- belongs to STEP 6a's scope and to a ticket about those workflows' coverage; it is not a
-- defect in this build.
--
-- WHY AN EMPTIED VID IS NOT A CONTRADICTION. The VID is in the VAT, so `gvs_all_ac` > 0, so
-- some carrier passed `FT` -- yet every NAMED carrier was excluded. Both hold if the passing
-- carriers sit at a different representation and were never in the delivered array, which is
-- what case 4 means. The magnitudes support it: the worst is 463,534 named carriers against a
-- 535,662 cohort, the hg38-carries-the-minor-allele class, where P4 showed the VAT entry is
-- carried by a handful of `1/2` genotypes rescued via `any_yes` on a partner allele.
--
-- CONFIRMED 2026-09-05 by reading the VAT for all 172. `gvs_all_ac` is small and nonzero
-- throughout -- mostly 1-5, largest 746 (`21-44410741-T-TG`), nothing above four figures --
-- and `ac_variants` is 1 on every VID, so the VAT's redundant per-transcript copies agree.
-- Neither failure mode fired: no zero AC (which would mean the VID does not belong in the
-- VAT), and no large AC (which would mean the delivered array was missing real carriers at
-- scale -- the one reading that would have indicted the build rather than the fixups).
--
-- `gvs_all_an` INVERTS HERE, AND THAT IS THE SHARPEST EVIDENCE FOR CASE 4. For the top-20
-- correctly-mapped extreme VIDs, `excluded + AN/2` sat at or just below the 535,662 cohort
-- every time, which is what ruled out over-exclusion there: a participant cannot both hold a
-- passing genotype and be excluded AT THE SAME SITE, so the two populations are disjoint and
-- the sum is bounded. Among these 172 the bound BREAKS -- on 3 of the top 20 by entry count:
--
--     18-57678586-C-CA    AN/2 500,531 + excluded 105,340 = 605,871   +70,209 over cohort
--     12-124467987-A-AC   AN/2 107,559 + excluded 461,748 = 569,307   +33,645 over cohort
--     4-76917565-A-AG     AN/2  72,297 + excluded 463,534 = 535,831      +169 over cohort
--
-- Under any reading where the mapping table and the VAT describe the same site that is a
-- contradiction. Under case 4 it is the expected result: `excluded` was computed from
-- `alt_allele` at the vid's coordinates, `AN` from the VDS at the vid's true left-aligned
-- representation, and if those are different sites the two counts are not competing for the
-- same 535,662 people, so overlap is unremarkable. Same test, opposite outcome, and the
-- outcome tracks exactly which VIDs are mismapped -- a sharper discriminator than the AC
-- magnitude. `4-76917565-A-AG` is the cleanest instance: AC 12 against AN 144,594, so 72,297
-- participants hold a passing call while all 463,534 NAMED carriers were excluded, which is
-- coherent only if the passing carriers sit at a representation the delivered array never
-- named.
--
-- AGGREGATED OVER ALL 172, 2026-09-05:
--
--     vids_in_vat 172   missing_from_vat 0   ac_zero 0
--     min_ac 1   median_ac 2   max_ac 746
--     over_cohort 3   worst_excess 70,209
--
-- `missing_from_vat` = 0 and `ac_zero` = 0 close the last inference: all 172 belong in the
-- VAT and none is a phantom. `median_ac` 2 against `max_ac` 746 puts the distribution even
-- further down than the top-30 view suggested.
--
-- READ `over_cohort` = 3 AS A FLOOR, NOT A COUNT. It matches the top-20 hand check exactly,
-- so no violation exists outside the top 20 -- but the bound can only break when BOTH the
-- excluded count and `AN/2` are large, and where either is small the sum stays under 535,662
-- for arithmetic reasons unrelated to whether the VID is mapped correctly. The test is
-- one-sided and mostly powerless across these 172: breaking it proves mismapping, not
-- breaking it proves nothing. The count of mismapped VIDs is 172, established by the 0-SNP
-- argument above, which does not depend on this at all.

-- =====================================================================================
-- STEP 6.6 -- FULL validation. 1 x T, ~$119. Optional, and recommended.
-- =====================================================================================
--
-- Reads the corrected copy only -- not both copies -- which is what keeps it at 1 x T. It
-- answers the half of the question STEP 6.4 cannot: that no excluded pair SURVIVED. STEP
-- 6.4 establishes that the right NUMBER of entries went away; this establishes that they
-- were the right ones.
--
-- The other half -- that nothing else went away -- stays at sample fidelity (STEP 6.2's
-- V2b), because checking it at full scale requires reading both copies at 2 x T. That
-- asymmetry is deliberate and worth stating in the delivery note rather than glossing:
-- over-removal is bounded by the sample rehearsal and by STEP 6.4's reconciliation, not
-- proven exhaustively.
--
-- The exclusions side is restricted to non-fixup VIDs to match what STEP 6.3 applied.
-- Without that, the 1,542 dropped-duplicate VIDs deliberately left uncorrected would each
-- report as a surviving excluded pair and the check would FAIL by design.

WITH fixup_vids AS (
  SELECT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`      GROUP BY vid
  UNION DISTINCT
  SELECT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01` GROUP BY vid
  -- MANUAL: the two case-4 VIDs restored by hand after the scrub. VERIFIED INERT 2026-09-05
  -- -- P1 check (a) returned zero rows, so `exclusions` cannot reach them. Leave commented;
  -- uncomment BOTH lines only if that check is ever re-run and returns rows.
  -- UNION DISTINCT
  -- SELECT vid FROM UNNEST(["17-72146919-GAGAAGAAGAAGAAGAAGAAGAAGAAGAAGA-G",
  --                         "1-143186828-T-TG"]) AS vid
),
applicable AS (
  SELECT DISTINCT e.vid, e.person_id
  FROM `foxtrot.exclusions` AS e
  LEFT JOIN fixup_vids AS f ON f.vid = e.vid
  WHERE f.vid IS NULL
)
SELECT "no applicable excluded pair survives" AS check_name,
       IF(NOT EXISTS(
            SELECT 1
            FROM `foxtrot.vid_to_participant_mapping_2026_09_05` AS c
            CROSS JOIN UNNEST(c.person_ids) AS p
            JOIN applicable AS a ON a.vid = c.vid AND a.person_id = p
          ), "pass", "FAIL") AS result
UNION ALL
SELECT "no empty person_ids array",
       IF(NOT EXISTS(
            SELECT 1 FROM `foxtrot.vid_to_participant_mapping_2026_09_05`
            WHERE ARRAY_LENGTH(person_ids) = 0
          ), "pass", "FAIL")
UNION ALL
SELECT "n_removed agrees with the arrays",
       IF(NOT EXISTS(
            SELECT 1 FROM `foxtrot.vid_to_participant_mapping_2026_09_05`
            WHERE n_removed != delivered - ARRAY_LENGTH(person_ids)
          ), "pass", "FAIL");

-- RESULT 2026-09-05 -- ALL THREE PASS. THE TABLE IS FULLY VALIDATED.
--
--     no applicable excluded pair survives   pass
--     no empty person_ids array              pass
--     n_removed agrees with the arrays       pass
--
-- Together with 6.4 this closes the build: 6.4 established that the right NUMBER of entries
-- went away, 6.6 that at full scale they were the right ONES. The residual over-reporting is
-- now exactly the 121 fixup VIDs deliberately carved out (order 1,800 pairs, 7.7e-07 of the
-- correction) plus the 172 case-4 VIDs found by 6.5 -- both STEP 6a's scope, both gating
-- STEP 8, neither a defect here.
--
-- 6.6R WAS NOT NEEDED. Job job_IaTHh288aJ7RGG3FQmZQUF9fWCW2 took a 32-stage plan against
-- 6.3's 54 and was 24 of 32 stages done at 24.4 min on 47,900 slot-min at 1,962 slots, so
-- BigQuery co-located on `vid` and never materialized the 2.56e12 entries the flat form
-- risks. The three `Output` stages are the three-way UNION ALL: checks 2 and 3 do not unnest
-- and retire early. Keep 6.6R anyway -- the plan choice was BigQuery's, not the query's, and
-- is not guaranteed to repeat on a differently-shaped table.


-- =====================================================================================
-- STEP 6.6R -- semi-join reducer form of 6.6. Same contract, different shape.
-- =====================================================================================
--
-- USE THIS INSTEAD OF 6.6 IF 6.6 RUNS LONG OR FAILS ON RESOURCES. It is not a correction --
-- 6.6 is right -- it is a rewrite of the same predicate to keep an unbounded term out of the
-- shuffle.
--
-- WHY. 6.6's first check joins `applicable` to the CROSS JOIN UNNEST of EVERY row in the
-- corrected table. That unnest side is 2,557,041,211,829 entries. `applicable` is ~2.34e9
-- pairs at ~34 B/record, so ~80 GB -- far above the measured broadcast threshold of somewhere
-- under 133 MiB -- which forces a hash join. If BigQuery partitions on the full
-- (vid, person_id) key rather than co-locating on `vid` alone, all 2.56 trillion entries
-- enter the shuffle. STEP 6.3 never did this: its stage table shows S00 reading `ex_pairs_t`
-- at 34 B/record while the four expensive stages ran 1.6e9 records at 14,417 B/record, i.e.
-- narrow pairs and INTACT arrays joined on `vid`. And it still spilled 56 TiB of an 86.2 TiB
-- shuffle. There is no memory headroom to absorb a heavier shape.
--
-- WHAT THE REDUCER BUYS, AND WHAT IT DOES NOT. `touched_vids` is ~5,437,495 VIDs (6.4's
-- 5,437,323 changed plus the 172 emptied). Affected VIDs carry roughly 9.8e10 delivered
-- entries -- 2.34e9 removed at 6.4's median inflation of 1.0244 -- so the unnest drops from
-- 2.56e12 to ~1e11, about 26x. It does NOT cut the bill. 5.4M VIDs is 0.34 % of
-- 1,601,242,026 rows scattered across the whole `vid` keyspace, so essentially every cluster
-- block holds at least one and CLUSTER BY vid (line 1946) prunes nothing. Expect ~18.7 TiB
-- and ~$117 either way; this is about whether it FINISHES, not what it costs.
--
-- Checks 2 and 3 are unchanged and never unnest, so they are left exactly as in 6.6.

WITH fixup_vids AS (
  SELECT vid FROM `foxtrot.unmapped_vid_mapping_2026_06_01`      GROUP BY vid
  UNION DISTINCT
  SELECT vid FROM `foxtrot.dropped_duplicate_mapping_2026_06_01` GROUP BY vid
),
applicable AS (
  SELECT DISTINCT e.vid, e.person_id
  FROM `foxtrot.exclusions` AS e
  LEFT JOIN fixup_vids AS f ON f.vid = e.vid
  WHERE f.vid IS NULL
),
touched_vids AS (
  SELECT vid FROM applicable GROUP BY vid
),
candidates AS (
  -- The ONLY rows that can carry a surviving excluded pair. Every other row is unreachable
  -- by `applicable` by construction, so declining to unnest it loses no coverage.
  SELECT c.vid, c.person_ids
  FROM `foxtrot.vid_to_participant_mapping_2026_09_05` AS c
  JOIN touched_vids AS t ON t.vid = c.vid
)
SELECT "no applicable excluded pair survives" AS check_name,
       IF(NOT EXISTS(
            SELECT 1
            FROM candidates AS c
            CROSS JOIN UNNEST(c.person_ids) AS p
            JOIN applicable AS a ON a.vid = c.vid AND a.person_id = p
          ), "pass", "FAIL") AS result
UNION ALL
SELECT "no empty person_ids array",
       IF(NOT EXISTS(
            SELECT 1 FROM `foxtrot.vid_to_participant_mapping_2026_09_05`
            WHERE ARRAY_LENGTH(person_ids) = 0
          ), "pass", "FAIL")
UNION ALL
SELECT "n_removed agrees with the arrays",
       IF(NOT EXISTS(
            SELECT 1 FROM `foxtrot.vid_to_participant_mapping_2026_09_05`
            WHERE n_removed != delivered - ARRAY_LENGTH(person_ids)
          ), "pass", "FAIL");


-- =====================================================================================
-- STEP 6.7 -- preserve the audit columns before dropping them. ~$0.30.
-- =====================================================================================
--
-- Reads `vid` and the two INT64 columns, never `person_ids`, so this is cents rather than
-- $119. It is the per-VID impact artifact the plan wanted from D2 -- "how wrong was my
-- cohort", answerable per VID forever after -- and it must be written BEFORE STEP 6.8,
-- which makes the columns unreadable.

-- THE 172 EMPTIED VIDs MUST BE UNIONED IN, AND THIS IS EASY TO MISS. They are NOT rows of
-- `vid_to_participant_mapping_2026_09_05` -- STEP 6.3 dropped them -- so `WHERE n_removed > 0`
-- against that table alone silently omits precisely the VIDs with the LARGEST impact: 100 %
-- removal, up to 463,534 entries. An artifact described as per-VID impact that drops its own
-- worst cases is worse than one that never claimed to be complete, so `fully_emptied` is a
-- column rather than an omission. For those rows `delivered` is taken from `excluded_pairs`,
-- which is exact: the array was emptied entirely, so every delivered entry matched an
-- exclusion pair, and the two can differ only by a VS-2010 duplicate (at most 932 genome-wide
-- across 924 VIDs, so ~0 expected among these 172). `inflation` is NULL there, not infinity.

CREATE OR REPLACE TABLE `foxtrot.mapping_correction_audit`
CLUSTER BY vid
OPTIONS (description = "Step 6: per-VID impact of the participant mapping correction, i.e. the delta from vid_to_participant_mapping_2026_08_28 to vid_to_participant_mapping_2026_09_05. delivered = entries in the table as shipped after the withdrawn/control scrub; removed = GQ 0 and FT-failing entries taken out; corrected = what ships now. fully_emptied marks the 172 VIDs whose every named carrier was excluded; those rows were DROPPED from the corrected mapping table rather than left with an empty array, their delivered is taken from the exclusion count, and their inflation is NULL rather than infinite. They are case-4 representation mismatches missed by GvsMapUnmappedVIDs/GvsMapDroppedDuplicateVIDs -- see Step 6a.")
AS
SELECT
  vid,
  delivered,
  n_removed                                          AS removed,
  delivered - n_removed                              AS corrected,
  SAFE_DIVIDE(delivered, delivered - n_removed)      AS inflation,
  FALSE                                              AS fully_emptied
FROM `foxtrot.vid_to_participant_mapping_2026_09_05`
WHERE n_removed > 0
UNION ALL
SELECT
  vid,
  excluded_pairs                                     AS delivered,
  excluded_pairs                                     AS removed,
  0                                                  AS corrected,
  CAST(NULL AS FLOAT64)                              AS inflation,
  TRUE                                               AS fully_emptied
FROM `foxtrot.emptied_vids_t`;

-- Restricted to changed VIDs: 5.4 million rows rather than 1.6 billion, and a VID with
-- `removed = 0` carries no information the delivered table does not already carry.


-- =====================================================================================
-- STEP 6.8 -- drop the audit columns, which promotes the table. DDL metadata; free.
-- =====================================================================================
--
-- The delivered table must have the schema the pipeline produces -- `vid` and `person_ids`,
-- nothing else -- so that the next callset's output is comparable and so that Step 8's
-- Parquet export does not depend on an explicit column list to stay correct.
--
-- Because there is no rename, THIS is the step that promotes `2026_09_05` from a candidate
-- to the current version of the table, and the audit columns are what mark the difference.
-- A `vid_to_participant_mapping_*` table carrying `delivered` and `n_removed` has not
-- passed STEP 6.4 through STEP 6.6 yet; one carrying just `vid` and `person_ids` has. Do
-- not run this until those are clean -- it is the only signal anyone looking at the dataset
-- later will have.
--
-- `ALTER TABLE DROP COLUMN` is a metadata operation and bills nothing. Storage for the
-- dropped columns is reclaimed asynchronously rather than immediately, which is fine at 26
-- GB. If BigQuery rejects it -- the usual cause is recent streaming activity on the table,
-- which does not apply here -- the fallback is
-- `CREATE OR REPLACE TABLE ... AS SELECT vid, person_ids FROM ...`, which costs another
-- 1 x T. Confirm the DDL succeeded with the P1 query before proceeding.
--
-- ***** HOLD 2026-09-05 -- DO NOT RUN YET. *****
--
-- Deferred at the user's request: collaborators are expected to ask about specific VIDs, and
-- a one-table clustered lookup carrying `delivered` and `n_removed` is easier to run in front
-- of someone than a join. STEP 8's Parquet export is the only hard deadline -- the columns
-- must be gone before it, or they change the delivered schema.
--
-- THAT DEADLINE IS DISCHARGED BY THE EXPORT MECHANISM, NOT BY DROPPING THESE COLUMNS.
-- STEP 8 is to use `EXPORT DATA OPTIONS (uri = ..., format = 'PARQUET', compression =
-- 'SNAPPY', overwrite = TRUE) AS SELECT vid, person_ids FROM ...`, which emits the two-column
-- artifact while leaving this table's schema untouched. It costs ~$117 -- it is a query, not
-- a free `bq extract`, and `person_ids` is essentially the whole 18.65 TiB.
--
-- DO NOT "OPTIMIZE" THIS BY DROPPING THE COLUMNS FIRST AND USING THE FREE EXTRACT. That was
-- the initial recommendation here and the timing in it is backwards: collaborator questions
-- FOLLOW the delivery -- they cannot ask about mappings they have not received -- so
-- "drop immediately before the export" is the same instant as "drop immediately before the
-- questions start," which is exactly when the columns are wanted.
--
-- The purchase is convenience, not data preservation, and the distinction matters if anyone
-- reopens this: `mapping_correction_audit` and `exclusions` both outlive the drop, and
-- `explain_vid_for_collaborators.sql` reads only those. What ~$117 buys is that a single
-- clustered lookup against the delivered table keeps working while someone is waiting.
--
-- See the Step 8 section of participant-mapping-over-reporting-plan.md.
--
-- THE "AUDIT COLUMNS MEAN UNVALIDATED" SIGNAL ABOVE IS NOW WRONG, and that is the cost of
-- holding. 6.4 and 6.6 both went green on 2026-09-05 while the columns are still present, so
-- anyone reading the dataset cold would draw the opposite conclusion from the one intended.
-- The table's OPTIONS description was rewritten the same day to say so explicitly; that text,
-- not the schema, is now the signal. Re-read it before relying on either.
--
-- KEEPING THEM IS A CONVENIENCE, NOT A DEPENDENCY. Nothing becomes unrecoverable at 6.8. For
-- a changed VID `mapping_correction_audit` carries `delivered`, `removed` and `corrected`;
-- for an unchanged VID `removed` is 0 and `delivered` is `ARRAY_LENGTH(person_ids)` by
-- construction. `explain_vid_for_collaborators.sql` is written against those two facts and
-- needs neither column, so it keeps working after the drop.

-- Run FIRST, so the description stops asserting something false while the columns remain.
ALTER TABLE `foxtrot.vid_to_participant_mapping_2026_09_05`
SET OPTIONS (description = "Step 6: Foxtrot participant mapping with GQ 0 no-calls and FT-failing genotypes removed. Derived from vid_to_participant_mapping_2026_08_28 by removing the (vid, person_id) pairs in the exclusions table. Predicate validated against the delivered VDS in Step 4; exclusions built in Step 5. VALIDATED 2026-09-05: Step 6.4 reconciled removals to within 930 of the 932-pair VS-2010 ceiling, and Step 6.6 confirmed at full scale that no applicable excluded pair survives. The `delivered` and `n_removed` audit columns are RETAINED for now to support per-VID questions from collaborators, and must be dropped before the Step 8 Parquet export; their presence no longer indicates an unvalidated table. The 5,776 GvsMapUnmappedVIDs/GvsMapDroppedDuplicateVIDs fixup VIDs are NOT corrected here, and 172 further case-4 VIDs were found emptied and dropped -- see Step 6a and `emptied_vids_t`.");

-- HELD -- see above. Run both, then re-run the P1 query to confirm the DDL took.
-- ALTER TABLE `foxtrot.vid_to_participant_mapping_2026_09_05` DROP COLUMN delivered;
-- ALTER TABLE `foxtrot.vid_to_participant_mapping_2026_09_05` DROP COLUMN n_removed;


-- =====================================================================================
-- STEP 6.9 -- housekeeping. Free.
-- =====================================================================================
--
-- There is nothing to swap. `vid_to_participant_mapping_2026_09_05` is the current version
-- the moment STEP 6.8 promotes it, and `vid_to_participant_mapping_2026_08_28` is left
-- exactly as delivered -- which makes it the rollback target for free, with no rename to
-- reverse and no window in which the canonical name points at an unvalidated table.
--
-- Run only after STEP 6.4 and STEP 6.5 are clean, STEP 6.2's rehearsal passed, and STEP 6.6
-- if you are running it.
--
-- The rehearsal tables and STEP 6.2b's scratch tables are pure throwaway; drop them.
-- Drop the `_t` tables LAST and only after STEP 6.6, since STEP 6.4's free `COUNT(*)`
-- cross-check on `ex_pairs_t` is the cheapest evidence available that the carve-out was
-- right, and it is unavailable once the table is gone.

DROP TABLE IF EXISTS `foxtrot.vid_to_participant_mapping_2026_08_28_sample`;
DROP TABLE IF EXISTS `foxtrot.vid_to_participant_mapping_2026_09_05_sample`;
DROP TABLE IF EXISTS `foxtrot.ex_vids_t`;
DROP TABLE IF EXISTS `foxtrot.light_vids_t`;
DROP TABLE IF EXISTS `foxtrot.light_ex_t`;
DROP TABLE IF EXISTS `foxtrot.heavy_ex_t`;
DROP TABLE IF EXISTS `foxtrot.ex_pairs_t`;

-- `emptied_vids_t` IS DELIBERATELY NOT DROPPED HERE. It is 172 rows, it costs nothing, and
-- it is the only surviving record of which VIDs left the table -- reconstructing it after
-- `ex_pairs_t` is gone is impossible. STEP 6.4 reads it. Fold it into STEP 6.7's audit side
-- table, and drop it only once that is in place.

-- NOTHING ELSE IS DROPPED. The `vid_to_participant_mapping%` family totals 57.07 TiB and
-- ~$1,168/month and this step adds ~$382, but the collaborators have never been sensitive to
-- BigQuery storage costs and what they want is the corrected mappings. `2026_08_28` is the
-- rollback target and stays; the older versions are harmless where they are. Revisit only if
-- someone asks.


-- =====================================================================================
-- WHAT REMAINS BEFORE STEP 8
-- =====================================================================================
--
-- **Step 6a is not optional.** The 5,776 fixup VIDs are carved out of everything above and
-- are still wrong. They are 0.0004% of the table, but they are the pathological VIDs the
-- fixup workflows exist for, and shipping a "corrected" table that quietly skips them is
-- worse than shipping the uncorrected one, because the skip is invisible. The route is in
-- the plan: re-run the SQL stages of `GvsMapUnmappedVIDs` and `GvsMapDroppedDuplicateVIDs`
-- with the Step 7 predicate, evaluating `FT` per `(vid, input_location, input_ref,
-- input_alt)` and unioning the survivors, against a table that is now already populated --
-- so `GvsMapUnmappedVIDs` needs a DELETE step it does not currently have.
--
-- Step 6a runs AFTER STEP 6.8, against `vid_to_participant_mapping_2026_09_05`, and its own
-- DELETE-then-INSERT replaces those VIDs' rows wholesale. Nothing above needs to be redone
-- for it. Note that this makes `2026_09_05` a table whose contents change after its date,
-- which is a small wart in the naming convention; the alternative -- a third dated version
-- for 5,776 rows out of 1.6 billion -- costs another 19 TiB of storage to avoid it.
--
-- Then Step 8 re-exports Parquet -- a third delivery to the same collaborator, worth
-- flagging to the PO in advance rather than after.
