-- =====================================================================================
-- Scrub withdrawn/control person ids from an already-delivered participant mapping table
-- =====================================================================================
--
-- WHY: all three mapping-table workflows join `alt_allele` to `sample_info` on sample_id
-- with no filter on `withdrawn`. `alt_allele` is append-only -- GvsPopulateAltAllele.wdl
-- extends it only for sample ids above MAX(sample_id) already present -- so it retains
-- rows for samples withdrawn after those rows were written (e.g. Echo samples withdrawn
-- during the Foxtrot cycle). Those samples are excluded from the VDS, so they surface as
-- person ids in the mapping table that have no counterpart in the VDS.
--
-- This scrubs the delivered table in place. The pipeline fix (adding the filter to all
-- three WDLs) is separate and only affects future runs.
--
-- Substitute the placeholders throughout, e.g.:
--   sed -e 's/<PROJECT>/aou-genomics-curation-prod/g' \
--       -e 's/<DATASET>/foxtrot/g' \
--       -e 's/<MAPPING_TABLE>/vid_to_person_id/g' scrub_participant_mapping_table.sql
--
-- -------------------------------------------------------------------------------------
-- COST MODEL. On Foxtrot, one full pass over the mapping table's `person_ids` column was
-- measured at ~19 TiB; call that T. Everything expensive here is expensive for exactly
-- one reason -- touching `person_ids`. Bytes are billed per column per STATEMENT, so
-- repeated references within a single statement cost nothing extra; that is why the
-- validation below is one UNION ALL rather than several statements.
--
--   STEP 1  precondition        reads `vid` only, not `person_ids`   -- small
--   STEP 2  impact numbers      COUNT(*) is metadata-only            -- free
--   STEP 2b reported person ids `sample_info` only                   -- negligible
--   STEP 3  build               1 x T                                -- unavoidable
--   STEP 4a validate (required) 1 x T   (reads the scrubbed copy only)
--   STEP 4b validate (optional) 2 x T   (reads both copies; it is its own statement, so
--                                        its bytes are billed independently of 4a. If you
--                                        want both, paste them as ONE statement joined by
--                                        UNION ALL and the pair costs 2 x T, not 3 x T.)
--   STEP 5  swap                metadata rename                      -- free
--   APPENDIX full impact report 1 x T                                -- optional
--   STEP 4a-SMOKE  sampled 4a   ~0.01 x T                            -- optional gate
--
-- So the floor is 2T (build + required validation) -- about 38 TiB on Foxtrot. This runs
-- on ON-DEMAND billing, so those bytes are money: at roughly $6/TiB (US multi-region;
-- confirm the current rate) each pass over the table is on the order of $120. Skipping
-- the appendix and STEP 4b is therefore worth about $360, not merely some slot time.
--
-- Two consequences that are easy to get backwards:
--   * Bytes are billed as READ, not as shuffled. STEP 3 and STEP 3-ALT read the same
--     columns and so cost the same dollars; 3-ALT is a wall-clock play only.
--   * Cancelling a running query does not reliably refund it -- BigQuery may still bill
--     bytes processed up to the point of cancellation. Killing a long job to "save money"
--     generally does not.
-- -------------------------------------------------------------------------------------
--
-- Semantics are verified by verify_participant_mapping_scrub.sql (synthetic fixture,
-- all mutants of the logic caught). Every statement here has additionally been dry-run
-- against a real GVS `sample_info` to confirm it compiles -- see the note above STEP 4b
-- for why that second check is not redundant. Nothing here has been run against a real
-- mapping table; run steps 1 and 2b and read them before step 3.
-- =====================================================================================


-- -------------------------------------------------------------------------------------
-- STEP 1 -- PRECONDITION. Must return zero rows. Cheap: reads `vid` only.
--
-- Step 3 groups by vid, which would silently merge two rows sharing a vid into one. The
-- pipeline should already guarantee uniqueness (each of the three writers groups by vid,
-- and the dropped-duplicate step deletes before it inserts), but the table carries no
-- constraint enforcing it, so confirm rather than assume. If this returns rows, stop:
-- step 3 needs reworking to preserve row identity.
-- -------------------------------------------------------------------------------------
SELECT vid, COUNT(*) AS row_count
FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>`
GROUP BY vid
HAVING COUNT(*) > 1
LIMIT 100;


-- -------------------------------------------------------------------------------------
-- STEP 2 -- IMPACT NUMBERS, free version. COUNT(*) is answered from metadata.
--
-- Run the first query now and the second after step 3; the difference is the number of
-- VIDs that lost every carrier and were dropped. That is the headline number for the
-- collaborator. The per-person-id breakdown costs a full pass and is in the appendix.
-- -------------------------------------------------------------------------------------
SELECT COUNT(*) AS vid_rows_before FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>`;
-- after step 3:
-- SELECT COUNT(*) AS vid_rows_after FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed`;

-- Context from sample_info alone (negligible cost): how many samples the fix excludes.
SELECT
  COUNTIF(withdrawn IS NOT NULL)                                        AS withdrawn_samples,
  COUNTIF(is_control)                                                   AS control_samples,
  COUNTIF(withdrawn IS NOT NULL AND SAFE_CAST(sample_name AS INT64) IS NOT NULL)
                                                                        AS withdrawn_with_numeric_person_id
FROM `<PROJECT>.<DATASET>.sample_info`;


-- -------------------------------------------------------------------------------------
-- STEP 2b -- Targeted check on the person ids the collaborator reported. Negligible cost.
--
-- Confirms `withdrawn IS NULL` is actually the predicate that catches them. If any comes
-- back with withdrawn IS NULL and is_control = false, this scrub will NOT remove it and
-- the diagnosis is incomplete -- see the caveat at the bottom of this file.
-- -------------------------------------------------------------------------------------
SELECT sample_name, sample_id, is_control, withdrawn
FROM `<PROJECT>.<DATASET>.sample_info`
WHERE sample_name IN UNNEST(['<REPORTED_PERSON_ID_1>', '<REPORTED_PERSON_ID_2>']);


-- -------------------------------------------------------------------------------------
-- STEP 3 -- BUILD the scrubbed copy alongside the original. 1 x T.
-- Does not touch the original.
--
-- The source table is created by a plain CREATE TABLE AS SELECT with no partitioning or
-- clustering, so this CTAS reproduces its physical layout.
-- -------------------------------------------------------------------------------------
CREATE TABLE `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` AS
WITH good_person_ids AS (
  -- A person id is good iff at least one sample_info row bearing that sample_name is
  -- neither withdrawn nor a control. Stated this way, a person re-ingested under a new
  -- sample_id survives even though an older withdrawn row for the same name exists.
  -- is_control is REQUIRED in the sample_info schema, so `= false` has no NULL trap.
  SELECT DISTINCT SAFE_CAST(sample_name AS INT64) AS person_id
  FROM `<PROJECT>.<DATASET>.sample_info`
  WHERE withdrawn IS NULL
    AND is_control = false
    AND SAFE_CAST(sample_name AS INT64) IS NOT NULL
),
exploded AS (
  SELECT m.vid AS vid, p AS person_id
  FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS m, UNNEST(m.person_ids) AS p
)
-- INNER JOIN drops bad person ids; the GROUP BY then drops any vid left with none, which
-- is what the fixed pipeline produces for such a vid (its join yields no rows, so it gets
-- no row at all, rather than a row with an empty array).
--
-- No ORDER BY: person_ids order is not preserved, deliberately. The pipeline's own
-- ARRAY_AGG has no ORDER BY either, so the delivered order was already arbitrary and
-- nothing downstream can depend on it. Preserving it would mean carrying a WITH OFFSET
-- column through the shuffle -- an extra INT64 per array element, which at this table's
-- element count is a large amount of additional shuffle -- plus a per-group sort, all to
-- make the diff prettier. Not worth it.
SELECT e.vid AS vid, ARRAY_AGG(e.person_id) AS person_ids
FROM exploded AS e
JOIN good_person_ids AS g ON e.person_id = g.person_id
GROUP BY e.vid;


-- -------------------------------------------------------------------------------------
-- STEP 3-ALT -- row-local formulation. Same output, different execution shape.
--
-- Verified equivalent to STEP 3 above on the shared fixture: same row count, same vids,
-- same array contents including duplicate multiplicity (compared as multisets, since
-- neither formulation fixes an order). Use it only as a retry if STEP 3 is
-- too slow or dies with "Resources exceeded during query execution".
--
-- WHY IT MIGHT BE FASTER: STEP 3 explodes the table into one row per (vid, person_id),
-- which repeats the vid STRING for every array element and then shuffles all of that by
-- vid to regroup it. At Foxtrot scale the vid strings alone dominate the shuffle -- far
-- more bytes than the person ids they accompany. This form never explodes or regroups:
-- each row is filtered in place and the good-id array is broadcast, so there is no
-- shuffle and no GROUP BY.
--
-- WHY IT MIGHT BE SLOWER -- UNVERIFIED AT SCALE, and the reason this is not the default:
-- `p IN UNNEST(g.ids)` tests each person id against an array of every good person id
-- (~450k on Foxtrot). If BigQuery builds a hash set for that array it is O(1) per element
-- and this wins outright; if it degrades to a linear scan it is O(450k) per element and
-- this is far worse than the hash join in STEP 3. Which one happens cannot be established
-- from a small fixture -- correctness was verified, performance was not.
--
-- MEASURED ON FOXTROT, and the reason STEP 3 is no longer the recommended path: the
-- first STEP 3 run held all 2000 on-demand slots for 82 minutes and did not finish. Its
-- scan, join and repartition cascade completed cleanly in 18 minutes; the output stage
-- then wrote 11,120 of 20,008 units in nine minutes and collapsed from ~1,233 units/min
-- to ~0.4, flat for the following 44 minutes with worker concurrency slowly bleeding off.
-- MEASURED SHAPE OF THE TABLE (from a 1% TABLESAMPLE, scaled x100):
--   ~1.57 billion vid rows          ~2.58 trillion person-id entries
--   mean array 1647   p50 2   p99 10,586   max 540,539
--
-- THESE ARE PRE-SCRUB FIGURES and do not describe the delivered table. They were taken
-- before this scrub removed withdrawn and control participants, which is why `max` exceeds
-- the 535,662 active non-control participants in `sample_info` -- that 540,539-element
-- array included people the scrub then took out. Post-scrub no array can exceed 535,662.
-- Do not carry `max 540,539` forward as a property of the delivered table
-- `vid_to_participant_mapping_2026_08_28`.
--
-- Half of all vids carry two person ids or fewer, while the mean is ~800x the median, so
-- nearly all entries sit in a small minority of rows. That distribution is why STEP 3
-- fails, but NOT via partition skew -- an earlier reading of this. Hashing 1.57 billion
-- vids into 20,008 partitions puts the hot vids so evenly across them that the heaviest
-- partition is under a percent above average. The problem is simply volume: exploding the
-- arrays turns 1.57 billion rows into 2.58 TRILLION, each carrying a repeated vid string,
-- for a shuffle on the order of 70 TB. The clean cliff at 56% is what exceeding shuffle
-- capacity and spilling looks like, not what one hot worker looks like.
--
-- 3-ALT never explodes: it processes one row per vid, ~1.57 billion, filtering each array
-- in place. That is ~1,600x fewer rows through the engine and no shuffle at all. The
-- bytes READ are identical, so this is not a cost difference -- only a tractability one.
--
-- A TEMPTING OPTIMIZATION THAT IS NOT EQUIVALENT -- do not adopt it blind. The open risk
-- in 3-ALT is `p IN UNNEST(g.ids)` scanning ~450k good ids per element if BigQuery does
-- not hash it. Testing against the BAD set instead (withdrawn + control, thousands rather
-- than hundreds of thousands) would cut that worst case by ~100x:
--
--   WITH bad AS (
--     SELECT ARRAY_AGG(person_id IGNORE NULLS) AS ids FROM (
--       SELECT SAFE_CAST(sample_name AS INT64) AS person_id
--       FROM `<PROJECT>.<DATASET>.sample_info` GROUP BY sample_name
--       HAVING LOGICAL_AND(withdrawn IS NOT NULL OR is_control)))
--   ...  WHERE p NOT IN UNNEST(b.ids)
--
-- Verified on the fixture: it does NOT agree with the good-set form. They differ on a
-- person id that appears in the mapping table but in NO sample_info row -- good-set drops
-- it, bad-set keeps it. The good-set form is the correct one, because the fixed pipeline
-- joins alt_allele to sample_info and so can never emit such an id. Adopt the bad-set
-- form only if the sample race below shows it is both meaningfully faster AND produces
-- identical output on the sample, which is also the cheapest way to learn whether any
-- such person ids exist at all. STEP 4a would catch the divergence if one slipped through
-- (it tests membership in the good set), but only after you have paid for the build.
--
-- IF THE MEMBERSHIP TEST IS THE BOTTLENECK, hoist the bad ids out of the query entirely
-- and inline them as a CONSTANT array literal, which BigQuery can hash at plan time
-- rather than rescanning a correlated array per element. Verified equivalent on the
-- fixture. Generate the literal with:
--
--   SELECT CONCAT('[', STRING_AGG(CAST(person_id AS STRING), ','), ']') AS bad_ids_literal,
--          COUNT(*) AS n_bad
--   FROM (SELECT SAFE_CAST(sample_name AS INT64) AS person_id
--         FROM `<PROJECT>.<DATASET>.sample_info` GROUP BY sample_name
--         HAVING LOGICAL_AND(withdrawn IS NOT NULL OR is_control))
--   WHERE person_id IS NOT NULL;
--
-- then paste it in place of the CROSS JOIN:
--
--   SELECT vid, person_ids FROM (
--     SELECT m.vid AS vid,
--            ARRAY(SELECT p FROM UNNEST(m.person_ids) AS p
--                  WHERE p NOT IN UNNEST(<literal>)) AS person_ids
--     FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS m
--   )
--   WHERE ARRAY_LENGTH(person_ids) > 0;
--
-- Two details in that statement are load-bearing:
--   * `AS person_ids` is REQUIRED. An unaliased ARRAY(...) expression gets an autogenerated
--     column name, so a CTAS without it produces a table whose second column is not called
--     person_ids -- which then breaks STEP 4a, STEP 4b and the STEP 5 swap. Check the
--     resulting schema before going further.
--   * Dropping the emptied rows via a wrapping ARRAY_LENGTH(...) > 0 rather than a
--     WHERE EXISTS(...) on the same predicate. The EXISTS form is correct but evaluates the
--     entire membership test a second time, doubling the per-element work that is the
--     suspected bottleneck in the first place.
--
-- This only works for the BAD set. BigQuery caps query text at 1 MB, so a few thousand
-- withdrawn ids inline comfortably while the ~450k good ids (~3.5 MB) cannot -- which is
-- the one concrete advantage the bad-set formulation has, and the reason to settle its
-- ghost-id equivalence on the sample rather than dismiss it.
--
-- TO DECIDE CHEAPLY, build a sample table with TABLESAMPLE (which really does read fewer
-- blocks, unlike LIMIT) and race the two formulations against it:
--   CREATE TABLE `<PROJECT>.<DATASET>.<MAPPING_TABLE>_sample` AS
--   SELECT * FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS m TABLESAMPLE SYSTEM (1 PERCENT);
--   (the alias, if any, must come BEFORE TABLESAMPLE -- the other order is a syntax error)
-- then compare total_slot_ms for each, from the diagnostic in the notes below.
-- -------------------------------------------------------------------------------------
-- CREATE TABLE `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` AS
-- WITH good AS (
--   SELECT ARRAY_AGG(DISTINCT SAFE_CAST(sample_name AS INT64) IGNORE NULLS) AS ids
--   FROM `<PROJECT>.<DATASET>.sample_info`
--   WHERE withdrawn IS NULL AND is_control = false
-- )
-- SELECT m.vid AS vid,
--        ARRAY(SELECT p FROM UNNEST(m.person_ids) AS p WHERE p IN UNNEST(g.ids)) AS person_ids
-- FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS m CROSS JOIN good AS g
-- -- reproduces the GROUP BY's drop semantics: a vid with no surviving carrier gets no row
-- WHERE EXISTS(SELECT 1 FROM UNNEST(m.person_ids) AS p WHERE p IN UNNEST(g.ids));


-- -------------------------------------------------------------------------------------
-- STEP 4a -- VALIDATE, required. 1 x T (reads the scrubbed copy only).
-- Both rows must read 'pass'.
--
-- These two catch the failure mode that is actually plausible here: the predicate being
-- wrong. They do not reference the original table, which is what keeps them at 1 x T.
-- -------------------------------------------------------------------------------------
-- Phrased as POSITIVE membership in the small bad-id literal, not absence from the large
-- good-id set. The good-set form -- `p NOT IN UNNEST((SELECT ids FROM good))` -- is the
-- same correlated-array pattern that made STEP 3-ALT take 12+ minutes on 1% of the table
-- before the literal replaced it; run at full scale it would not return. The literal is a
-- constant BigQuery can hash at plan time. It is also NULL-safe by construction, since
-- positive IN has none of the NOT IN NULL trap, and the generator excludes NULL ids.
-- This substitution is only valid because the ghost-id check returned 0: with no person
-- ids outside sample_info, "is a bad id" and "is not a good id" are the same predicate.
SELECT 'no withdrawn/control person ids remain' AS check_name,
       IF(NOT EXISTS(
            SELECT 1
            FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` AS m, UNNEST(m.person_ids) AS p
            WHERE p IN UNNEST(<BAD_IDS_LITERAL>)
          ), 'pass', 'FAIL') AS result
UNION ALL
SELECT 'no rows left with an empty person_ids array',
       IF(NOT EXISTS(
            SELECT 1 FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed`
            WHERE ARRAY_LENGTH(person_ids) = 0
          ), 'pass', 'FAIL')
UNION ALL
-- Free (COUNT(*) is metadata). Report the delta rather than asserting equality.
--
-- The near-invariant: a vid reaches the mapping table from the VAT, the VAT is built from
-- the VDS, and the VDS excludes withdrawn samples -- so the VDS carrier that put the vid
-- in the VAT is normally still there after scrubbing and the row survives. That holds for
-- rows written by GvsCreateParticipantMappingTable, whose carriers come from an exact
-- (location, ref, alt) join to alt_allele.
--
-- It does NOT hold for rows written by GvsMapUnmappedVIDs or GvsMapDroppedDuplicateVIDs.
-- Those steps exist precisely because the vid did not match alt_allele in that exact form;
-- they re-normalize and attach whatever carriers the realigned representation finds, which
-- need not include the VDS carrier that justified the vid's presence in the VAT. Such a row
-- can therefore consist entirely of withdrawn samples and legitimately disappear here.
--
-- So: 0 is expected, a handful is plausible and worth identifying individually (see the
-- note below on locating them cheaply), and a large number means the scrub is removing
-- more than it should -- investigate before swapping.
SELECT 'row count delta (expect 0 or a very small positive number)',
       CAST((SELECT COUNT(*) FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>`)
          - (SELECT COUNT(*) FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed`) AS STRING);

-- To identify dropped vids WITHOUT paying for a pass over person_ids: anti-join on vid
-- alone. Bytes are billed per column, and the vid columns are a rounding error next to
-- person_ids, so this costs a fraction of a percent of a full pass.
--   SELECT o.vid
--   FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS o
--   LEFT JOIN `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` AS s ON o.vid = s.vid
--   WHERE s.vid IS NULL;
-- Do NOT add o.person_ids to that SELECT -- referencing the column bills the whole 19 TiB
-- regardless of how few rows come back.


-- -------------------------------------------------------------------------------------
-- STEP 4a-SMOKE -- optional cheap gate to run BEFORE the full 4a. ~1% of a pass.
--
-- TABLESAMPLE really does read fewer blocks (unlike LIMIT), so this costs on the order of
-- a dollar. It is a smoke test, NOT a proof: a systematically wrong predicate would leave
-- violations in essentially every block and be caught here, but a handful of stragglers
-- would slip through. Its value is sequencing -- if this fails you have learned the scrub
-- is wrong without paying for the full pass.
-- Treat a 'pass' here as permission to run the real 4a, never as a substitute for it.
-- -------------------------------------------------------------------------------------
SELECT
  COUNTIF(p IN UNNEST(<BAD_IDS_LITERAL>))          AS violations_in_sample,
  COUNT(*)                                         AS person_id_entries_sampled,
  COUNT(DISTINCT m.vid)                            AS vids_sampled
FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` AS m TABLESAMPLE SYSTEM (1 PERCENT),
     UNNEST(m.person_ids) AS p;
-- Read all three numbers, not just the first. violations_in_sample must be 0, but a 0
-- alongside a tiny person_id_entries_sampled means the sample was too thin to conclude
-- anything -- the gate only carries weight when it actually looked at a lot of entries.
--
-- Two TABLESAMPLE restrictions worth knowing, both of which reject the query outright:
-- the sampled table may not be referenced elsewhere in the same query, and sampling may
-- not appear inside an IN/EXISTS subquery. Both are satisfied above because `good` reads
-- sample_info while the sample is taken of the scrubbed mapping table -- but if you
-- adapt this, note that pointing both at the same table will not run.


-- -------------------------------------------------------------------------------------
-- STEP 4b -- VALIDATE, optional cross-check against the original. Adds 1 x T.
--
-- Guards against the CTAS inventing or substituting data rather than only removing it.
--
-- DO NOT RUN THIS AT FOXTROT SCALE AS WRITTEN. It still uses the large correlated good-id
-- array (`p IN UNNEST(g.ids)`), the pattern that proved pathological in STEP 3-ALT before
-- the constant literal replaced it. Unlike 4a it cannot simply be switched to the bad-id
-- literal, because it reconstructs the expected array rather than testing membership.
-- 4a plus the free row-count invariant is the right place to stop on a table this size.
--
-- Note the shape of the content check. The natural way to write it -- a correlated
-- ARRAY subquery whose WHERE contains `IN (SELECT ... FROM sample_info)` -- parses fine
-- but is rejected at execution against real tables with "Correlated subqueries that
-- reference other tables are not supported unless they can be de-correlated". A fixture
-- built from CTE literals does NOT reproduce that error, so this form is what has to be
-- dry-run against a real table. Hence `IN UNNEST(g.ids)` over a CROSS JOINed array.
-- -------------------------------------------------------------------------------------
WITH good AS (
  SELECT ARRAY_AGG(DISTINCT SAFE_CAST(sample_name AS INT64) IGNORE NULLS) AS ids
  FROM `<PROJECT>.<DATASET>.sample_info`
  WHERE withdrawn IS NULL AND is_control = false
)
SELECT 'no vid gained a row' AS check_name,
       IF(NOT EXISTS(
            SELECT s.vid FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` AS s
            LEFT JOIN `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS o ON s.vid = o.vid
            WHERE o.vid IS NULL
          ), 'pass', 'FAIL') AS result
UNION ALL
-- Compared as multisets -- both sides sorted -- because STEP 3 no longer fixes an order.
-- Sorting rather than deduplicating keeps duplicate person ids significant, so a scrub
-- that dropped or added a repeat of an id still fails this.
SELECT 'every surviving array is exactly the original minus its bad person ids',
       IF(NOT EXISTS(
            SELECT 1
            FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` AS s
            JOIN `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS o ON s.vid = o.vid
            CROSS JOIN good AS g
            WHERE TO_JSON_STRING(ARRAY(SELECT p FROM UNNEST(s.person_ids) AS p ORDER BY p))
               != TO_JSON_STRING(ARRAY(SELECT p FROM UNNEST(o.person_ids) AS p
                                       WHERE p IN UNNEST(g.ids) ORDER BY p))
          ), 'pass', 'FAIL');


-- -------------------------------------------------------------------------------------
-- STEP 4c -- COVERAGE. Every active non-control sample should still have mapping entries.
-- ~0.01 x T: run it against the SCRUBBED 1% SAMPLE, not the full table.
--
-- Zero missing on the sample is conclusive for the whole table, and this asymmetry is the
-- point. The scrub is row-local, so scrubbing a subset equals a subset of scrubbing, which
-- makes the sampled scrubbed table a true subset of the full scrubbed table. Any person id
-- present in the sample is therefore present in the full table. Only a NON-zero result is
-- ambiguous -- it cannot distinguish "absent everywhere" from "absent from this 1%" -- and
-- only that case needs the full 1 x T pass.
--
-- It is also sensitive despite being 1%: every sample carries millions of variants and so
-- appears in millions of rows, of which a random 1% still leaves tens of thousands. A
-- sample missing here is genuinely unrepresented, which is the failure mode worth catching.
--
-- Like the ghost-id check, the UNNEST explodes but DISTINCT over a domain bounded at the
-- sample count pre-aggregates locally, so almost nothing reaches the shuffle -- nothing
-- like the 1.6-billion-group GROUP BY that made STEP 3 fail.
--
-- Run it BEFORE dropping the sample tables. Foxtrot result: 0.
-- -------------------------------------------------------------------------------------
WITH present AS (
  SELECT DISTINCT p AS person_id
  FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>_sample_scrubbed` AS m, UNNEST(m.person_ids) AS p
)
SELECT COUNT(*) AS active_samples_with_no_mapping_entry
FROM `<PROJECT>.<DATASET>.sample_info` si
LEFT JOIN present ON present.person_id = SAFE_CAST(si.sample_name AS INT64)
WHERE si.withdrawn IS NULL
  AND si.is_control = false
  AND present.person_id IS NULL;


-- -------------------------------------------------------------------------------------
-- STEP 5 -- SWAP, retaining the pre-scrub table as a backup. Free (metadata rename).
--
-- Run only after step 4a is all 'pass'. Drop the backup once the collaborator confirms.
-- -------------------------------------------------------------------------------------
ALTER TABLE `<PROJECT>.<DATASET>.<MAPPING_TABLE>`          RENAME TO `<MAPPING_TABLE>_prescrub_backup`;
ALTER TABLE `<PROJECT>.<DATASET>.<MAPPING_TABLE>_scrubbed` RENAME TO `<MAPPING_TABLE>`;


-- =====================================================================================
-- APPENDIX -- full per-person-id impact report. 1 x T. Purely diagnostic; nothing
-- depends on it. Run only if the collaborator wants the detailed breakdown.
-- =====================================================================================
-- WITH good AS (
--   SELECT ARRAY_AGG(DISTINCT SAFE_CAST(sample_name AS INT64) IGNORE NULLS) AS ids
--   FROM `<PROJECT>.<DATASET>.sample_info`
--   WHERE withdrawn IS NULL AND is_control = false
-- ),
-- classified AS (
--   SELECT m.vid AS vid, p AS person_id, p IN UNNEST((SELECT ids FROM good)) AS is_good
--   FROM `<PROJECT>.<DATASET>.<MAPPING_TABLE>` AS m, UNNEST(m.person_ids) AS p
-- )
-- SELECT
--   COUNT(DISTINCT IF(NOT is_good, person_id, NULL)) AS distinct_bad_person_ids,
--   COUNTIF(NOT is_good)                             AS bad_person_id_entries,
--   COUNT(DISTINCT IF(NOT is_good, vid, NULL))       AS vids_touched
-- FROM classified;


-- =====================================================================================
-- BLOCKING PRE-SWAP CHECK -- reconcile sample_info.withdrawn against the VDS
--
-- This scrub removes every person id whose sample_info rows are all withdrawn or control.
-- That is correct ONLY IF `withdrawn IS NOT NULL` implies the sample is actually absent
-- from the delivered VDS. Nothing in the code enforces that: GvsWithdrawSamples stamps a
-- timestamp in BigQuery, while GvsMergeAndRescoreVDSes removes samples from the VDS by
-- research id off a separately maintained TSV. The two can disagree in both directions.
--
-- The direction that matters for correctness is a sample marked withdrawn in BigQuery but
-- still present in the VDS. For such a sample the mapping rows are RIGHT, and this scrub
-- deletes them -- turning a fix into data loss across every entry it touches, not just the
-- rows that vanish outright.
--
-- What surfaced this on Foxtrot: two vids lost their only person id and disappeared
-- (1-143186828-T-TG and 17-72146919-GAGAAG...-G, each a single withdrawn carrier). Both
-- are in the VAT and both match alt_allele exactly at the expected location, so this is
-- NOT the normalization/pseudo-vid phenomenon documented under
-- scripts/variantstore/scripts/variant_annotation_table/left_alignment_fixups/.
-- And it should not happen: hail_gvs_util.py's filter_samples_and_remove_monomorphic_rows
-- removes samples with remove_dead_alleles=True and then drops rows with no remaining
-- non-ref calls, so a site whose only carrier was removed cannot survive in the VDS -- and
-- its vid should not be in a VAT built from that VDS.
--
-- So before swapping, confirm the withdrawn set is a subset of what the VDS excludes:
--   1. Cheap first pass -- are the person ids of any vids dropped by this scrub on the VDS
--      sample-removal TSV (foxtrot_samples_to_withdraw_from_echo_*.tsv per
--      AOU_DELIVERABLES.md)? Locate those ids with the vid anti-join above.
--   2. The real check -- every sample_name with withdrawn IS NOT NULL should be absent from
--      the delivered VDS's column list.
-- If any withdrawn sample is still in the VDS, do NOT swap: the scrub is over-removing.
--
-- FOXTROT, checked 2026-08. The ledger reconciled exactly, with every sample_info row
-- accounted for and nothing left over:
--     535,662  active non-control  == VDS column count, and every VDS column joined to a
--                                     sample_info row (so the join key is sample_name)
--       4,875  numeric all-withdrawn names == the size of the bad-ids literal
--           8  control samples with non-numeric names
--     -------
--     540,545  sample_info rows
-- The 8 cannot appear in person_ids at all, since those are INT64 and the names do not
-- cast -- which is the same incidental protection that made controls a non-issue before
-- `is_control = false` was added to the WDLs. No duplicate/re-ingested sample_names exist
-- in this dataset, so the "at least one active row" phrasing of good_person_ids never
-- diverged from a simpler "not withdrawn" test here; keep the phrasing anyway, since a
-- future callset that re-ingests a withdrawn person would need it.
--
-- STEP R1 on Foxtrot returned 0 / 0 / 0: no withdrawn sample present in the VDS, no active
-- sample missing from it, no control in it. The VDS column set is exactly the active
-- non-control set, so every person id this scrub removes belongs to a sample genuinely
-- absent from the delivered VDS. Gate passed; the scrub is correct, not over-removing.
-- =====================================================================================


-- =====================================================================================
-- RELATED FINDING -- person_ids can be WRONG, not merely stale (Foxtrot, 2026-08)
--
-- Traced while explaining the two vids this scrub dropped entirely. It looked at first
-- like an independent, broader defect; it is not. It is a consequence of the same stale
-- alt_allele rows, and the withdrawn/control filter added to the three WDLs fixes it --
-- see SCOPE below. This scrub does not fix it, because the scrub only removes the bad
-- person id; it cannot add the correct one.
--
-- 1-143186828-T-TG had person_ids = [<PERSON_C>], a withdrawn sample absent from the VDS.
-- The variant itself is real and IS in the VDS -- as chr1:143186829 G>GG, the same
-- single-G insertion shifted one base right, carried by ONE active sample. Both forms
-- exist in alt_allele because different input gVCFs normalize indels differently.
--
-- GvsCreateParticipantMappingTable joins a VAT vid to alt_allele on the vid's OWN
-- coordinates. VAT vids are left-aligned; GVS representations often are not. So the join
-- matched the withdrawn sample (whose gVCF happened to use the left-aligned form) and
-- never saw the active carrier one base to the right. The delivered mapping was not
-- incomplete -- it named the wrong person.
--
-- GvsMapUnmappedVIDs cannot repair this: it only sees vids that matched NOTHING, and a
-- partial match is still a match. The vid was never classified unmapped.
--
-- SCOPE -- narrower than it first appears, and NOT a second independent defect. Writing
-- R_LA for the left-aligned form (always the VAT vid, since bcftools norm left-aligns)
-- and R_x for a shifted one:
--
-- Let A be the set of representations used by the ACTIVE (in-VDS) carriers. The base join
-- matches every alt_allele row at R_LA, from any sample. The split is on whether any
-- active carrier uses the left-aligned form:
--
--   1. A = {R_LA}                     -> join finds all active carriers; correct
--   2. R_LA in A, |A| > 1             -> join finds only the R_LA subset, but the VDS holds
--                                        a record for each form, so normalization collides
--                                        them onto one vid and GvsMapDroppedDuplicateVIDs
--                                        repairs it
--   3. R_LA not in A, and no non-VDS  -> join finds nothing; the vid is unmapped and
--      sample sits at R_LA               GvsMapUnmappedVIDs repairs it
--   4. R_LA not in A, but a WITHDRAWN -> join matches only that non-VDS sample. A row is
--      or CONTROL sample sits at R_LA    written, so the vid never looks unmapped and
--                                        nothing repairs it.
--
-- Withdrawn samples appear only in case 4. Elsewhere their rows either coincide with the
-- active carriers' (ordinary contamination, which this scrub removes) or sit at a
-- representation the join never touches, where they are harmless.
--
-- Cases 3 and 4 are identical from the active carriers' point of view and differ only in
-- whether one stale row exists -- which is why removing stale rows turns 4 into 3.
--
-- Only case 4 is unhandled, and the stale alt_allele row is its CAUSE rather than an
-- incidental detail: it manufactures a spurious match that makes an otherwise-unmapped vid
-- look mapped, suppressing the detection that would have repaired it. Without that sample
-- case 4 is case 2.
--
-- CONSEQUENCE: the withdrawn/control filter added to the three WDLs fixes THIS TRIGGER of
-- case 4. The join yields no rows, so no mapping row is written, so
-- GvsMapUnmappedVIDs.wdl:77 (vid NOT IN (SELECT vid FROM mapping)) selects the vid and the
-- synonym search maps the real carriers. The is_control half is load-bearing, not merely
-- defensive: a control-only match produces a row with an EMPTY person_ids array
-- (ARRAY_AGG ... IGNORE NULLS over all-NULL names), which still counts as mapped and
-- suppresses detection identically.
--
-- CASE 4 HAS A SECOND TRIGGER THAT THE WDL FILTER DOES NOT FIX. The spurious R_LA row need
-- not belong to a withdrawn or control sample. It is enough that R_LA was excluded from the
-- VAT while R_x was kept, and VETS can do exactly that: filter_set_info is keyed by
-- location/ref/alt, so synonymous representations are scored independently and R_LA can
-- fail calibration sensitivity while R_x passes. GvsCreateVATFromVDS.wdl normalizes
-- `filtered_sites_only.bcf`, so filtering precedes normalization -- only R_x is normalized,
-- one record results, no duplicate is detected, and GvsMapDroppedDuplicateVIDs does not
-- fire. alt_allele retains R_LA regardless of filtering and its carrier is an ordinary
-- active sample, so the withdrawn/control predicate passes them through and the vid is
-- mapped to a carrier whose variant failed filtering while the passing carrier is missed.
--
-- THE GENERAL STATEMENT: the base join must be restricted to exactly the population the VAT
-- represents, on two axes -- sample eligibility (active, non-control; added by the commit)
-- and variant quality (passing the callset's calibration-sensitivity threshold; nothing
-- does this today). A mismatch on either axis manufactures a spurious match and suppresses
-- unmapped detection identically.
--
-- Unlike the withdrawn trigger this one is DETECTABLE: join mapping rows back to
-- filter_set_info and look for source alleles that fail. Expensive at full scale, but a
-- bounded sample establishes whether it occurs in practice.
--
-- STILL OPEN: case 3 depends entirely on GvsMapDroppedDuplicateVIDs, whose DELETE-then-
-- INSERT (lines 133, 146) re-inserts carriers joined on dup.input_* only. Whether that
-- preserves carriers found at the SURVIVING representation depends on whether
-- dropped_duplicate_mappings.tsv contains every synonym in the cluster or only the dropped
-- ones -- unverified. If only the dropped ones, case 3 is broken and the WDL filter does
-- not address it.
-- =====================================================================================


-- =====================================================================================
-- CAVEAT
--
-- This scrub trusts `sample_info.withdrawn` to be the complete record of what was kept
-- out of the VDS. GvsMergeAndRescoreVDSes / merge_and_rescore_vdses.py remove samples
-- from the VDS by research id, and nothing in the CODE guarantees the two agree in either
-- direction -- a sample marked withdrawn but left in the VDS, or removed from the VDS but
-- never marked withdrawn.
--
-- RESOLVED FOR FOXTROT by the reconciliation above: STEP R1 returned 0 / 0 / 0, so both
-- directions are clean and the assumption holds for this callset. It is an assumption
-- about data, not an invariant enforced by code, so re-run the reconciliation on the next
-- callset rather than inheriting this result.
-- =====================================================================================
