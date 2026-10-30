-- =====================================================================================
-- Step 5 -- the one full pass over `alt_allele`, sharded per chromosome
-- =====================================================================================
--
-- Builds `exclusions`, a table of `(vid, person_id)` pairs naming every entry the
-- delivered participant mapping table holds but should not: a carrier whose genotype
-- was a GQ 0 no-call, or whose genotype fails the filter model's `FT`.
--
-- The predicate is the one validated in Step 4 against the delivered VDS -- 13,785,178
-- comparable genotypes, 100.0000% agreement, and 238,792 adjudications all resolving in
-- favor of it. Nothing here is new science; this is the same arithmetic applied to the
-- whole genome instead of to 2 Mb of chr21.
--
-- SHAPE OF THE JOB
--
--   Section A   creates the empty destination once.
--   Section B   is the per-chromosome shard. Run it 24 times, editing `chrom`.
--   Section C   is a bash driver for those 24 runs.
--   Section D   reconciles the result and emits per-VID relative inflation.
--
-- About 80 TiB scanned in total, roughly $501 at on-demand rates, dominated by
-- `alt_allele`. Sharding does not change that -- the same bytes are read either way --
-- but it buys restartability. The preflight measured `alt_allele` at 2.84 trillion rows
-- and chr21 at 1.769% of them, so scaling Step 3's 22 minutes puts the largest shard
-- (chr2) near 1.6 hours and the whole run near 21 sequential; a single whole-genome job
-- would be that same commitment with no partial progress on failure. Each shard here is
-- `DELETE` then `INSERT` over one partition; BigQuery's `INSERT ... SELECT` is atomic, so
-- a failed shard leaves that chromosome empty rather than half-written, and re-running it
-- is safe with no manual cleanup.
--
-- Expect roughly 2.5 billion rows. Scale by chr21's share of VAT-restricted MAPPING
-- ENTRIES (~1.40%: 35.85 billion of ~2.56 trillion), not by its 1.769% share of
-- `alt_allele` ROWS -- chr21 has a lower VAT-restricted fraction than the genome average,
-- so the two shares differ by a factor of 1.26 and using the wrong one understates this
-- table by a third. Row share governs bytes and wall clock; entry share governs row counts.
-- chr21 is the shard to check that against, and the reason Section C runs it first.
--
-- WITHDRAWN AND CONTROL SAMPLES ARE FILTERED AT THE START
--
-- Step 3 deliberately did not filter them, to stay comparable with the original
-- measurement, and QUERY E showed the consequence: 540,545 distinct samples carry chr21
-- rows, more than the delivered callset. `populate_alt_allele_table.py:36,42` excludes
-- withdrawn samples only as of the moment each vet table was processed, so anyone
-- withdrawn later is still present, and controls are never excluded at all.
--
-- Step 5 writes real `(vid, person_id)` rows against a table that has already been
-- scrubbed of those people, so rows naming them are wasted work. Unlike the VAT
-- restriction of defect 3, this filter has no apply-it-after-`FT` trap: it is a predicate
-- on SAMPLES, and no sample's `FT` depends on another's.
--
-- The `si.withdrawn IS NULL AND si.is_control = false` join condition below is copied
-- from `GvsCreateParticipantMappingTable.wdl:66-67` rather than reinvented, so the
-- population here is by construction the population the table is supposed to contain.
-- `SAFE_CAST(si.sample_name AS INT64)` is likewise that workflow's own derivation of
-- `person_id`; a name that does not cast produced no entry there, so it must produce no
-- exclusion here.
--
-- THE VAT RESTRICTION COMES AFTER `FT`, NOT BEFORE
--
-- This is defect 3, the one that moved the chr21 rate by 2.77x, and the single thing
-- most worth not getting wrong again. `FT` is a property of a GENOTYPE and is computed
-- over every called non-ref allele of it. VAT membership is a property of an ALLELE and
-- is settled afterwards -- an allele whose every carrier fails reaches `AC = 0` and never
-- becomes a VID at all. So `any_yes` and `all_ok` below are aggregated over all of a
-- genotype's `alt_allele` rows with no reference to `vid`, and `WHERE s.vid IS NOT NULL`
-- appears only at the point where rows are emitted.
--
-- WHY AN ANALYTIC WINDOW RATHER THAN GROUP BY PLUS JOIN
--
-- Step 3 grouped to one row per genotype and joined that back to the rows -- two
-- shuffles. `LOGICAL_OR(...) OVER (PARTITION BY sample_id, location)` gets the same
-- answer in one, which at chr1 scale is the difference worth having. If BigQuery rejects
-- `LOGICAL_OR`/`LOGICAL_AND` as analytic functions in your region, the mechanical
-- fallback is `MAX(CAST(is_yes AS INT64)) OVER w = 1` and
-- `MIN(CAST(is_ok AS INT64)) OVER w = 1`, which are exactly equivalent over non-null
-- inputs -- and both `is_yes` and `is_ok` are non-null by construction, being `IFNULL`ed.
-- Failing that, lift the `genotypes` CTE from `step3_rescore_chr21.sql` verbatim.
--
-- Only het-var genotypes actually need the partitioning -- single-alt `FT` reduces to the
-- allele's own status and is row-local -- and restricting the window to multi-row
-- genotypes would cut the partitioned set by roughly 31x. That optimization is NOT taken
-- here: Step 3 ran the readable form at chr21 scale without trouble, and the readable
-- form is the one whose agreement with the VDS was measured. Reach for it only if a
-- chromosome stalls.
--
-- Note also that the arity test is "more than one `alt_allele` row at this
-- `(sample_id, location)`", never a `call_GT` literal list. The two agree on Foxtrot as
-- far as the ingest code goes, but the row-count test does not depend on that list
-- staying in sync, and per defect 4 there are 47,632 chr21 entries where they disagree
-- for a reason not yet established.
--
-- WHAT THIS DOES NOT COVER
--
-- The 5,776 VIDs rewritten by `GvsMapUnmappedVIDs` and `GvsMapDroppedDuplicateVIDs` carry
-- persons found at a DIFFERENT allele representation, so exclusions computed at the VID's
-- own coordinates do not reach them. That is Step 6a, and it is deliberately not folded
-- in here. If you would rather fold it in, the route is described in the plan: replace
-- those VIDs' rows in `allele_status` with one row per `(vid, input_location, input_ref,
-- input_alt)` before running this, using the `input_*` coordinates ONLY.
--
-- -------------------------------------------------------------------------------------


-- =====================================================================================
-- PREFLIGHT -- confirm the chromosome range covers everything, and size the shards.
-- Metadata only; free.
-- =====================================================================================
--
-- The shard loop runs 1..24 (23 = X, 24 = Y). Anything outside that -- chrM, or a
-- location encoded from a contig this assumes away -- would be silently skipped, and a
-- silently skipped partition is exactly the failure this sharding is meant to make
-- impossible.
--
-- `alt_allele` is partitioned `RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000,
-- 1000000000000))` (GvsPopulateAltAllele.wdl:253) -- one partition per chromosome, the
-- same scheme Section A gives `exclusions`. So the partition IDs *are* the answer, and
-- `INFORMATION_SCHEMA.PARTITIONS` reads them from metadata rather than from 80 TiB of
-- data. Out-of-range locations land in `__UNPARTITIONED__` and nulls in `__NULL__`;
-- neither casts to INT64, so `chrom IS NULL` catches both without depending on which one
-- BigQuery chooses.
--
-- Two conditions to check. Every row's `chrom` is between 1 and 24, and there are exactly
-- 24 of them -- a missing chromosome is as much a problem as an extra one, since the loop
-- would run a shard that reads nothing and report success.
--
-- The row and byte counts are the reason to prefer this over the data-scanning version
-- even setting cost aside: they size every shard in advance. chr21 took 22 minutes in
-- Step 3, so `total_rows / chr21_rows` predicts each chromosome's wall clock directly,
-- and `SUM(total_logical_bytes)` reconciles against the 80 TiB / ~$501 estimate before a
-- cent is spent rather than after.

SELECT
  SAFE_CAST(partition_id AS INT64) / 1000000000000     AS chrom,
  partition_id                                         AS raw_partition_id,
  total_rows,
  ROUND(total_logical_bytes / POW(1024, 4), 2)         AS logical_tib,
  ROUND(100 * total_rows
        / SUM(total_rows) OVER (), 3)                  AS pct_of_rows,
  last_modified_time
FROM `foxtrot.INFORMATION_SCHEMA.PARTITIONS`
WHERE table_name = "alt_allele"
ORDER BY chrom;

-- RUN 2026-09-04. Clean: 24 partitions, chrom 1 through 24, nothing unpartitioned or
-- null, so the shard loop covers the table exactly. Three things came out of it that
-- change the numbers used elsewhere in this file.
--
--   Total: 2,840,988,469,953 rows across 428.1 TiB logical.
--
-- 1. chr21 is 1.769% of rows, not the 1.3% the plan assumed -- the genome is 56.5 chr21s,
--    not 77. This corrects the per-shard byte and runtime projections and NOTHING ELSE. It
--    does not resize `exclusions`: that is an entry count, and chr21's share of
--    VAT-restricted mapping entries is ~1.40%, not 1.769%, because a lower fraction of
--    chr21's `alt_allele` rows are VAT alleles. Scaling D1's measured 35,052,223 by 1.40%
--    gives ~2.5 billion, the plan's original figure. An earlier revision DOWN to ~1.84
--    billion applied the row share to an entry count and was wrong.
--
-- 2. 428 TiB is the whole table, every column, and it does NOT contradict the estimate for
--    this job, which is column-pruned to the five columns Section B reads. SETTLED by
--    `--dry_run` on the chr21 shard, 2026-09-04: 1,558,142,894,853 bytes = 1.417 TiB. At
--    chr21's 1.769% row share that scales to ~80.1 TiB genome-wide, ~$501 at $6.25/TiB --
--    18.7% of the table, and under the ~$575 the plan budgeted. The pruning factor is the
--    check worth reading: near 100% would have meant a predicate defeating partition
--    pruning, and 18.7% is what five narrow-ish columns of a wide table should cost.
--
-- 3. Every partition carries the same `last_modified_time` to the millisecond,
--    2025-07-31 10:35:45.192 UTC, and this is expected rather than notable. `alt_allele`
--    is partitioned by LOCATION but populated by SAMPLE -- one `INSERT` per vet table
--    (`populate_alt_allele_table.py:30`), each carrying a batch of samples' variants
--    genome-wide -- so every insert writes into all 24 partitions and they all end up
--    stamped with the last one. Partition freshness therefore says nothing about how the
--    table was built or which vet tables landed when. What it does give is a single
--    as-of date for the whole table: two days after the `foxtrot_v4_2025_07_29` filter
--    set and about ten months before the `2026_06_01` mapping tables, so anyone withdrawn
--    in that window is still present. That is the withdrawal lag, dated -- and the reason
--    the `si.withdrawn IS NULL` join below reads `sample_info` at query time rather than
--    trusting what `alt_allele` was filtered on at write time.
--
-- Per-chromosome wall clock, scaling chr21's 22 minutes in Step 3 by row count: chr2 and
-- chr1 near 1.6 hours each, chr22 and chrY under half an hour, about 20.7 hours for all
-- 24 run sequentially. Treat that as an upper bound -- Section B does one shuffle where
-- Step 3 did two.
--
-- Note what this does NOT check, since the version it replaces did not either: that a
-- location *within* a chromosome's range is sane. Partition metadata cannot see inside a
-- partition. That question belonged to F2 in `step4_followups.sql`, which has since run
-- and closed it: chrX and chrY carry `0/1 0|1 1 1/1 1/2 1|0 1|1` and nothing else, so
-- there is no haploid `2` and no member of the `0/2` family anywhere on the sex
-- chromosomes. The `alt[0]` misattribution at `alt_allele_positions.sql:2-3` is therefore
-- unreachable, and shards 23 and 24 need no special handling.


-- =====================================================================================
-- SECTION A -- create the destination. Run once.
-- =====================================================================================
--
-- Partitioned by chromosome so a shard can replace its own slice, and clustered by `vid`
-- because Step 6 groups by it and the delivered mapping table is clustered on it too
-- (`GvsCreateParticipantMappingTable.wdl:89`).
--
-- `location` is carried for the partitioning and for per-chromosome reconciliation; it is
-- redundant with `vid`, which encodes it. `is_gq0` and `is_ft_fail` are carried so the
-- decomposition can be audited later without re-reading `alt_allele` -- at 80 TiB, being
-- able to answer a follow-up question from the output rather than from the input is worth
-- two boolean columns.
--
-- `CREATE TABLE IF NOT EXISTS`, not `CREATE OR REPLACE`: re-running this section by
-- accident partway through the shard loop must not silently discard finished shards.

CREATE TABLE IF NOT EXISTS `foxtrot.exclusions`
(
  vid         STRING  NOT NULL,
  person_id   INT64   NOT NULL,
  location    INT64   NOT NULL,
  is_gq0      BOOL    NOT NULL,
  is_ft_fail  BOOL    NOT NULL
)
PARTITION BY RANGE_BUCKET(location, GENERATE_ARRAY(0, 25000000000000, 1000000000000))
CLUSTER BY vid
OPTIONS (description = "Step 5: (vid, person_id) entries the Foxtrot participant mapping table over-reports -- GQ 0 no-calls and FT-failing genotypes. Predicate validated against the delivered VDS in Step 4.");


-- =====================================================================================
-- SECTION B -- one chromosome. Run 24 times, editing `chrom`.
-- =====================================================================================
--
-- Run these two statements together as one script so the `DECLARE`s apply to both. The
-- `DELETE` makes the shard idempotent; on a first run it removes nothing and costs
-- nothing, since it prunes to a single empty partition.

DECLARE chrom        INT64 DEFAULT 1;     -- 1..22, 23 = X, 24 = Y
DECLARE window_start INT64 DEFAULT chrom       * 1000000000000;
DECLARE window_end   INT64 DEFAULT (chrom + 1) * 1000000000000 - 1;

DELETE FROM `foxtrot.exclusions`
WHERE location BETWEEN window_start AND window_end;

INSERT INTO `foxtrot.exclusions` (vid, person_id, location, is_gq0, is_ft_fail)

WITH
-- Every deliverable sample's calls on this chromosome, scored allele by allele. The
-- LEFT JOIN to `filter_set_info` is what makes the coalesce defaults reachable: a missing
-- row means the allele was never scored, which `import_gvs.py:356-363` treats as passing,
-- and Step 4 confirmed Hail's `as_vets.get` misses exactly the same alleles.
calls AS (
  SELECT
    aa.sample_id,
    aa.location,
    aa.call_GQ,
    SAFE_CAST(si.sample_name AS INT64)                                   AS person_id,
    st.vid                                                               AS vid,

    -- import_gvs.py:356 -- missing filter row defaults to YES.
    IFNULL(f.yng_status = "Y", TRUE)                                     AS is_yes,
    -- import_gvs.py:358-363 -- missing filter row defaults to OK. Thresholds confirmed
    -- from the delivered VDS globals, not assumed from the script's defaults.
    IFNULL(f.calibration_sensitivity <=
           IF(LENGTH(aa.ref) = 1 AND LENGTH(aa.allele) = 1, 0.997, 0.990),
           TRUE)                                                         AS is_ok

  FROM `foxtrot.alt_allele` AS aa

  -- Deliverable samples only, and `person_id` derived exactly as the mapping table
  -- derives it (GvsCreateParticipantMappingTable.wdl:60,66-67).
  JOIN `foxtrot.sample_info` AS si
    ON  si.sample_id  = aa.sample_id
    AND si.withdrawn  IS NULL
    AND si.is_control = false

  LEFT JOIN `foxtrot.filter_set_info` AS f
    ON  f.location = aa.location
    AND f.ref      = aa.ref
    AND f.alt      = aa.allele
    AND f.filter_set_name = "foxtrot_v4_2025_07_29"
    AND f.location BETWEEN window_start AND window_end

  -- VAT membership. Tagged here, applied at the very end -- this is defect 3.
  LEFT JOIN `foxtrot.allele_status` AS st
    ON  st.location = aa.location
    AND st.ref      = aa.ref
    AND st.alt      = aa.allele
    AND st.location BETWEEN window_start AND window_end

  WHERE aa.location BETWEEN window_start AND window_end
    -- A sample name that does not cast to INT64 produced no mapping entry, so it can
    -- produce no exclusion. Mirrors the workflow's ARRAY_AGG(... IGNORE NULLS).
    AND SAFE_CAST(si.sample_name AS INT64) IS NOT NULL
),

-- Hail's FT, per genotype, over ALL of the genotype's alleles. The `~any_no` term of
-- `import_gvs.py:381` is dropped as vacuous: `yng_status` is only ever `G` or `Y`. No
-- arity special case is needed -- for a single-alt genotype `any_yes` reduces to that
-- allele's YES and `all_ok` to its OK.
scored AS (
  SELECT
    vid,
    person_id,
    location,
    call_GQ,
    LOGICAL_OR(is_yes) OVER genotype OR LOGICAL_AND(is_ok) OVER genotype  AS ft
  FROM calls
  WINDOW genotype AS (PARTITION BY sample_id, location)
)

SELECT
  vid,
  person_id,
  location,
  call_GQ = 0                                                            AS is_gq0,
  NOT ft                                                                 AS is_ft_fail
FROM scored
WHERE vid IS NOT NULL              -- the VAT restriction, after FT and not before
  AND (call_GQ = 0 OR NOT ft);

-- Not deduplicated on purpose. The only source of a repeated `(vid, person_id)` is the
-- duplicated ingestion of VS-2010 -- 40,000 byte-identical vet rows across 4 samples on
-- chr21, extent elsewhere unmeasured -- since QUERY E established that no
-- `(sample_id, location)` is two genotypes. A `DISTINCT` here would cost a full extra
-- shuffle over ~2.5 billion rows per chromosome to remove them. Step 6 aggregates with
-- `ARRAY_AGG(DISTINCT person_id)` and removes by value, so duplicates are inert there;
-- Section D uses `COUNT(DISTINCT ...)` so they do not inflate the reported figures
-- either. The one place to remember it is any raw `COUNT(*)` on this table.


-- =====================================================================================
-- SECTION C -- driver for the 24 shards
-- =====================================================================================
--
-- Sequential, not parallel. These are large slot-hungry jobs and running them one at a
-- time keeps the reservation predictable and the failure attributable; the whole point of
-- sharding was restartability, not wall clock. Re-running a single chromosome is
-- `for c in 7` -- the DELETE makes it clean.
--
-- **chr21 first, then chr2, then the rest descending.** Not 1..24. chr21 is the only
-- chromosome whose answer is already known two independent ways -- Step 3's 0.0978% and
-- Step 4's 0.0987% on a different slice -- so running it first turns 22 minutes into a
-- calibration: run D1, compare, and only then commit the remaining 20 hours. chr2 second
-- because it is the largest, so if the query survives chr21 but not scale, that surfaces
-- at hour two rather than hour twenty.
--
-- Before any of it, dry-run one chromosome to get the real billed bytes. Free, and it
-- settles the estimate against the 428 TiB the preflight reported for the whole table:
--
--   bq query --dry_run --nouse_legacy_sql --project_id="${PROJECT}" < /tmp/step5_chr21.sql
--
-- Note the DELETE makes the script DML, and a dry run reports bytes for the statements it
-- can plan; if it balks at the multi-statement form, dry-run the INSERT's SELECT alone --
-- the DELETE is partition-metadata work and scans nothing.
--
-- RUN 2026-09-04 against `aou-genomics-curation-prod`: it planned the multi-statement form
-- without complaint and reported 1,558,142,894,853 bytes = 1.417 TiB for the chr21 shard.
-- That is preflight item 2 settled -- ~80.1 TiB and ~$501 genome-wide.
--
--   #!/usr/bin/env bash
--   set -o errexit -o nounset -o pipefail
--
--   for c in 21 2 1 4 3 6 5 7 8 10 11 12 9 13 14 23 16 15 17 18 19 20 22 24; do
--     echo "=== chromosome ${c} -- $(date -u +%FT%TZ)"
--     sed "s/DECLARE chrom        INT64 DEFAULT 1;/DECLARE chrom        INT64 DEFAULT ${c};/" \
--       step5_section_b.sql > "/tmp/step5_chr${c}.sql"
--     bq query --nouse_legacy_sql --project_id="${PROJECT}" \
--       --label=step5:exclusions --label=chrom:"chr${c}" \
--       < "/tmp/step5_chr${c}.sql"
--   done
--
-- Extract Section B into `step5_section_b.sql` for that to work -- the `sed` needs the
-- DECLARE line to be the only match in the file, which it is not while the header
-- comments above are attached.
--
-- The `--label`s are worth the keystrokes: `cost_observability` will not see these, so
-- job labels are the only way to reconstruct what this cost afterwards.


-- =====================================================================================
-- SECTION D -- reconcile, and measure per-VID inflation
-- =====================================================================================
--
-- Reads `exclusions` and the delivered mapping table only. Cheap relative to Section B,
-- but the mapping table is 1.6 billion rows, so not free.

-- D1 -- per-chromosome shape. Run after each shard, or once at the end.
--
-- The check that matters is `rate`, against Step 3's chr21 figure of 0.0978% and Step 4's
-- independent 0.0987% on a different slice. A chromosome landing far off that is a
-- finding, not noise -- and chr21's own row is the direct comparison, since it is the one
-- already measured two other ways.
--
-- `entries` is the denominator the mapping table actually holds for this chromosome, so
-- it is recomputed from `alt_allele` rather than taken from the mapping table, which has
-- already had the withdrawn/control scrub applied and would otherwise mix two
-- populations. This costs a second pass; if that is unwelcome, drop `entries` and `rate`
-- and compare absolute counts against Step 3 scaled by allele count instead.

SELECT
  DIV(e.location, 1000000000000)                        AS chrom,
  COUNT(*)                                              AS exclusion_rows,
  COUNT(DISTINCT FORMAT("%s|%d", e.vid, e.person_id))   AS distinct_entries,
  COUNT(DISTINCT e.vid)                                 AS affected_vids,
  COUNTIF(e.is_gq0)                                     AS gq0_rows,
  COUNTIF(e.is_ft_fail AND NOT e.is_gq0)                AS ft_only_rows
FROM `foxtrot.exclusions` AS e
GROUP BY ROLLUP(chrom)
ORDER BY chrom;


-- D2 -- per-VID relative inflation.
--
-- The headline rate is an average over 1.6 billion VIDs and says nothing about what a
-- researcher actually experiences. Step 3 found the damage concentrated in 0.47% of VIDs
-- at ~359 bad entries each; this is that measurement at genome scale, per VID rather than
-- in aggregate, and it is the artifact to hand to anyone asking "how wrong was my
-- cohort".
--
-- `inflation` is the multiplicative overstatement: 1.0 means the row is correct, 2.0
-- means half its listed carriers are not carriers. A VID whose entries are ALL excluded
-- gives `corrected = 0` and a null inflation -- those are the rows Step 6 empties
-- entirely, counted separately below rather than divided by zero.
--
-- Written to a table because it is worth keeping and because the ORDER BY at the end
-- would otherwise be a single-slot sort over the whole result.

CREATE OR REPLACE TABLE `foxtrot.exclusion_inflation`
CLUSTER BY vid
AS
WITH per_vid AS (
  SELECT
    vid,
    COUNT(DISTINCT person_id) AS excluded
  FROM `foxtrot.exclusions`
  GROUP BY vid
)
SELECT
  m.vid,
  ARRAY_LENGTH(m.person_ids)                                       AS delivered,
  IFNULL(x.excluded, 0)                                            AS excluded,
  ARRAY_LENGTH(m.person_ids) - IFNULL(x.excluded, 0)               AS corrected,
  SAFE_DIVIDE(ARRAY_LENGTH(m.person_ids),
              ARRAY_LENGTH(m.person_ids) - IFNULL(x.excluded, 0))  AS inflation
FROM `foxtrot.vid_to_participant_mapping_2026_08_28` AS m
LEFT JOIN per_vid AS x USING (vid);

-- D3 -- read the distribution off it. Pennies.
--
-- `corrected < 0` must be zero. The correction is purely subtractive -- the new predicate
-- is the old one plus `GQ != 0` plus `FT` -- so a corrected count below zero would mean
-- `exclusions` names a person the delivered row does not contain, which can only happen
-- if the two were built from different populations or different VAT versions. It is the
-- cheapest possible check on the whole exercise and it gates Step 6.
--
-- `orphan_exclusions` is the same question from the other side: exclusion rows whose VID
-- has no row in the delivered table at all. A small nonzero count is expected and benign
-- -- a VID all of whose carriers were withdrawn produced an empty array and may have been
-- dropped -- but a large one means `allele_status` and the delivered table disagree about
-- which VIDs exist, which must be understood before patching anything.

SELECT
  COUNT(*)                                              AS vids,
  COUNTIF(excluded > 0)                                 AS vids_affected,
  SAFE_DIVIDE(COUNTIF(excluded > 0), COUNT(*))          AS affected_rate,
  SUM(delivered)                                        AS entries_delivered,
  SUM(excluded)                                         AS entries_excluded,
  SAFE_DIVIDE(SUM(excluded), SUM(delivered))            AS entry_rate,
  COUNTIF(corrected = 0 AND delivered > 0)              AS vids_emptied,
  COUNTIF(corrected < 0)                                AS vids_negative,
  APPROX_QUANTILES(IF(excluded > 0, inflation, NULL), 100)[OFFSET(50)] AS median_inflation_affected,
  APPROX_QUANTILES(IF(excluded > 0, inflation, NULL), 100)[OFFSET(99)] AS p99_inflation_affected,
  MAX(inflation)                                        AS max_inflation
FROM `foxtrot.exclusion_inflation`;

SELECT COUNT(*) AS orphan_exclusion_vids
FROM (SELECT DISTINCT vid FROM `foxtrot.exclusions`) AS e
LEFT JOIN `foxtrot.vid_to_participant_mapping_2026_08_28` AS m USING (vid)
WHERE m.vid IS NULL;
