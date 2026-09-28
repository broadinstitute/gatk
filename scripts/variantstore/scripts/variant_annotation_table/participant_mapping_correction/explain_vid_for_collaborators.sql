-- =====================================================================================
-- EXPLAIN ONE VID -- the standing answer to "why did the carrier count for X change?"
-- =====================================================================================
--
-- WHAT THIS IS FOR. Collaborators reviewing the corrected Foxtrot participant mapping table
-- are, as far as we can tell, comparing a participant count taken from the SOFT-filtered VDS
-- against one taken from the HARD-filtered VAT, without consulting `FT`. Those two numbers
-- are not supposed to agree, and the gap between them is exactly what the correction removes.
-- This query decomposes that gap for a single VID so the reconciliation can be shown rather
-- than asserted.
--
-- THE IDENTITY IT DEMONSTRATES:
--
--     delivered  =  corrected  +  gq0_only  +  ft_only  +  both
--
-- `delivered` is what the mapping table named before the correction -- the number a
-- soft-filtered, FT-unaware read of the VDS reproduces. `corrected` is what it names now.
-- The three middle terms are the named participants who do not actually hold a passing
-- genotype: GQ 0 no-calls, `FT`-failing calls, and calls that are both. A collaborator
-- looking at the soft-filtered VDS sees `delivered`; the VAT's allele count reflects
-- something close to `corrected`; the difference is not an error in either artifact.
--
-- COST. `exclusions` is PARTITION BY RANGE_BUCKET(location) per chromosome and CLUSTER BY
-- vid, and the partition bound is computed from the VID's own chromosome below, so that read
-- prunes to one partition and a few blocks. `mapping_correction_audit` is ~5.4M rows
-- clustered by vid. The VAT read is the expensive term -- it has no vid clustering we rely
-- on here, so budget a few dollars per invocation, and batch several VIDs into one run
-- (see the BATCH FORM at the bottom) rather than looping this.
--
-- WHAT IT DOES NOT NEED. Not the 18.65 TiB mapping tables, and not the `delivered` /
-- `n_removed` audit columns on `vid_to_participant_mapping_2026_09_05`. Everything here
-- survives STEP 6.8. The audit columns are a convenience, not a dependency: for a changed VID
-- `mapping_correction_audit` carries the same numbers, and for an unchanged VID `removed` is
-- 0 and `delivered` is `ARRAY_LENGTH(person_ids)` by construction.

DECLARE target_vid   STRING DEFAULT '4-76917565-A-AG';   -- <- edit this
DECLARE chrom        INT64;
DECLARE window_start INT64;
DECLARE window_end   INT64;

SET chrom = (
  SELECT CASE SPLIT(target_vid, '-')[OFFSET(0)]
           WHEN 'X' THEN 23
           WHEN 'Y' THEN 24
           ELSE CAST(SPLIT(target_vid, '-')[OFFSET(0)] AS INT64)
         END
);
SET window_start =  chrom      * 1000000000000;
SET window_end   = (chrom + 1) * 1000000000000 - 1;

WITH excl AS (
  -- The reason breakdown. `is_gq0` and `is_ft_fail` are not mutually exclusive -- a call can
  -- be both -- so counting them as three disjoint buckets is what makes the identity close.
  SELECT
    COUNTIF(is_gq0 AND NOT is_ft_fail)  AS gq0_only,
    COUNTIF(is_ft_fail AND NOT is_gq0)  AS ft_only,
    COUNTIF(is_gq0 AND is_ft_fail)      AS both,
    COUNT(*)                            AS excluded_rows,
    COUNT(DISTINCT person_id)           AS excluded_people
  FROM `foxtrot.exclusions`
  WHERE location BETWEEN window_start AND window_end
    AND vid = target_vid
),
aud AS (
  -- Absent for a VID the correction did not touch, which is why every reference below is
  -- wrapped in IFNULL against the unchanged case.
  SELECT delivered, removed, corrected, fully_emptied
  FROM `foxtrot.mapping_correction_audit`
  WHERE vid = target_vid
),
vat AS (
  SELECT ANY_VALUE(gvs_all_ac) AS ac, ANY_VALUE(gvs_all_an) AS an
  FROM `foxtrot.foxtrot_v4_2025_07_29_vat_v9_r2_p2`
  WHERE vid = target_vid
)
SELECT
  target_vid                                                    AS vid,
  (SELECT delivered FROM aud)                                   AS delivered,
  (SELECT corrected FROM aud)                                   AS corrected,
  IFNULL((SELECT removed FROM aud), 0)                          AS removed,
  IFNULL((SELECT fully_emptied FROM aud), FALSE)                AS fully_emptied,
  -- the decomposition
  (SELECT gq0_only FROM excl)                                   AS gq0_only,
  (SELECT ft_only  FROM excl)                                   AS ft_only,
  (SELECT both     FROM excl)                                   AS both,
  -- the identity: 0 means delivered = corrected + gq0_only + ft_only + both.
  -- A small POSITIVE residual is expected where VS-2010 duplicated a vet batch, since the
  -- delivered array can hold a person twice while `exclusions` counts the pair once. It is
  -- bounded by 932 genome-wide across 924 VIDs, so anything above 2 here is not that.
  (SELECT delivered FROM aud)
    - IFNULL((SELECT corrected FROM aud), 0)
    - (SELECT excluded_rows FROM excl)                          AS identity_residual,
  -- the VAT side, i.e. what a hard-filtered read reports
  (SELECT ac FROM vat)                                          AS gvs_all_ac,
  (SELECT an FROM vat)                                          AS gvs_all_an,
  DIV((SELECT an FROM vat), 2)                                  AS participants_called;

-- HOW TO READ IT.
--
-- `identity_residual` = 0 is the headline: every participant dropped from this VID's array is
-- accounted for by a named reason, and none was dropped for any other cause.
--
-- `gvs_all_ac` against `corrected` is the cross-artifact check, and it is APPROXIMATE by
-- construction -- do not present it as an equality. AC counts ALLELES, so a homozygous
-- carrier contributes 2 while the mapping table names them once; AC is computed on the VDS
-- while the mapping table derives from `alt_allele`; and the VAT is hard filtered at the site
-- level by the max-over-class rule, which is not the same predicate as the per-genotype `FT`
-- this correction applies. Expect the same order of magnitude, not the same number. What
-- matters is that `corrected` is close to `gvs_all_ac` while `delivered` is not -- often by
-- four or five orders of magnitude -- which is the whole demonstration.
--
-- `fully_emptied = TRUE` means this is one of the 172 case-4 VIDs, and the reconciliation
-- WILL NOT close against the VAT. Every named carrier was excluded while `gvs_all_ac` is
-- nonzero, because the passing carriers sit at a different allele representation and were
-- never in the delivered array at all. For those, `excluded + participants_called` can even
-- exceed the 535,662 active non-control cohort -- 3 of the 172 do -- which is not a
-- contradiction, because the two counts are measured at different sites. Do not walk a
-- collaborator through one of these as an introductory example; they are the hardest case.
-- Use a `fully_emptied = FALSE` VID with a large `ft_only` first.
--
-- A NULL `delivered` means the correction did not touch this VID: nothing was over-reported
-- there, `removed` is 0, and its array is unchanged from what shipped.

-- =====================================================================================
-- BATCH FORM -- same decomposition for a list of VIDs, one VAT read.
-- =====================================================================================
--
-- Prefer this when preparing examples ahead of a meeting. The single-VID form above is for
-- answering a question live. Partition pruning is lost here unless every VID is on one
-- chromosome, so this scans `exclusions` -- still cheap next to the VAT read.

-- DECLARE target_vids ARRAY<STRING> DEFAULT [
--   '4-76917565-A-AG',
--   '12-124467987-A-AC'
-- ];
--
-- WITH excl AS (
--   SELECT vid,
--          COUNTIF(is_gq0 AND NOT is_ft_fail) AS gq0_only,
--          COUNTIF(is_ft_fail AND NOT is_gq0) AS ft_only,
--          COUNTIF(is_gq0 AND is_ft_fail)     AS both,
--          COUNT(*)                           AS excluded_rows
--   FROM `foxtrot.exclusions`
--   WHERE vid IN UNNEST(target_vids)
--   GROUP BY vid
-- ),
-- vat AS (
--   SELECT vid, ANY_VALUE(gvs_all_ac) AS ac, ANY_VALUE(gvs_all_an) AS an
--   FROM `foxtrot.foxtrot_v4_2025_07_29_vat_v9_r2_p2`
--   WHERE vid IN UNNEST(target_vids)
--   GROUP BY vid
-- )
-- SELECT t AS vid, a.delivered, a.corrected, a.fully_emptied,
--        e.gq0_only, e.ft_only, e.both,
--        a.delivered - IFNULL(a.corrected, 0) - e.excluded_rows AS identity_residual,
--        v.ac AS gvs_all_ac, v.an AS gvs_all_an
-- FROM UNNEST(target_vids) AS t
-- LEFT JOIN `foxtrot.mapping_correction_audit` AS a ON a.vid = t
-- LEFT JOIN excl AS e ON e.vid = t
-- LEFT JOIN vat  AS v ON v.vid = t
-- ORDER BY a.delivered DESC;
