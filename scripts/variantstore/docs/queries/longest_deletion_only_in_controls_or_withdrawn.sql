-- =============================================================================
-- Query: Loci where the longest deletion is carried exclusively by
--        control or withdrawn samples (superpartition 104, chromosome 20)
-- =============================================================================
--
-- Background
-- ----------
-- In GVS the `vet_NNN` tables store one row per variant per sample. At a
-- given `location` (chromosome × 1,000,000,000,000 + position), different
-- samples may have different REF allele lengths because they carry deletions
-- of different depths.  A "deletion" is any row where LENGTH(ref) > LENGTH(alt).
--
-- This query finds locations on chromosome 20 (within superpartition 104)
-- where:
--   1.  At least one sample carries a deletion.
--   2.  The deepest deletion at that locus (longest REF allele among all
--       deletion carriers) is carried ONLY by control or withdrawn samples.
--   3.  No non-control, non-withdrawn sample carries a deletion of equal
--       depth (i.e. equal REF length) at the same locus.
--
-- The result confirms the condition described in the user's statement:
--   "the reference allele for these [control/withdrawn] samples would be
--    longer than the reference allele for any non-control, non-withdrawn
--    samples at the locus."
--
-- Schema references
-- -----------------
-- sample_info  : sample_id (INT), sample_name (STRING), is_control (BOOL),
--                withdrawn (TIMESTAMP – NULL means not withdrawn)
-- vet_NNN      : sample_id (INT), location (INT), ref (STRING), alt (STRING),
--                call_GT (STRING), …
--                partitioned by sample_id, clustered by location
--
-- Location encoding (SchemaUtils.encodeLocation)
-- -----------------------------------------------
--   location = chromosome_index * 1_000_000_000_000 + genomic_position
--   chromosome 20 → index 20
--   chr20 location range: [20_000_000_000_000, 21_000_000_000_000)
--
-- Superpartition 104 sample_id range
-- ------------------------------------
--   PARTITION_START = (104 - 1) * 4000 + 1 = 412_001
--   PARTITION_END   = 412_001 + 4000       = 416_001  (exclusive upper bound)
--   Usable sample_ids: 412_001 – 416_000
--
-- Multi-allelic note
-- ------------------
-- Rows with comma-separated ALT values (e.g. 'A,ACG' for a 1/2 hetvar) are
-- excluded from the deletion check because LENGTH(alt) would include the
-- comma character, yielding an incorrect REF-vs-ALT length comparison.
-- If multi-allelic deletions are important, a separate UNNEST/SPLIT pass
-- over the alt field would be required.
--
-- Usage
-- -----
-- Replace <project> and <dataset> with the actual GCP project ID and
-- BigQuery dataset name before running.
-- =============================================================================

WITH

-- -------------------------------------------------------------------------
-- Step 1 – Collect simple deletion rows for chr20, superpartition 104
-- -------------------------------------------------------------------------
-- Filters applied here:
--   • location range  → only chromosome 20
--   • sample_id range → only superpartition 104 (enables partition pruning)
--   • LENGTH(ref) > LENGTH(alt) → only deletion alleles
--   • NOT CONTAINS_SUBSTR(alt, ',') → exclude multi-allelic rows (see note)
--   • call_GT excludes homozygous-reference calls (should not normally exist
--     in vet, but guarded here for safety)
chr20_deletions AS (
    SELECT
        v.sample_id,
        v.location,
        v.ref,
        v.alt,
        LENGTH(v.ref)                       AS ref_length,
        LENGTH(v.alt)                       AS alt_length,
        LENGTH(v.ref) - LENGTH(v.alt)       AS deleted_bases,
        si.is_control,
        si.withdrawn IS NOT NULL            AS is_withdrawn,
        (si.is_control OR si.withdrawn IS NOT NULL)
                                            AS is_control_or_withdrawn
    FROM
        `<project>.<dataset>.vet_104`   AS v
        INNER JOIN `<project>.<dataset>.sample_info` AS si
            ON v.sample_id = si.sample_id
    WHERE
        -- ── Chromosome 20 location range ─────────────────────────────────
        v.location >= 20000000000000       -- 20 * 1,000,000,000,000
        AND v.location  < 21000000000000   -- 21 * 1,000,000,000,000

        -- ── Superpartition 104 sample_id range (partition pruning) ────────
        AND v.sample_id >= 412001
        AND v.sample_id  < 416001

        -- ── Simple deletions only (REF longer than ALT) ───────────────────
        AND LENGTH(v.ref) > LENGTH(v.alt)
        AND NOT CONTAINS_SUBSTR(v.alt, ',')   -- exclude multi-allelic rows

        -- ── Exclude ref-only calls (safety guard) ─────────────────────────
        AND v.call_GT NOT IN ('0/0', '0|0')
),

-- -------------------------------------------------------------------------
-- Step 2 – For each locus, find the maximum REF length among all deletions
-- -------------------------------------------------------------------------
max_ref_per_location AS (
    SELECT
        location,
        MAX(ref_length) AS max_ref_length
    FROM
        chr20_deletions
    GROUP BY
        location
),

-- -------------------------------------------------------------------------
-- Step 3 – Identify which samples carry the longest deletion at each locus
-- -------------------------------------------------------------------------
carriers_of_max_deletion AS (
    SELECT
        d.location,
        d.sample_id,
        d.ref,
        d.alt,
        d.deleted_bases,
        d.is_control,
        d.is_withdrawn,
        d.is_control_or_withdrawn
    FROM
        chr20_deletions                AS d
        INNER JOIN max_ref_per_location AS m
            ON  d.location  = m.location
            AND d.ref_length = m.max_ref_length
)

-- -------------------------------------------------------------------------
-- Final – Keep only loci where every carrier of the deepest deletion is
--         a control or withdrawn sample
-- -------------------------------------------------------------------------
SELECT
    -- Decoded genomic position (chromosome 20 assumed)
    MOD(location, 1000000000000)            AS position,
    location,

    -- Deletion depth (in bases deleted)
    MAX(deleted_bases)                      AS max_deleted_bases,

    -- Representative REF / ALT alleles from the deepest deletion
    ANY_VALUE(ref)                          AS example_ref,
    ANY_VALUE(alt)                          AS example_alt,

    -- Carrier counts
    COUNT(*)                                AS total_max_deletion_carriers,
    COUNTIF(is_control_or_withdrawn)        AS control_or_withdrawn_carrier_count,
    COUNTIF(NOT is_control_or_withdrawn)    AS non_control_withdrawn_carrier_count,  -- always 0 here

    -- Comma-separated carrier sample IDs with labels
    STRING_AGG(
        CONCAT(
            CAST(sample_id AS STRING),
            CASE
                WHEN is_control   THEN ' [ctrl]'
                WHEN is_withdrawn THEN ' [withdrawn]'
                ELSE ''
            END
        ),
        ', ' ORDER BY sample_id
    )                                       AS max_deletion_carrier_sample_ids

FROM
    carriers_of_max_deletion
GROUP BY
    location

-- The key predicate: no non-control, non-withdrawn sample carries the
-- deepest deletion at this locus
HAVING
    COUNTIF(NOT is_control_or_withdrawn) = 0

ORDER BY
    position
;

