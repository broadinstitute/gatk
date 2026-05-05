WITH valid_samples AS (
    -- 1. Get all non-withdrawn, non-control samples.
    -- SAFE_CAST prevents query failures if any alphanumeric names bypass the is_control filter.
    SELECT SAFE_CAST(sample_name AS INT64) AS person_id
    FROM `foxtrot.sample_info`
    WHERE withdrawn IS NULL
      AND is_control = FALSE
),


     mapping_unnested AS (
         -- 2. Flatten the mapping table and filter down to only valid samples.
         SELECT DISTINCT
             p_id AS person_id,
             m.vid
         FROM `foxtrot.vid_to_participant_mapping_2025_10_28` m
                  CROSS JOIN UNNEST(m.person_ids) AS p_id
                  INNER JOIN valid_samples vs
                             ON p_id = vs.person_id
                                 -- Ensure the SAFE_CAST didn't result in a NULL person_id
                                 AND vs.person_id IS NOT NULL
     )


-- 3. Find valid sample/vid pairs where the VID is in Classic but NOT in r2.
SELECT
    m.person_id,
    m.vid
FROM mapping_unnested m
WHERE EXISTS (
    -- Confirm this VID appears in the classic VAT
    SELECT 1
    FROM `foxtrot.foxtrot_v4_2025_07_29_vat` classic_vat
    WHERE classic_vat.vid = m.vid
)
  AND NOT EXISTS (
    -- Confirm this VID does NOT appear in the r2 VAT
    SELECT 1
    FROM `foxtrot.foxtrot_v4_2025_07_29_vat_v9_r2` r2_vat
    WHERE r2_vat.vid = m.vid
);
